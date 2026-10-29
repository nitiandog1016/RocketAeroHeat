"""将逐态重建的真实连续试算接回原面工作台，严格绑定方向、前驱与累计历史。"""
from pathlib import Path
import json, os, shutil, sys, re
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT, sha

OUT = ROOT / 'outputs/v87_validation'
STUDY = OUT / 'continuous01/active_continuation.json'


def configure(output_folder):
    """显式选择本轮证据目录；缺省仍重现v87，不改写历史目录。"""
    global OUT, STUDY
    folder = Path(output_folder)
    OUT = (folder if folder.is_absolute() else ROOT / folder).resolve()
    OUT.relative_to(ROOT.resolve())
    assert OUT.is_dir()
    STUDY = OUT / 'continuous01/active_continuation.json'


def stamp_viewer_urls(folder, data):
    """按实际文件内容给内嵌页面和脚本加版本，防止父页面刷新后仍沿用旧嵌套资源。"""
    for case in data['cases']:
        relative = case['resultsURL'].split('?', 1)[0]
        file = folder/relative
        file.resolve().relative_to(folder.resolve())
        assert file.is_file()
        case['resultsURL'] = relative+'?view='+sha(file)[:12]
    raw = json.dumps(data, ensure_ascii=False, separators=(',', ':'))
    (folder/'wall-velocity-study-data.json').write_bytes(raw.encode('utf8'))
    script = folder/'wall-velocity-study-data.js'
    script.write_bytes(('window.ROCKET_WALL_VELOCITY_STUDY='+raw+';\n').encode('utf8'))
    page = folder/'index.html'; html = page.read_text('utf8')
    for name in ('wall-velocity-study-data.js', 'rocket-wall-velocity.js'):
        html, count = re.subn(r'src="'+re.escape(name)+r'(?:\?[^"]*)?"',
                             'src="'+name+'?v='+sha(folder/name)[:12]+'"', html)
        assert count == 1
    name = 'rocket-sensitivity.css'
    html, count = re.subn(r'href="'+re.escape(name)+r'(?:\?[^\"]*)?"',
                         'href="'+name+'?v='+sha(folder/name)[:12]+'"', html)
    assert count == 1
    page.write_bytes(html.encode('utf8'))


def aliases():
    """仅建立已保存端点的逐字节别名；不执行计算，不增加接受次数。"""
    final = STUDY.is_file()
    data = json.loads((STUDY if final else OUT / 'continuous01/progress.json').read_text('utf8'))
    steps = data['steps'] if final else data
    if isinstance(steps, dict):
        steps = [steps]
    entries = []
    previous = Path(steps[0]['sourceBinding']['sourceFile']).parent.name
    match = re.fullmatch(r'krylov(\d+)', previous)
    assert match, '显示编号必须绑定本轮实际首个前驱，不能沿用旧版本编号。'
    initial_number = int(match.group(1))
    for step in steps:
        if not step['main'] or step['main']['acceptedNonlinearIterations'] == 0:
            continue
        source = Path(step['resultFolder'])
        target = ROOT / 'outputs/v83_validation' / f"krylov{initial_number + step['step']}"
        target.mkdir(exist_ok=True)
        for name in ('wall_velocity_krylov_evidence.mat', 'wall_velocity_krylov_summary.json'):
            if not (target / name).exists():
                os.link(source / name, target / name)
            assert sha(source / name) == sha(target / name)
        entries.append(dict(step=step['step'], sourceFolder=str(source), displayCase=target.name,
                            predecessorCase=previous, SHA256=sha(target / 'wall_velocity_krylov_evidence.mat'),
                            RMS=step['main']['rmsMerit'], maximum=step['main']['maximumResidual']))
        previous = target.name
    # 没有接受更新时显示原检查点，诊断仍可接回工作台；不伪造新态或新接受编号。
    selected = entries[-1] if entries else dict(step=steps[-1]['step'], sourceFolder=str(Path(steps[0]['sourceBinding']['sourceFile']).parent),
        displayCase=Path(steps[0]['sourceBinding']['sourceFile']).parent.name, predecessorCase=previous,
        SHA256=sha(steps[0]['sourceBinding']['sourceFile']), displayedStateUnchanged=True)
    result = dict(status='byte-identical-actual-continuous-current-state-display-aliases',
                  studyTerminal=final, entries=entries, selected=selected,
                  newFlowAcceptances=0, newProductionSolvedAreaM2=0)
    if final:
        target = OUT / 'view_source_binding.json'
        assert not target.exists()
        target.write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
    print(json.dumps(result))


def bind(study_file, candidate_file, state_file, step_number):
    """验证云图所示端点确实使用该当前状态生成的完整方向。"""
    study_file, candidate_file, state_file = map(Path, (study_file, candidate_file, state_file))
    study = json.loads(study_file.read_text('utf8'))
    steps = study['steps']
    if isinstance(steps, dict):
        steps = [steps]
    step = next(item for item in steps if item['step'] == step_number)
    accepted = bool(step['main'] and step['main']['acceptedNonlinearIterations'] == 1)
    assert sha(candidate_file) == step['directionFileSHA256'] == sha(step['directionFile'])
    candidate = loadmat(candidate_file, simplify_cells=True,
                        variable_names=['value', 'directionReport', 'sourceBinding', 'evidence',
                                        'preparationFile', 'preparationSHA256'])
    state = loadmat(state_file, simplify_cells=True, variable_names=['result', 'report'])
    actual = state['result']
    assert candidate['sourceBinding'] == step['sourceBinding']
    assert sha(candidate['preparationFile']) == candidate['preparationSHA256'] == step['preparationSHA256']
    assert sha(candidate['sourceBinding']['sourceFile']) == candidate['sourceBinding']['sourceFileSHA256']
    if accepted:
        line = actual['linearLog']; line = line if isinstance(line, dict) else line[0]
        np.testing.assert_array_equal(candidate['value'], line['steadyDirectionEvidence']['directionU'])
        probe = loadmat(actual['options']['steadyInitialGuess']['file'], simplify_cells=True,
                        variable_names=['sourceBinding', 'directionU'])
        assert probe['sourceBinding']['sha256'] == candidate['sourceBinding']['sourceFileSHA256']
        np.testing.assert_array_equal(probe['directionU'], candidate['value'])
        assert step['main']['totalNewEquationAcceptedIterations'] == state['report']['totalNewEquationAcceptedIterations']
        assert step['main']['totalNewEquationResidualEvaluations'] == state['report']['totalNewEquationResidualEvaluations']
        assert sha(state_file) == sha(Path(step['resultFolder']) / 'wall_velocity_krylov_evidence.mat')
    else:
        assert sha(state_file) == candidate['sourceBinding']['sourceFileSHA256'], '拒绝试算只能展示实际当前态。'
    report = candidate['directionReport']
    coordinates = report.get('repreview', {}).get('trialStateCoordinates',
                       report.get('lookaheadSelection', {}).get('trialStateCoordinates', 'conservative'))
    assert coordinates in ('conservative', 'primitive-tangent')
    assert report.get('lookaheadSelection', {}).get('trialStateCoordinates', 'conservative') == coordinates
    if accepted:
        assert actual['options'].get('trialStateCoordinates', 'conservative') == coordinates, '前瞻路径与实际接受路径不符。'
    assert (report['bandDirectionFeasible'] or not accepted) and not report['oldActiveSetInherited']
    detail = report['linearReport']
    cycles = candidate['evidence']['cycles']
    cycles = [cycles] if isinstance(cycles, dict) else list(np.atleast_1d(cycles).ravel())
    projections = []
    for cycle in cycles:
        p = cycle['projector']
        q = int(p.get('independentConstraintCount', np.size(p['normalizedQRDiagonal'])))
        representation = dict(originalRows=int(p['constraintCount']), distinctConditions=q,
            compactStorage='basisCoordinateIndices' in p, storedCoordinates=int(np.size(p.get('basisCoordinateIndices', []))),
            basisArrayBytes=int(p.get('basisArrayBytes', np.asarray(p['basis']).nbytes)),
            equivalentFullBasisArrayBytes=int(p.get('equivalentFullBasisArrayBytes', np.asarray(p['basis']).nbytes)),
            allOriginalRowsExactlyRepresented=bool(p.get('allOriginalRowsExactlyRepresented',p['constraintCount']==q)),
            rankTruncationUsed=bool(p.get('rankTruncationUsed',False)))
        if 'constraintRowPreconditioning' in p:
            exact = p['constraintRowPreconditioning']
            assert exact['allOriginalConditionsRetained'] and not exact['rankBasedRowSelection'] and not exact['rankTruncationUsed']
            representation['rowPreconditioning'] = dict(method=exact['kind'],
                transformedRows=int(np.size(exact['transformedRowIndices'])),
                allOriginalConditionsRetained=True, rankTruncationUsed=False,
                exactBinaryInverseVerified=bool(exact['binaryInverseReconstructionVerified']))
        projections.append(representation)
    target = detail.get('internalCTarget', detail.get('options', {}).get('relativeTolerance'))
    assert target == .001 and (not accepted or step['probe']['originalNSDirectionTargetMet'])
    preparation = dict(kind='current-band-active-temperature-v1',
                methodLabel=f"当前温域双侧筛选；自由方向及{report['activeSolveCount']}次活跃行投影重算",
                studyFile=str(study_file), studyFileSHA256=sha(study_file),
                candidateFile=str(candidate_file), candidateFileSHA256=sha(candidate_file),
                step=step_number, baseCase=Path(candidate['sourceBinding']['sourceFile']).parent.name,
                usedTotalDirectionExactlyVerified=accepted, proposedTotalDirectionBoundToCurrentState=True,
                displayedStateUnchanged=not accepted, hardUpdateLimitClaimed=False,
                firstSegmentIterations=report['freeLinearReport']['nativeIterations'],
                restartedSegmentIterations=report['totalNativeIterations'] - report['freeLinearReport']['nativeIterations'],
                internalCRelativeResidual=detail['relativeCResidual'], internalCTarget=target,
                internalCTargetMet=bool(detail['internalCTargetMet']),
                completeModelRelativeDefect=step['completeSelectedModelRelativeDefect'],
                originalNSDirectionDefect=line['fullSystemRelativeDefect'] if accepted else
                    (step['probe']['originalNSDirectionDefect'] if step['probe'] else None),
                originalNSDirectionTarget=actual['options']['linearRelativeTolerance'] if accepted else .1,
                bandFaceIDs=step['directionReport']['bandFaceIDs'],
                activeFaceIDs=step['directionReport']['activeFaceIDs'],
                temperatureDerivativeK=step['directionReport']['temperatureDerivativeK'],
                oldActiveSetInherited=False, physicalTemperatureBoundaryAdded=False,
                nonlinearTemperatureBoundGuaranteed=False, constrainedOptimumClaimed=False,
                projectionRepresentation=projections, activeSolveCount=int(report['activeSolveCount']))
    # 复用已绑定JSON记录中的数组列表，避免把MATLAB ndarray直接写入离线JSON。
    selection = step['directionReport'].get('lookaheadSelection')
    if selection:
        # 前瞻集合可能包含原10K之外的面，不能再将全部选择数写成10K带内面数。
        preparation['lookaheadSelection'] = selection
        preparation['methodLabel'] = f"当前自由候选前瞻；保留原温域带并补全部预测越温单元；{report['activeSolveCount']}次受限重算"
        if selection.get('selectionPolicy') == 'original-and-predicted-cell-and-interface-faces':
            preparation['methodLabel'] = f"当前单元与Riemann界面共同前瞻；保留原温域带、单元所需面与可检查的越温界面；{report['activeSolveCount']}次受限重算"
    if 'repreview' in report:
        preparation['repreview'] = step['directionReport']['repreview']
        preparation['repreviewStopReason'] = report['stopReason']
        preparation['bandDirectionFeasible'] = bool(report['bandDirectionFeasible'])
        preparation['methodLabel'] = f"投影后逐轮检查指定步长的内部温域；累计保留必要条件；{report['activeSolveCount']}次受限重算"
    if 'freeCheckpointReplay' in report:
        preparation['freeCheckpointReplay'] = step['directionReport']['freeCheckpointReplay']
        preparation['newFreeNativeIterations'] = report['freeCheckpointReplay']['newFreeNativeIterations']
        preparation['reusedFreeNativeIterations'] = report['reusedFreeNativeIterations']
        preparation['newProjectedNativeIterations'] = report['newNativeIterations']
    if coordinates != 'conservative':
        preparation['proposedTrialStateCoordinates'] = coordinates
        preparation['displayedStateTrialStateCoordinates'] = actual['options'].get('trialStateCoordinates', 'conservative')
        preparation['methodLabel'] += '；采用原始变量切向候选路径'
    if not accepted:
        preparation['unacceptedAttemptNSCounts'] = dict(preparation=step['newPreparationActualNS'],
            probe=step['newProbeActualNS'], main=step['newMainActualNS'])
    return preparation


def build():
    binding = json.loads((OUT / 'view_source_binding.json').read_text('utf8'))
    selected = binding['selected']
    study = json.loads(STUDY.read_text('utf8'))
    steps = study['steps'] if isinstance(study['steps'], list) else [study['steps']]
    step = next(x for x in steps if x['step'] == selected['step'])
    preparation = bind(STUDY, step['directionFile'], ROOT / 'outputs/v83_validation' /
                       selected['displayCase'] / 'wall_velocity_krylov_evidence.mat', selected['step'])
    source, folder = OUT / 'workbench_base', OUT / 'workbench01'
    assert not folder.exists()
    shutil.copytree(source, folder)
    shutil.copy2(ROOT / 'web/rocket-wall-velocity.js', folder / 'rocket-wall-velocity.js')
    data = json.loads((folder / 'wall-velocity-study-data.json').read_text('utf8'))
    version = OUT.name.split('_')[0] if re.fullmatch(r'v\d+_validation', OUT.name) else '当前研究'
    data['cases'][2]['label'] = (f"{version} 第{selected['step']}步方向未接受，保留当前态（未收敛）" if selected.get('displayedStateUnchanged') else
                               f"{version} 当前态重建：第{selected['step']}步实际接受态（未收敛）")
    data['cases'][2]['solverDiagnostics']['directionPreparation'] = preparation
    labels = []
    for item in steps:
        m = item['main']
        labels.append(f"第{item['step']}步：" +
                      (f"接受{m['acceptedNonlinearIterations']}次，RMS {m['rmsMerit']:.3f} /s" if m else '未通过方向检查'))
    data['temperatureActiveResearch'] = dict(selectedCase=selected['displayCase'],
        comparisonLabel='；'.join(labels),
        scopeLabel='逐态重建当前方程、方向和温域行；未收敛研究，新增合格面积0',
        constraintLabel='保留原温域行；指定有限步越域所需面显式激活；每轮检查不能代替开放边界、完整NS或稳态验收',
        reportURL='temperature-active-diagnostics.html', physicalSolution=False,
        newProductionSolvedAreaM2=0, continuation=study)
    raw = json.dumps(data, ensure_ascii=False, separators=(',', ':'))
    (folder / 'wall-velocity-study-data.json').write_text(raw, encoding='utf8')
    # 固定离线脚本字节，不让Windows文本换行转换破坏JSON/JS一致性。
    (folder / 'wall-velocity-study-data.js').write_bytes(('window.ROCKET_WALL_VELOCITY_STUDY=' + raw + ';\n').encode('utf8'))
    stamp_viewer_urls(folder, data)
    # 诊断页逐步列出真实原行/不同条件及表示大小，不把存储节省写成物理或加速收益。
    representation_rows = []
    for item in steps:
        saved = loadmat(item['directionFile'], simplify_cells=True, variable_names=['evidence'])['evidence']
        cycles = saved['cycles']; cycles = [cycles] if isinstance(cycles, dict) else list(np.atleast_1d(cycles).ravel())
        if cycles:
            p = cycles[-1]['projector']; compact = 'basisCoordinateIndices' in p
            stored = int(p.get('basisArrayBytes', np.asarray(p['basis']).nbytes))
            equivalent = int(p.get('equivalentFullBasisArrayBytes', np.asarray(p['basis']).nbytes))
            cells = [str(item['step']), str(np.size(item['directionReport']['bandFaceIDs'])), str(p['constraintCount']),
                     str(p.get('independentConstraintCount',np.size(p['normalizedQRDiagonal']))), '完整局部支撑' if compact else '原完整坐标',
                     f'{stored/1024:.3f}', f'{equivalent/1024:.3f}']
        else:
            cells = [str(item['step']), str(np.size(item['directionReport']['bandFaceIDs'])), '0', '0', '自由方向', '0', '0']
        representation_rows.append('<tr>'+''.join('<td>'+x+'</td>' for x in cells)+'</tr>')
    representation_table = ('<h2>方向约束与表示核对</h2><p>所选数值约束中的每个原面均保留；只有逐项相同的导数行共享表示。'
        '存储大小仅指投影基，不是进程峰值，也不是计算加速或新增物理面积。</p><table><thead><tr>'
        '<th>步骤</th><th>所选原面</th><th>激活原行</th><th>不同条件</th><th>表示</th>'
        '<th>基存储 KiB</th><th>等价完整基 KiB</th></tr></thead><tbody>' + ''.join(representation_rows)+'</tbody></table>')
    conditioned = [p['rowPreconditioning'] for p in preparation['projectionRepresentation'] if 'rowPreconditioning' in p]
    if conditioned:
        representation_table += ('<p>全部不同条件均保留；同支撑近重复行采用精确差值与二进制幂缩放，'
            '逐轮变换行数：'+ '、'.join(str(p['transformedRows']) for p in conditioned)+
            '。未按秩删行；仅证明保存浮点矩阵的可逆表示，不证明微小差异具有独立物理意义。</p>')
    lookahead_rows = []
    interface_policy = False
    for item in steps:
        selection = item['directionReport'].get('lookaheadSelection')
        if not selection:
            continue
        size = lambda x: int(np.size(x))
        preview = selection.get('interfacePreview')
        preview = preview if isinstance(preview, dict) else None
        interface_policy = interface_policy or preview is not None
        cells = [str(item['step']), str(size(selection['originalBandFaceIDs'])),
                 str(size(selection['predictedViolatingCellIDs'])), str(size(selection['requiredOriginalStateFaceIDs'])),
                 str(size(selection.get('predictedViolatingInterfaceFaceIDs', preview['predictedViolatingInterfaceFaceIDs'] if preview else []))) if preview is not None else '未启用',
                 str(size(preview['blockedFaceIDs'])) if preview is not None else '未启用',
                 str(size(selection['selectedBandFaceIDs'])), str(size(selection.get('candidateExtendedBandFaceIDs', selection['selectedBandFaceIDs'])))]
        lookahead_rows.append('<tr>'+''.join('<td>'+x+'</td>' for x in cells)+'</tr>')
    if lookahead_rows:
        policy_text = ('选择集合为原10 K温域带、全部预测越温单元所需原态面及可检查的预测越温Riemann界面。'
                       '候选单元温度超域导致无法检查的相邻面单列为阻挡，不计为已检查。' if interface_policy else
                       '选择集合为原10 K温域带加全部预测越温单元所需原态面。')
        representation_table += ('<h2>当前候选前瞻核对</h2><p>'+policy_text+
            '候选扩展带仅用于来源和距离诊断；选面不按秩删行，也不保证非线性候选合法。'
            '所有候选仍须原完整NS、80–2000 K温域和下降检查。</p><table><thead><tr>'
            '<th>步骤</th><th>原10 K面</th><th>预测越温单元</th><th>所需原态面</th><th>预测越温界面</th><th>阻挡面</th><th>实际选择面</th><th>候选扩展带面</th>'
            '</tr></thead><tbody>'+''.join(lookahead_rows)+'</tbody></table>')
    repreview_rows = []
    for item in steps:
        r = item['directionReport'].get('repreview')
        if not r:
            continue
        history = r['history']; history = [history] if isinstance(history, dict) else history
        for h in history:
            counts = lambda key: str(np.size(h[key]))
            values = [str(item['step']), str(h['round']), counts('predictedViolatingCellIDs'),
                      counts('predictedViolatingInterfaceFaceIDs'), str(h['assessedInterfaceCount']), str(h['blockedInterfaceCount']),
                      counts('selectedFaceIDs'), counts('forcedActiveFaceIDs'),
                      '通过（仅内部温域）' if h['targetCandidateThermallyFeasible'] else ('未通过' if h['previewCellStateEvaluable'] else '未执行完整检查')]
            repreview_rows.append('<tr>'+''.join('<td>'+v+'</td>' for v in values)+'</tr>')
    if repreview_rows:
        representation_table += ('<h2>投影后逐轮温域复查</h2><p>第0轮是自由候选，随后逐轮重新检查投影方向。'
            '累计必要条件包含先前暴露的全部面，旧约束及旧活跃面不删除；有限步越界面即使一阶朝内也需显式激活。'
            '通过仅表示指定步长的单元与内部Riemann温域；开放边界、原128eps方向精度、完整NS和稳态另行检查。</p>'
            '<table><thead><tr><th>步骤</th><th>预览轮</th><th>越温单元</th><th>预测越温界面</th><th>已检查面</th><th>阻挡面</th>'
            '<th>累计选择面</th><th>累计必要面</th><th>内部温域</th></tr></thead><tbody>'+''.join(repreview_rows)+'</tbody></table>'
            '<img src="temperature_repreview_rounds.png" alt="逐轮内部温域检查及累计必要条件">')
    page = '''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>当前温域与连续推进</title>
<style>body{font:16px/1.7 system-ui,sans-serif;max-width:1350px;margin:28px auto;padding:0 24px;color:#213747;background:#f7f9fb}img{width:100%;background:white}a{color:#146da2}.note{padding:16px;border-left:4px solid #ae7c24;background:#fff3d9}table{border-collapse:collapse;width:100%;background:#fff;font-size:14px}th,td{border-bottom:1px solid #d8e0e7;padding:8px;text-align:right}th{color:#324b60}</style>
<a href="index.html" target="_top">返回当前结果工作台</a><h1>真实20°箭体：当前温域与连续推进</h1>
<p class="note">压力、热流与邻壁气体温度为未收敛研究值。当前试算尚未新增合格求解面积。</p>
<p>每步重新建立当前方程与界面温度导数，自由方向先行，必要时再投影重算。给定壁温仍为300 K；采样温度方向限制不是物理边界条件，也不保证非线性候选在物性范围内。</p>
<img src="temperature_active_continuation.png" alt="各步真实残差、原NS缺陷、接受步长和温域方向状态">
<p><a href="temperature_active_continuation.pdf">科研图 PDF</a> · <a href="temperature_active_record.json">连续推进原始记录</a></p></html>'''
    page = page.replace('<img src="temperature_active_continuation.png"', representation_table+'<img src="temperature_active_continuation.png"')
    (folder / 'temperature-active-diagnostics.html').write_text(page, encoding='utf8')
    (folder / 'temperature_active_record.json').write_text(json.dumps(study, ensure_ascii=False, indent=2), encoding='utf8')
    for ext in ('png', 'pdf'):
        shutil.copy2(OUT / f'temperature_active_continuation.{ext}', folder / f'temperature_active_continuation.{ext}')
        if repreview_rows:
            shutil.copy2(OUT/f'temperature_repreview_rounds.{ext}', folder/f'temperature_repreview_rounds.{ext}')
    same = []
    for file in (source / 'cases').glob('*/results.json'):
        relative = file.relative_to(source)
        assert file.read_bytes() == (folder / relative).read_bytes()
        same.append(str(relative))
    receipt = dict(status='bound-current-state-active-temperature-continuation-and-original-facets',
                   selected=selected, directionPreparation=preparation,
                   originalFacetDatasetsByteIdentical=same, physicalSolution=False,
                   newProductionSolvedAreaM2=0, browserVisualOrWebGLRuntimeVerified=False)
    (folder / 'temperature_active_view_binding.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    # 全部面号已写入回执；控制台只输出摘要，避免849336面前瞻生成数MB重复日志。
    print(json.dumps(dict(status=receipt['status'], selected=selected,
                         originalFacetDatasetsByteIdentical=same, physicalSolution=False,
                         newProductionSolvedAreaM2=0)))


if __name__ == '__main__':
    if len(sys.argv) > 2:
        configure(sys.argv[2])
    {'aliases': aliases, 'build': build}[sys.argv[1]]()
