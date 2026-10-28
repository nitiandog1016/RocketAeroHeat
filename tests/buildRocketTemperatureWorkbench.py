"""将同态自由/温度投影研究与实际选用方向接回既有五步骤科研工作台。"""
from pathlib import Path
import json,os,shutil,sys
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,sha
OUT=ROOT/'outputs/v86_validation';COMP=OUT/'comparison01'


def aliases():
    report=json.loads((COMP/'temperature_projection_comparison.json').read_text('utf8'));available=[];entries=[]
    for branch,number in zip(report['branches'],(27,28)):
        if not branch['main'] or branch['main']['acceptedNonlinearIterations']==0:continue
        source=Path(branch['resultFolder']);target=ROOT/'outputs/v83_validation'/f'krylov{number}'
        assert not target.exists();target.mkdir()
        for name in ('wall_velocity_krylov_evidence.mat','wall_velocity_krylov_summary.json'):
            os.link(source/name,target/name);assert sha(source/name)==sha(target/name)
        entry=dict(method=branch['method'],sourceFolder=str(source),displayCase=target.name,
                   SHA256=sha(target/'wall_velocity_krylov_evidence.mat'),RMS=branch['main']['rmsMerit'],maximum=branch['main']['maximumResidual'])
        entries.append(entry);available.append(entry)
    assert available,'没有新接受态时须保留Q26并采用独立诊断展示，不伪造新云图。'
    chosen=min(available,key=lambda a:(a['RMS'],a['maximum']))
    binding=dict(status='byte-identical-display-aliases-of-independent-actual-temperature-study-branches',
                 branches=entries,selected=chosen,parallelHistoriesMerged=False,newFlowAcceptances=0,newProductionSolvedAreaM2=0)
    (OUT/'view_source_binding.json').write_text(json.dumps(binding,indent=2)+'\n',encoding='utf8');print(json.dumps(binding))


def bind(study_file,candidate_file,state_file,method):
    r=json.loads(study_file.read_text('utf8'));branch=next(b for b in r['branches'] if b['method']==method)
    c=loadmat(candidate_file,simplify_cells=True,variable_names=['value','detail','sourceBinding','preparationFile','preparationSHA256'])
    s=loadmat(state_file,simplify_cells=True,variable_names=['result','report']);actual=s['result'];lines=actual['linearLog'];line=lines if isinstance(lines,dict) else lines[0]
    assert c['sourceBinding']==r['sourceBinding'] and c['preparationSHA256']==r['preparationSHA256']==sha(c['preparationFile'])
    assert sha(c['sourceBinding']['sourceFile'])==c['sourceBinding']['sourceFileSHA256']
    np.testing.assert_array_equal(c['value'],line['steadyDirectionEvidence']['directionU'])
    probe=loadmat(actual['options']['steadyInitialGuess']['file'],simplify_cells=True,variable_names=['sourceBinding','directionU'])
    np.testing.assert_array_equal(c['value'],probe['directionU']);assert probe['sourceBinding']['sha256']==c['sourceBinding']['sourceFileSHA256']
    assert branch['main']['residualEvaluations']==s['report']['residualEvaluations'] and branch['completeSelectedModelRelativeDefect']<=.1
    detail=c['detail'];projected=method=='projected';expected='complete-parent-C-fixed-primitive-penalty-'+('temperature-row-nullspace-QR-LSQR' if projected else 'augmented-patch-QR-LSQR')
    assert detail['method']==expected
    target=detail['internalCTarget'] if projected else detail['options']['relativeTolerance']
    return dict(kind='temperature-comparison-QR-LSQR-v1',branchMethod=method,
                methodLabel=('分块QR＋界面温度一阶齐次投影' if projected else '分块QR自由方向')+f"；零初猜{detail['nativeIterations']}步LSQR",
                studyFile=str(study_file),studyFileSHA256=sha(study_file),candidateFile=str(candidate_file),candidateFileSHA256=sha(candidate_file),
                baseCase='krylov26',usedTotalDirectionExactlyVerified=True,hardUpdateLimitClaimed=False,
                firstSegmentIterations=detail['nativeIterations'],restartedSegmentIterations=0,
                internalCRelativeResidual=detail['relativeCResidual'],internalCTarget=target,internalCTargetMet=bool(detail['internalCTargetMet']),
                completeModelRelativeDefect=branch['completeSelectedModelRelativeDefect'],originalNSDirectionDefect=line['fullSystemRelativeDefect'],
                originalNSDirectionTarget=actual['options']['linearRelativeTolerance'],constraintApplied=projected,
                constraintFaceIDs=r['faceIDs'],constraintDerivativeK=branch['constraintDerivativeK'],
                physicalTemperatureBoundaryAdded=False,nonlinearTemperatureBoundGuaranteed=False)


def build():
    binding=json.loads((OUT/'view_source_binding.json').read_text('utf8'));selected=binding['selected'];method=selected['method']
    source=OUT/'workbench_base';folder=OUT/'workbench01';assert not folder.exists()
    prep=bind(COMP/'temperature_projection_comparison.json',COMP/method/'direction.mat',
              ROOT/'outputs/v83_validation'/selected['displayCase']/'wall_velocity_krylov_evidence.mat',method)
    shutil.copytree(source,folder);shutil.copy2(ROOT/'web/rocket-wall-velocity.js',folder/'rocket-wall-velocity.js')
    d=json.loads((folder/'wall-velocity-study-data.json').read_text('utf8'));d['cases'][2]['label']='v86 '+('温度投影' if method=='projected' else '自由')+'实际接受态（未收敛）'
    d['cases'][2]['solverDiagnostics']['directionPreparation']=prep
    report=json.loads((COMP/'temperature_projection_comparison.json').read_text('utf8'))
    labels=[]
    for b in report['branches']:
        m=b['main'];labels.append(f"{b['method']}："+(f"本段{m['acceptedNonlinearIterations']}次接受，RMS {m['rmsMerit']:.3f} /s，最大 {m['maximumResidual']:.3f} /s" if m else '未通过方向门槛，无流场接受'))
    d['temperatureProjectionResearch']=dict(selectedCase=selected['displayCase'],comparisonLabel='；'.join(labels),
        constraintLabel='原内部面19207的采样温度一阶方向约束；不固定物理温度，不保证非线性温域',
        scopeLabel='两分支从相同Q26独立推进，历史未合并；没有新增合格面积',reportURL='temperature-projection-diagnostics.html',
        comparison=report,physicalSolution=False,newProductionSolvedAreaM2=0)
    raw=json.dumps(d,ensure_ascii=False,separators=(',',':'));(folder/'wall-velocity-study-data.json').write_text(raw,encoding='utf8')
    (folder/'wall-velocity-study-data.js').write_bytes(('window.ROCKET_WALL_VELOCITY_STUDY='+raw+';\n').encode('utf8'))
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>界面温度方向研究</title>
<style>body{font:16px/1.7 system-ui,sans-serif;max-width:1350px;margin:28px auto;padding:0 24px;color:#213747;background:#f7f9fb}img{width:100%;background:white}a{color:#146da2}.note{padding:16px;border-left:4px solid #ae7c24;background:#fff3d9}</style>
<a href="index.html" target="_top">返回当前结果工作台</a><h1>真实20°箭体：界面温度切向约束</h1>
<p>在工作台内浏览时，可直接点击左侧“查看结果”回到三维云图。</p>
<p class="note">方向研究及一次实际接受不能证明整箭收敛。压力、热流和温度均为未收敛研究值，新增合格面积为0。</p>
<p>相同Q26、完整C、冻结W和penalty，各128步。受限分支显式增加H·方向=0；它不是给定壁温，也不固定最终界面温度。原NS、正性、80–2000 K及稳态门槛均保持。</p>
<img src="temperature_projection_comparison.png" alt="同预算方向缺陷、实际残差改善和内部界面非线性温度">
<p><a href="temperature_projection_comparison.pdf">科研图 PDF</a> · <a href="temperature_projection_record.json">完整比较记录</a></p></html>'''
    (folder/'temperature-projection-diagnostics.html').write_text(page,encoding='utf8')
    (folder/'temperature_projection_record.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    for ext in ('png','pdf'):shutil.copy2(OUT/f'temperature_projection_comparison.{ext}',folder/f'temperature_projection_comparison.{ext}')
    same=[]
    for file in (source/'cases').glob('*/results.json'):
        rel=file.relative_to(source);assert file.read_bytes()==(folder/rel).read_bytes();same.append(str(rel))
    out=dict(status='bound-temperature-comparison-and-actual-total-direction-in-five-step-workbench',selected=binding['selected'],
             directionPreparation=prep,originalFacetDatasetsByteIdentical=same,parallelHistoriesMerged=False,
             physicalSolution=False,newProductionSolvedAreaM2=0,browserVisualOrWebGLRuntimeVerified=False)
    (folder/'temperature_view_binding.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf8');print(json.dumps(out,ensure_ascii=False))


if __name__=='__main__':{'aliases':aliases,'build':build}[sys.argv[1]]()
