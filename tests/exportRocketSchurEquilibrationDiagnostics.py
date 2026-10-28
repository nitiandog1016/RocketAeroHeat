# -*- coding: utf-8 -*-
"""将尺度均衡的真实终态接回同源工作台，未知指标保持null，不增加已解面积。"""
from pathlib import Path
import os, json, hashlib
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
target = VAL / 'block_linear_research_summary_v4.json'
assert not target.exists()

def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            result.update(block)
    return result.hexdigest()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

case = 'schurscalereal01'
receipt = VAL / case / 'pivoted_schur_real_summary.json'
process = read(VAL / (case + '_process.json'))
assert process['status'] == 'terminal' and all(process['sourceStableAfter'].values())
audit_file = VAL / (case + '_independent_equilibration.json')
audit = read(audit_file)
assert audit['status'].startswith('passed-') and not audit['coefficientsAddedOrDeleted']
out = read(VAL / 'block_linear_research_summary_v3.json')
new_actions = 0
if not receipt.exists():
    receipt = VAL / case / 'agent_stop_record.json'
    r = read(receipt)
    assert process['exitCode'] != 0 and r['processID'] in process['engineProcessIDs']
    assert r['status'] == 'agent-resource-stopped-before-native-factorization-return-not-numerically-judged'
    assert not r['usableFactorProduced'] and not r['numericalFailureOrZeroPivotClaimed']
    assert not (VAL / case / 'pivoted_schur_factor_evidence.mat').exists()
    source_key = 'outputs\\v83_validation\\krylov17\\wall_velocity_krylov_evidence.mat'
    assert process['sourcesBefore'][source_key] == out['sourceEvidenceSHA256']
    state = '资源停止（数值未定）'
    note = ('完整S的可逆尺度及原单位还原已独立核对。库分解尚未返回，用户态CPU约' +
            format(r['userCPUSeconds'], '.1f') + '秒，私人提交约' + format(r['privateBytes']/1024**3, '.2f') +
            'GiB时因资源压力主动停止。没有完整因子；未判数值失败或零主元；右端、GMRES和原NS尚未开始。')
elif read(receipt)['candidateRejected']:
    r = read(receipt)
    assert process['exitCode'] == 0 and r['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
    assert not r['wholePreconditionerReady'] or r['stage'] != 'bounded-fixed-actions-complete'
    state = '尺度均衡后仍拒绝'
    note = ('全部原S系数的可逆行列尺度及原单位还原已独立核对。真实阶段：' + r['stage'] +
            '；错误：' + r['errorIdentifier'] + '。没有放宽门槛、补对角或更新物理场。')
else:
    r = read(receipt)
    assert process['exitCode'] == 0 and r['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
    assert r['factorReady'] and r['wholePreconditionerReady']
    complete = read(VAL / (case + '_independent_pivoted_schur.json'))
    assert complete['status'].startswith('passed-')
    state = '仅代数动作完成'
    new_actions = 3
    note = ('已生成完整五量预条件器并核对三条原C作用；同Q剩余的单次相对C残差为' +
            format(r['summaries'][0]['relativeCResidual'], '.6g') +
            '。尚未开始GMRES或原NS方向验收，不代表流场收敛或新增物理解面积。')
assert r['newActualNSEvaluations'] == r['newProductionSolvedAreaM2'] == 0
out['studies'].append(dict(
    label='可逆行列尺度 / 行主元标量Schur', caseID=case, linearSystem='C=-Jgraph',
    arnoldiSteps=0, CApplications=new_actions,
    ILUApplications=new_actions,
    innerRelativeCResidual=None, basisOrthogonalityError=None, modelRelationRelativeDifference=None,
    orthogonalityThreshold=1e-8, modelRelationThreshold=1e-8, candidateUsable=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, evidenceFile=str(receipt),
    evidenceSHA256=sha(receipt), displayState=state, note=note))
out['reversibleEquilibrationDiagnostics'] = dict(
    auditFile=str(audit_file), auditSHA256=sha(audit_file),
    rowScaleRange=audit['rowScaleRange'], columnScaleRange=audit['columnScaleRange'],
    allInputBlocks=audit['allInputBlocks'], originalScalarNonzeros=audit['originalScalarNonzeros'],
    scaledScalarNonzeros=audit['scaledScalarNonzeros'], numericalFullRankClaimed=False,
    zeroPivotNativePositionKnown=False, originalNSDirectionRetested=False)
out['status'] = 'source-bound-equilibrated-Schur-real-termination-with-unchanged-physical-field'
out['diagnosticBuilderProcessID'] = os.getpid()
out['builderSHA256'] = sha(__file__)
target.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], studies=len(out['studies']), displayState=state), ensure_ascii=False))
