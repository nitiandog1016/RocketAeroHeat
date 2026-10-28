# -*- coding: utf-8 -*-
"""接入完整F与重分组对照，所有F残差与固定C原矩阵诊断保持分开。"""
from pathlib import Path
import os, json, hashlib
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
target = VAL/'block_linear_research_summary_v10.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


out = read(VAL/'block_linear_research_summary_v9.json')
for case, expected in [('fullFunit01', 0), ('fullFreal01', 0), ('fullFstrengthunit01', 0), ('fullFstrengthunit02', 0), ('fullFstrengthreal01', 1), ('fullFstrengthreal02', 0)]:
    p = read(VAL/(case+'_process.json'))
    assert p['status'] == 'terminal' and p['exitCode'] == expected and all(p['sourceStableAfter'].values())
for case, filename, audit_name, label, state, note in [
    ('fullFreal01', 'full_F_real_summary.json', 'fullFreal01_independent_full_F.json',
     '完整A_FF / 全部跨组F连接', '完整F仍未达标',
     '此行指标属于完整F子系统，原C列保持未知。187749F/1285801块全部保留，32步F残差0.278621（原局部值0.397038）；固定旧C值回代原C为9.65848，不是原C迭代或原NS通过。'),
    ('fullFstrengthreal02', 'full_F_strength_summary.json', 'fullFstrengthreal02_independent_full_F_strength.json',
     '完整A_FF / F内强连接分组', 'F内重分组未改善',
     '同系统/右端32步，F残差0.283053，较旧0.278621变差。27461组、1组CV对角替代；有效纳入145730个旧跨组块，数量增加不代表收敛改善。固定C回代原C为13.8211，旧近似S未冒称一致。')]:
    file, audit_file = VAL/case/filename, VAL/audit_name
    s, audit = read(file), read(audit_file)
    assert audit['status'].startswith('passed-') and s['innerRunComplete'] and s['parentDiagnosticComplete'] and not s['candidateRejected']
    assert not s['fullFInternalTargetMet'] and s['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
    cycle = s['innerFullFReport']['krylov']['cycles'][0]
    out['studies'].append(dict(label=label, caseID=case, linearSystem='exact-complete-spatial-A_FF-with-separate-fixed-C-parent-diagnostic',
        arnoldiSteps=32, CApplications=1, ILUApplications=0, innerRelativeCResidual=None,
        basisOrthogonalityError=cycle['basisOrthogonalityError'], modelRelationRelativeDifference=cycle['arnoldiRelationRelativeDifference'],
        orthogonalityThreshold=1e-8, modelRelationThreshold=1e-8, candidateUsable=False,
        newActualNSEvaluations=0, newProductionSolvedAreaM2=0, evidenceFile=str(file), evidenceSHA256=sha(file), displayState=state, note=note))
    out['fullFDiagnostics' if case == 'fullFreal01' else 'fullFStrengthDiagnostics'] = dict(auditFile=str(audit_file), auditSHA256=sha(audit_file),
        fullFResidual=audit['finalFullFResidual'], parentOriginalCResidual=s['parentDiagnostic']['relativeParentOriginalCResidual'],
        fullFInternalTargetMet=False, originalCAndFResidualDenominatorsKeptDistinct=True,
        fullFArnoldiBasisSaved=False, independentAllRealFullFArnoldiColumnsClaimed=False,
        originalNSOrRegionalQualificationProven=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0)
strength = read(VAL/'fullFstrengthreal02_independent_full_F_strength.json')
out['fullFStrengthDiagnostics']['strengthCoverage'] = strength['strengthCoverage']
out['fullFStrengthDiagnostics']['fullFResidualImproved'] = False
failure_file = VAL/'fullFstrengthreal01_shape_failure.json'
out['fullFStrengthDiagnostics']['retainedShapeFailure'] = dict(evidenceFile=str(failure_file), evidenceSHA256=sha(failure_file),
    driverExitCode=1, iterationHadStarted=False, numericalFactorWasPersisted=False, repairedCoverageRowAndColumnTestPassed=True)
preflight_file = VAL/'fullFoverlappreflight01/full_F_overlap_preflight.json'
p = read(preflight_file)
assert not p['numericalOverlapInversesConstructed'] and not p['convergenceOrRegionalQualificationProven']
out['restrictedOverlapPreflight'] = dict(evidenceFile=str(preflight_file), evidenceSHA256=sha(preflight_file),
    maximumAddedFNeighbors=2, maximumOverlapCells=10, factorArrayBoundBytes=p['conservativeFactorArrayBoundBytes'],
    proposedFactorBudgetBytes=p['proposedFactorBudgetBytes'], completeConstructionOrProcessBudgetProven=False,
    numericalOverlapInversesConstructed=False, physicalQualificationProven=False)
out['status'] = 'source-bound-complete-F-and-strength-repartition-unresolved-with-restricted-overlap-preflight'
out['diagnosticBuilderProcessID'] = os.getpid()
out['builderSHA256'] = sha(__file__)
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], researchRows=len(out['studies']), bytes=target.stat().st_size)))
