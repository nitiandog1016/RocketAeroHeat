# -*- coding: utf-8 -*-
"""接入完整C空间消元及分离的S/F诊断，防止不同残差分母混列。"""
from pathlib import Path
import os, json, hashlib
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
target = VAL/'block_linear_research_summary_v9.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


out = read(VAL/'block_linear_research_summary_v8.json')
summary_file = VAL/'spatialreductionreal01/spatial_reduction_real_summary.json'
audit_file = VAL/'spatialreductionreal01_independent_spatial_reduction.json'
s, a = read(summary_file), read(audit_file)
assert a['status'].startswith('passed-') and s['innerRunComplete'] and not s['internalCTargetMet']
assert s['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
for name in ('spatialreductionunit01', 'spatialreductionreal01', 'spatialschurreal01'):
    process = read(VAL/(name+'_process.json'))
    assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
cycle = s['innerReport']['krylov']['cycles'][0]
out['studies'].append(dict(label='完整C五变量 / 空间局部F消元', caseID='spatialreductionreal01',
    linearSystem='original-full-C-with-spatial-local-F-approximate-Schur-block-map', arnoldiSteps=32,
    CApplications=33, ILUApplications=0, innerRelativeCResidual=a['finalOriginalCResidual'],
    basisOrthogonalityError=cycle['basisOrthogonalityError'], modelRelationRelativeDifference=cycle['arnoldiRelationRelativeDifference'],
    orthogonalityThreshold=1e-8, modelRelationThreshold=1e-8, candidateUsable=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, evidenceFile=str(summary_file), evidenceSHA256=sha(summary_file),
    displayState='空间回代有限，内部未收敛',
    note='全部217272CV和五量保留，187749F/29523C。近似S=D-E*M_F*B全部254825块独立核对；32步完整C残差0.983921，未达1e-3。M_F不是完整F逆，跨组F-F仍在原C中。'))
out['spatialReductionDiagnostics'] = dict(auditFile=str(audit_file), auditSHA256=sha(audit_file),
    originalCVs=217272, FOriginalCVs=187749, coarseOriginalCVs=29523,
    approximateSchurBlocks=254825, completeCoefficientAndInverseAuditPassed=True,
    maximumSchurProductAndSumRoundingBudgetRatio=a['maximumSchurRoundingBudgetRatio'],
    FInverseDecimal80Checks=a['FLocalInverseChecks']['decimal80CheckedBlocks'],
    coarseInverseDecimal80Checks=a['coarseLocalInverseChecks']['decimal80CheckedBlocks'],
    originalCResidual=a['finalOriginalCResidual'], internalCTarget=1e-3, internalCTargetMet=False,
    numericalStudyOriginalCApplications=35, numericalStudyFLocalApplications=69,
    numericalStudyCoarseLocalApplications=33, numericalStudyRectangularApplications=75, numericalStudyPackedSchurProbeApplications=3,
    reductionArrayBytes=385297642, coarseFactorArrayBytes=45073794, combinedArrayBytes=430371436,
    processPeakMemoryClaimed=False, fullFFInverseOrExactSchurClaimed=False,
    realFullArnoldiBasisSaved=False, independentAllRealArnoldiColumnAuditClaimed=False,
    originalNSConvergenceOrRegionalQualificationProven=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0)
schur_summary_file = VAL/'spatialschurreal01/spatial_schur_krylov_summary.json'
schur_audit_file = VAL/'spatialschurreal01_independent_spatial_schur.json'
t, u = read(schur_summary_file), read(schur_audit_file)
assert u['status'].startswith('passed-') and t['numericalRunComplete'] and not u['approximateSTargetMet'] and not u['originalCInternalTargetMet']
assert t['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
out['studies'].append(dict(label='近似S / 单独32步及原C回代诊断', caseID='spatialschurreal01',
    linearSystem='spatial-approximate-Schur-S-with-separate-parent-C-diagnostic', arnoldiSteps=32,
    CApplications=1, ILUApplications=0, innerRelativeCResidual=None,
    basisOrthogonalityError=u['basisOrthogonalityError'], modelRelationRelativeDifference=u['finalSCombinationRelation'],
    orthogonalityThreshold=1e-8, modelRelationThreshold=1e-8, candidateUsable=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, evidenceFile=str(schur_summary_file), evidenceSHA256=sha(schur_summary_file),
    displayState='S未达标，F跨组误差显著',
    note='此行正交/组合指标属于近似S，不是完整C迭代。S残差0.232219→0.0239202；单次代回原C残差10.8747，99.9243%剩余平方范数在F行。该向量跨组F-F贡献10.8706、局部逆缺陷1.30e-12，未更新物理场。'))
out['spatialSchurDiagnostic'] = dict(auditFile=str(schur_audit_file), auditSHA256=sha(schur_audit_file),
    originalSAndCResidualDenominatorsKeptDistinct=True, fullSArnoldiBasisSaved=True,
    all32ApproximateSPhysicalDirectionsIndependentlyVerified=True,
    singleApproximateSResidual=u['singleApproximateSResidual'], finalApproximateSResidual=u['finalApproximateSResidual'],
    parentOriginalCResidual=u['finalOriginalCResidual'], FRowsRemainingSquaredNormFraction=u['FRowsRemainingSquaredNormFraction'],
    localFInverseEquationDefectRelativeToOriginalRHS=u['localFInverseEquationDefectRelativeToOriginalRHS'],
    crossPatchFFImageRelativeToOriginalRHS=u['crossPatchFFImageRelativeToOriginalRHS'],
    FRowExactDecompositionDifference=u['FRowExactDecompositionDifference'],
    attributionScope=u['attributionScope'], approximateSTargetMet=False, originalCInternalTargetMet=False,
    newFInverseConstructions=0, newSchurAssemblies=0, numericalStudyPackedSApplications=35,
    numericalStudyOriginalCDiagnosticApplications=1, newActualNSEvaluations=0, newProductionSolvedAreaM2=0)
out['status'] = 'source-bound-spatial-reduction-and-separated-Schur-F-error-diagnostics-unresolved-original-C'
out['diagnosticBuilderProcessID'] = os.getpid()
out['builderSHA256'] = sha(__file__)
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], researchRows=len(out['studies']), bytes=target.stat().st_size)))
