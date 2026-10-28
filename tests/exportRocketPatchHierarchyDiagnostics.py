# -*- coding: utf-8 -*-
"""把真实七层未收敛结果接回科研工作台，不修改原面物理场。"""
from pathlib import Path
import os, json, hashlib
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
target = VAL/'block_linear_research_summary_v8.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


out = read(VAL/'block_linear_research_summary_v7.json')
summary_file = VAL/'patchhierarchyreal01/patch_hierarchy_real_summary.json'
audit_file = VAL/'patchhierarchyreal01_independent_patch_hierarchy.json'
s, a = read(summary_file), read(audit_file)
assert a['status'].startswith('passed-') and s['innerRunComplete'] and not s['candidateRejected']
assert not s['internalCTargetMet'] and not a['fullCResidualImprovedComparedWithLocalOnly']
assert s['sourceBinding']['sourceFileSHA256'] == out['sourceEvidenceSHA256']
for name in ('patchhierarchyunit01', 'patchhierarchyunit02', 'patchhierarchyreal01'):
    process = read(VAL/(name+'_process.json'))
    assert process['status'] == 'terminal' and process['exitCode'] == 0 and all(process['sourceStableAfter'].values())
cycle = s['innerReport']['krylov']['cycles'][0]
out['studies'].append(dict(label='完整C五变量 / 强连接七层V循环', caseID='patchhierarchyreal01',
    linearSystem='C=-Jgraph-with-fixed-all-five-field-unsmoothed-aggregation-V-cycle', arnoldiSteps=32,
    CApplications=33, ILUApplications=0, innerRelativeCResidual=a['finalFullCResidual'],
    basisOrthogonalityError=cycle['basisOrthogonalityError'], modelRelationRelativeDifference=cycle['arnoldiRelationRelativeDifference'],
    orthogonalityThreshold=1e-8, modelRelationThreshold=1e-8, candidateUsable=False,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, evidenceFile=str(summary_file), evidenceSHA256=sha(summary_file),
    displayState='粗校正放大，内部未收敛',
    note='全部原CV保留，7层217272→3。所有粗层系数和新局部逆独立通过；单V循环原C残差76393，32步后0.999657，较局部块0.986936变差。未达1e-3，未执行原NS或更新云图。'))
out['patchHierarchyDiagnostics'] = dict(auditFile=str(audit_file), auditSHA256=sha(audit_file),
    levelCellCounts=a['levelCellCounts'], allOriginalCVsAndFiveFieldsKept=True,
    allSixCompleteCoefficientRCPTransfersVerified=True, rootGraphAndLocalCacheExactlyUnchanged=True,
    maximumRCPSourceActionDifference=max(t['RCPActionRelativeDifference'] for t in s['hierarchyReport']['transfers']),
    coefficientAggregationRoundingBoundsIndependentlyChecked=True,
    secondAlgebraicLevelOriginalCVDiagonalFallbackPatches=4,
    singleVcycleFullCResidual=a['singleFullCResidual'], internalFullCResidual=a['finalFullCResidual'],
    localOnlyReferenceFullCResidual=a['previousLocalOnlyFullCResidual'], fullCResidualImproved=False,
    internalCTargetMet=False, solverRootCApplications=33, preconditionerRootCApplications=32,
    constructionRootCApplications=1, singleProbeRootCApplications=1, extraRootCDiagnosticApplications=2,
    allStudyRootCApplications=69, allStudyCApplicationsAcrossLevels=245,
    realFullArnoldiBasisSaved=False, independentAllRealArnoldiColumnAuditClaimed=False,
    sameRhsTwoImageDiagnostic=a['sameRhsTwoImageDiagnostic'],
    originalNSConvergenceOrRegionalQualificationProven=False, newActualNSEvaluations=0, newProductionSolvedAreaM2=0)
preflight_file = VAL/'spatialreductionpreflight01/spatial_reduction_preflight.json'
p = read(preflight_file)
assert not p['numericalLocalFFInversesConstructed'] and not p['coarseInverseOrConvergenceProven']
out['spatialReductionPreflight'] = dict(evidenceFile=str(preflight_file), evidenceSHA256=sha(preflight_file),
    coarseOriginalCVs=p['coarseOriginalCVs'], FOriginalCVs=p['FOriginalCVs'],
    retainedSchurStructuralBlocks=p['retainedSchurStructuralBlocks'],
    partialCacheEstimateBytes=p['inversePlusCoarseCSRBytes'], completeMemoryBudgetClaimed=False,
    numericalFactorsConstructed=False, convergenceOrPhysicalQualificationProven=False)
out['status'] = 'source-bound-verified-seven-level-coarse-correction-amplification-with-unresolved-original-C'
out['diagnosticBuilderProcessID'] = os.getpid()
out['builderSHA256'] = sha(__file__)
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], researchRows=len(out['studies']), bytes=target.stat().st_size)))
