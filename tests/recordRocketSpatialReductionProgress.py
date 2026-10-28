# -*- coding: utf-8 -*-
"""登记真实空间消元与粗S/F行分解的实际证据，保留完整目标及原场量身份。"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, subprocess
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
target = ROOT/'docs/validation/v83/progress_spatial_reduction.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def identity(path):
    path = Path(path)
    return dict(path=str(path.relative_to(ROOT)), sha256=sha(path), bytes=path.stat().st_size)


processes = {}
for name in ('spatialreductionunit01', 'spatialreductionreal01', 'spatialschurreal01', 'spatialreductionui01'):
    file = VAL/(name+'_process.json')
    p = read(file)
    assert p['status'] == 'terminal' and p['exitCode'] == 0 and all(p['sourceStableAfter'].values())
    processes[name] = dict(exitCode=0, sourceStable=True, engineProcessIDs=p['engineProcessIDs'], receipt=identity(file))
small = read(VAL/'spatialreductionunit01_independent_spatial_reduction.json')
real = read(VAL/'spatialreductionreal01_independent_spatial_reduction.json')
schur = read(VAL/'spatialschurreal01_independent_spatial_schur.json')
assert all(a['status'].startswith('passed-') for a in (small, real, schur))
assert not real['internalCTargetMet'] and not schur['approximateSTargetMet'] and not schur['originalCInternalTargetMet']
for audit in (small, real, schur):
    for path, digest in audit['inputSHA256'].items():
        assert sha(ROOT/path) == digest
summary_file = VAL/'spatialreductionreal01/spatial_reduction_real_summary.json'
schur_file = VAL/'spatialschurreal01/spatial_schur_krylov_summary.json'
assert summary_file.stat().st_size < 128*1024 and schur_file.stat().st_size < 128*1024
summary, schur_summary = read(summary_file), read(schur_file)
assert summary['innerRunComplete'] and schur_summary['numericalRunComplete'] and not summary['candidateRejected'] and not schur_summary['candidateRejected']
assert not summary['physicalSolution'] and not schur_summary['physicalSolution']
assembly = summary['assemblyReport']
original_ff = assembly['originalBlocks']-assembly['CFBlocks']-assembly['FCBlocks']-assembly['CCBlocks']
within_ff = summary['FFactorReport']['includedOriginalFFBlocks']
assert original_ff == 1285801 and within_ff == 603745 and original_ff-within_ff == 682056
ui = read(VAL/'workbench25/independent_spatial_reduction_study_view.json')
node = read(VAL/'workbench25/node_wall_velocity_summary.json')
assert ui['status'].startswith('passed-') and ui['researchRows'] == 14
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not ui['browserVisualOrWebGLRuntimeVerified'] and not node['browserVisualOrWebGLRuntimeVerified']
for item in ui['sourceFacetPayload']:
    assert sha(ROOT/item['path']) == item['sha256']
old = read(ROOT/'docs/validation/v83/progress_patch_hierarchy.json')
production = old['production']
for key in ('latestEndpoint', 'originalNSSource'):
    assert sha(ROOT/production[key]['path']) == production[key]['sha256']
assert subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip() == '42601c9'
input_file = VAL/'spatialreductionpreflight01/spatial_reduction_input.mat'
assert sha(input_file) == 'f7a581f6467423799b0d771923718d56a678040b01b65e4afb9ce4aa83e5df20'
backend = ['buildRocketSpatialSplit.m', 'factorRocketSpatialFBlocks.m', 'solveRocketSpatialFBlocks.m',
    'applyRocketRectangularCoupling.m', 'buildRocketSpatialReduction.m', 'applyRocketSpatialSchur.m',
    'applyRocketSpatialReduction.m', 'solveRocketSpatialReductionKrylov.m', 'compactRocketSpatialReductionReport.m']
tests = ['testRocketSpatialReduction.m', 'studyRocketSpatialReduction.m', 'studyRocketSpatialSchurKrylov.m',
    'prepareRocketSpatialReductionInput.py', 'auditRocketSpatialReductionUnit.py', 'auditRocketSpatialReductionReal.py',
    'auditRocketSpatialSchurKrylov.py', 'exportRocketSpatialReductionDiagnostics.py', 'auditRocketSpatialReductionWorkbench.py']
out = dict(status='tested-real-five-field-spatial-reduction-and-verified-cross-patch-F-error-on-saved-Schur-candidate',
    version=83, recordedUTC=datetime.now(timezone.utc).isoformat(), processID=os.getpid(),
    previousGoalTurnClassification='progress-real-seven-level-test-independent-audits-and-fixed-spatial-CF-structural-preflight',
    currentGoalTurnClassification='progress-new-complete-MATLAB-spatial-reduction-real-test-and-S-F-error-attribution-with-updated-workbench',
    processes=processes, smallDomainAudit=small, fullSpatialReductionRealAudit=real, separateSchurKrylovAudit=schur,
    compactNumericalSummary=summary, compactSchurDiagnosticSummary=schur_summary,
    compactSummaryReceipts=[identity(summary_file), identity(schur_file)],
    frozenPreflightMATLABInput=identity(input_file), numericalFailuresInThisSliceClaimed=False,
    FBlockTopology=dict(allOriginalFFBlocks=original_ff, withinOriginalPatchFFBlocks=within_ff,
        crossPatchFFBlocksKeptInCompleteC=original_ff-within_ff, originalPhysicalCVs=217272,
        FOriginalCVs=187749, coarseOriginalCVs=29523, physicalGridChanged=False),
    codeSourceIdentity=[identity(ROOT/'internalFlow/cutCell'/name) for name in backend]+[identity(ROOT/'tests'/name) for name in tests],
    methodDocumentation=identity(ROOT/'docs/SPATIAL_REDUCTION_V83.md'), userInterface=ui, nodeDataVerification=node,
    production=production, delivery=old['delivery'],
    workAccounting=dict(mainSpatialStudyOriginalCApplications=35, separateSStudyOriginalCDiagnosticApplications=1,
        combinedOriginalCApplications=36, combinedPackedSApplications=38, combinedFLocalApplications=71,
        combinedCoarseLocalApplications=66, combinedRectangularApplications=77,
        newOriginalMatrixAssemblies=0, newFInverseConstructions=1, newApproximateSchurAssemblies=1,
        separateSStudyNewFInverseConstructions=0, separateSStudyNewSchurAssemblies=0,
        newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0),
    nextAction='Reuse the audited spatial F split, local F inverse, approximate S, frozen coarse factor and saved S/F candidate. The saved diagnostic F rows are dominated by cross-patch F-F coupling (10.870568 relative to original RHS) while the within-patch F equation defect is only 1.30e-12; this is evidence for that vector, not a global cause or rank proof. Build/check the complete all-five-field A_FF sub-operator from the original C, retaining all 682056 cross-patch FF blocks; first dense small reference, then one bounded same-RHS test using the existing fixed local F inverse. Keep it a diagnostic until consistency with the condensed S and original full-C back-substitution is established. Assess overlap or F-only strength grouping from actual full-F coupling/resource evidence if needed. Do not merely extend coarse-S iterations, silently replace M_F while claiming the old S identity, rerun old group-constant hierarchies/global momentum ILU/stopped scalar ILUTP, or relax original gates. All original CVs and original full C remain. Require complete selected-J and original-NS total-direction verification with original norm(R)/0.1 and state/temperature/descent/spatial/steady/flux conditions before any new accepted state; regional physical validation and scientific UI remain mandatory.',
    goalStatus='active', wholeGoalComplete=False, recorderSHA256=sha(__file__))
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], manifest=identity(target), researchRows=14)))
