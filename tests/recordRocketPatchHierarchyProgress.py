# -*- coding: utf-8 -*-
"""登记实际七层失败收益、独立核对和空间消元预检，不夸大为新增区域能力。"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, subprocess
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
target = ROOT/'docs/validation/v83/progress_patch_hierarchy.json'
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
for name in ('patchhierarchyunit01', 'patchhierarchyunit02', 'patchhierarchyreal01', 'patchhierarchyui01'):
    file = VAL/(name+'_process.json')
    p = read(file)
    assert p['status'] == 'terminal' and p['exitCode'] == 0 and all(p['sourceStableAfter'].values())
    processes[name] = dict(exitCode=0, sourceStable=True, engineProcessIDs=p['engineProcessIDs'], receipt=identity(file))
small = [read(VAL/(name+'_independent_patch_hierarchy.json')) for name in ('patchhierarchyunit01', 'patchhierarchyunit02')]
real = read(VAL/'patchhierarchyreal01_independent_patch_hierarchy.json')
assert all(a['status'].startswith('passed-') for a in small+[real])
assert not real['internalCTargetMet'] and not real['fullCResidualImprovedComparedWithLocalOnly']
for audit in small+[real]:
    for path, digest in audit['inputSHA256'].items():
        assert sha(ROOT/path) == digest
summary_file = VAL/'patchhierarchyreal01/patch_hierarchy_real_summary.json'
assert summary_file.stat().st_size < 128*1024
summary = read(summary_file)
assert summary['innerRunComplete'] and not summary['candidateRejected'] and not summary['physicalSolution']
assert summary['newRootLocalFactorConstructions'] == summary['newActualNSEvaluations'] == summary['newFlowAcceptances'] == 0
report = summary['innerReport']
construction_C = sum(t['fineCVerificationApplications']+t['coarseCVerificationApplications'] for t in summary['hierarchyReport']['transfers'])
all_C = report['solverCApplications']+sum(report['preconditionerCApplicationsByLevel'])+sum(summary['singleWork']['CApplicationsByLevel'])+summary['extraRootCDiagnosticApplications']+construction_C
root_C = report['totalRootCApplications']+summary['singleWork']['CApplicationsByLevel'][0]+summary['extraRootCDiagnosticApplications']+1
assert all_C == 245 and root_C == 69
preflight_file = VAL/'spatialreductionpreflight01/spatial_reduction_preflight.json'
preflight = read(preflight_file)
assert not preflight['numericalLocalFFInversesConstructed'] and not preflight['numericalSchurCoefficientsConstructed'] and not preflight['coarseInverseOrConvergenceProven']
assert sha(ROOT/preflight['topologyEvidenceFile']) == preflight['topologyEvidenceSHA256']
for path, digest in preflight['inputSHA256'].items():
    assert sha(ROOT/path) == digest
failure_file = VAL/'patchhierarchyunit02_first_audit_failure.json'
failure = read(failure_file)
assert failure['exitCode'] == 1 and not failure['numericEvidenceRecomputed']
old_audit_epoch = VAL/'source_epochs/passed_patchhierarchyunit01_audit_before_deep.py'
assert sha(old_audit_epoch) == small[0]['scriptSHA256']
assert sha(ROOT/'tests/auditRocketPatchHierarchyUnit.py') == small[1]['scriptSHA256']
ui = read(VAL/'workbench24/independent_patch_hierarchy_study_view.json')
node = read(VAL/'workbench24/node_wall_velocity_summary.json')
assert ui['status'].startswith('passed-') and ui['researchRows'] == 12
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not ui['browserVisualOrWebGLRuntimeVerified'] and not node['browserVisualOrWebGLRuntimeVerified']
for item in ui['sourceFacetPayload']:
    assert sha(ROOT/item['path']) == item['sha256']
old = read(ROOT/'docs/validation/v83/progress_coupled_patches.json')
production = old['production']
for key in ('latestEndpoint', 'originalNSSource'):
    assert sha(ROOT/production[key]['path']) == production[key]['sha256']
assert subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip() == '42601c9'
backend = ['restrictRocketPatchResidual.m', 'prolongRocketPatchCorrection.m', 'coarsenRocketPackedGraph.m',
    'buildRocketPatchHierarchy.m', 'applyRocketPatchHierarchy.m', 'solveRocketPatchHierarchyKrylov.m', 'compactRocketPatchHierarchyReport.m']
tests = ['testRocketPatchHierarchy.m', 'studyRocketPatchHierarchy.m', 'auditRocketPatchHierarchyUnit.py',
    'auditRocketPatchHierarchyReal.py', 'preflightRocketSpatialReduction.py', 'exportRocketPatchHierarchyDiagnostics.py',
    'auditRocketPatchHierarchyWorkbench.py']
epochs = ['pre_deep_patch_hierarchy_testRocketPatchHierarchy.m', 'pre_single_guard_studyRocketPatchHierarchy.m',
    'passed_patchhierarchyunit01_audit_before_deep.py', 'pre_singleton_auditRocketPatchHierarchyUnit.py']
out = dict(status='tested-seven-level-complete-coefficient-hierarchy-with-real-coarse-amplification-and-spatial-reduction-preflight',
    version=83, recordedUTC=datetime.now(timezone.utc).isoformat(), processID=os.getpid(),
    previousGoalTurnClassification='progress-real-local-patch-study-independent-audits-and-workbench23',
    currentGoalTurnClassification='progress-new-MATLAB-hierarchy-real-run-independent-evidence-and-next-candidate-structural-preflight',
    processes=processes, smallDomainAudits=small, realAudit=real, compactNumericalSummary=summary,
    compactSummaryReceipt=identity(summary_file),
    retainedIndependentReaderFailure=failure, readerFailureReceipt=identity(failure_file),
    readerRepairDidNotRecomputeNumerics=True, sourceEpochs=[identity(VAL/'source_epochs'/name) for name in epochs],
    codeSourceIdentity=[identity(ROOT/'internalFlow/cutCell'/name) for name in backend]+[identity(ROOT/'tests'/name) for name in tests],
    methodDocumentation=identity(ROOT/'docs/PATCH_HIERARCHY_V83.md'),
    officialDocumentationAndAuthorAbstractReview=identity(ROOT/'docs/validation/v83/reduction_method_sources.json'),
    spatialCFStructuralResourcePreflight=preflight, userInterface=ui, nodeDataVerification=node,
    production=production, delivery=old['delivery'],
    workAccounting=dict(newRootLocalFactorConstructions=0, newNonterminalCoarseLocalFactors=5,
        newRCPTransfers=6, constructionCApplicationsAcrossLevels=construction_C,
        solverRootCApplications=33, preconditionerCApplicationsAcrossLevels=192,
        singleProbePreconditionerCApplicationsAcrossLevels=6, rootDiagnosticCApplications=2,
        allStudyRootCApplications=root_C, allStudyCApplicationsAcrossLevels=all_C,
        independentAuditFixedVcycleApplications=1, independentAuditRootCDiagnostics=2,
        spatialPreflightNumericalInverseConstructions=0, spatialPreflightNewCApplications=0,
        newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0),
    nextAction='Do not repeat the tested group-constant seven-level V-cycle, adjust two scalar weights and rerun, or blindly add iterations, old normal scans, old projection spaces, global momentum ILU or stopped scalar ILUTP. Official MGR/PCAIR guides and the lAIR author abstract motivate a spatial C/F local-elimination candidate, not a claimed reproduction. Reuse the saved fixed original-CV representatives and complete Schur structural pattern from spatialreductionpreflight01. First implement and check a small dense five-field spatial split, certified within-patch F inverse, explicit approximate Schur D-C*MF*B and complete back-substitution, accounting for ignored cross-patch F-F connections only in the preconditioner. Preserve all original CVs and outer C. The 310074792-byte partial estimate excludes original C, CF/FC caches, temporaries, deeper levels and Krylov basis; establish full construction/storage bounds before real coefficients. Numerical F inverses and Schur coefficients remain unbuilt. Full selected-J/original-NS total direction, 0.1/norm(R), state/temperature/descent/spatial/steady/flux gates and physical regional validation remain mandatory.',
    goalStatus='active', wholeGoalComplete=False, recorderSHA256=sha(__file__))
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], manifest=identity(target), researchRows=12)))
