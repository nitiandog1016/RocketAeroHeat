# -*- coding: utf-8 -*-
"""登记完整F与重新分组实测、统计形状修复及尚未求逆的有限重叠预检。"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, subprocess
import numpy as np
from scipy.io import loadmat
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
target = ROOT/'docs/validation/v83/progress_full_F.json'
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
names = ('fullFunit01', 'fullFreal01', 'fullFstrengthunit01', 'fullFstrengthreal01', 'fullFstrengthunit02', 'fullFstrengthreal02', 'fullFui01')
for name in names:
    file = VAL/(name+'_process.json')
    p = read(file)
    expected = 1 if name == 'fullFstrengthreal01' else 0
    assert p['status'] == 'terminal' and p['exitCode'] == expected and all(p['sourceStableAfter'].values())
    processes[name] = dict(exitCode=expected, sourceStable=True, engineProcessIDs=p['engineProcessIDs'], receipt=identity(file))
small = read(VAL/'fullFunit01_independent_full_F.json')
small_strength = [read(VAL/(name+'_independent_full_F_strength.json')) for name in ('fullFstrengthunit01', 'fullFstrengthunit02')]
real = read(VAL/'fullFreal01_independent_full_F.json')
strength = read(VAL/'fullFstrengthreal02_independent_full_F_strength.json')
assert all(a['status'].startswith('passed-') for a in [small, real, strength]+small_strength)
assert not real['fullFInternalTargetMet'] and not strength['fullFInternalTargetMet'] and not strength['fullFResidualImproved']
for audit in [small, real, strength]+small_strength:
    for path, digest in audit['inputSHA256'].items():
        assert sha(ROOT/path) == digest
first = loadmat(VAL/'fullFstrengthunit01/full_F_strength_unit_evidence.mat', simplify_cells=True)
second = loadmat(VAL/'fullFstrengthunit02/full_F_strength_unit_evidence.mat', simplify_cells=True)
for a, b in zip(first['cases'], second['cases']):
    for key in ('P', 'A', 'rhs', 'single', 'value', 'complete', 'parentImage'):
        np.testing.assert_array_equal(a[key], b[key])
    for aa, bb in zip(a['factor']['buckets'], b['factor']['buckets']):
        for key in ('inverse', 'members', 'patchIDs', 'certifiedBounds', 'fallback'):
            np.testing.assert_array_equal(aa[key], bb[key])
assert sha(VAL/'source_epochs/passed_fullFstrengthunit01_audit_before_coverage.py') == small_strength[0]['scriptSHA256']
summary_files = [VAL/'fullFreal01/full_F_real_summary.json', VAL/'fullFstrengthreal02/full_F_strength_summary.json']
assert all(path.stat().st_size < 128*1024 for path in summary_files)
summary, strength_summary = [read(path) for path in summary_files]
assert summary['parentDiagnosticComplete'] and strength_summary['parentDiagnosticComplete'] and not summary['candidateRejected'] and not strength_summary['candidateRejected']
failure_file = VAL/'fullFstrengthreal01_shape_failure.json'
failure = read(failure_file)
assert failure['exitCode'] == 1 and not failure['localFactorPersistedBeforeFailure'] and failure['newFullFMatrixApplications'] == 0
assert not (VAL/'fullFstrengthreal01/full_F_strength_factor_evidence.mat').exists()
assert not (VAL/'fullFstrengthreal01/full_F_strength_inner_evidence.mat').exists()
assert 'RECYCLED_GMRES' not in (VAL/'fullFstrengthreal01_stdout.log').read_text(encoding='utf8', errors='replace')
preflight_file = VAL/'fullFoverlappreflight01/full_F_overlap_preflight.json'
preflight = read(preflight_file)
assert preflight['factorBoundFitsProposedBudget'] and not preflight['constructionAndWholeRuntimeBudgetProven'] and not preflight['numericalOverlapInversesConstructed']
assert sha(ROOT/preflight['topologyMATFile']) == preflight['topologyMATFileSHA256']
for path, digest in preflight['inputSHA256'].items():
    assert sha(ROOT/path) == digest
ui = read(VAL/'workbench26/independent_full_F_study_view.json')
node = read(VAL/'workbench26/node_wall_velocity_summary.json')
assert ui['status'].startswith('passed-') and ui['researchRows'] == 16
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
for item in ui['sourceFacetPayload']:
    assert sha(ROOT/item['path']) == item['sha256']
old = read(ROOT/'docs/validation/v83/progress_spatial_reduction.json')
production = old['production']
for key in ('latestEndpoint', 'originalNSSource'):
    assert sha(ROOT/production[key]['path']) == production[key]['sha256']
assert subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip() == '42601c9'
backend = ['extractRocketFullFGraph.m', 'solveRocketFullFKrylov.m', 'compactRocketFullFReport.m', 'rocketFullFStrengthCoverage.m']
tests = ['testRocketFullFGraph.m', 'studyRocketFullFGraph.m', 'testRocketFullFStrength.m', 'studyRocketFullFStrength.m',
    'auditRocketFullFUnit.py', 'auditRocketFullFReal.py', 'auditRocketFullFStrengthUnit.py', 'auditRocketFullFStrengthReal.py',
    'preflightRocketFullFOverlap.py', 'exportRocketFullFDiagnostics.py', 'auditRocketFullFWorkbench.py']
epochs = ['pre_F_strength_coverage_shape_studyRocketFullFStrength.m', 'pre_F_strength_coverage_testRocketFullFStrength.m', 'passed_fullFstrengthunit01_audit_before_coverage.py']
out = dict(status='tested-exact-full-FF-and-F-only-repartition-with-unresolved-residuals-and-bounded-overlap-preflight', version=83,
    recordedUTC=datetime.now(timezone.utc).isoformat(), processID=os.getpid(),
    previousGoalTurnClassification='progress-spatial-reduction-real-run-and-independent-cross-F-error-attribution',
    currentGoalTurnClassification='progress-exact-full-F-real-test-strength-repartition-comparison-shape-repair-and-overlap-preflight',
    processes=processes, smallFullFAudit=small, smallStrengthAudits=small_strength, realFullFAudit=real, realStrengthAudit=strength,
    compactNumericalSummary=summary, compactStrengthSummary=strength_summary, compactSummaryReceipts=[identity(path) for path in summary_files],
    retainedShapeFailure=failure, retainedShapeFailureReceipt=identity(failure_file),
    shapeRepair=dict(rowAndColumnInputTestPassed=True, smallNumericCachesAndSolutionsExactlyUnchanged=True,
        failedRequestNotClaimedAsActualAllocation=True, failedRunFIterationNotStarted=True, failedFactorNotPersisted=True),
    sourceEpochs=[identity(VAL/'source_epochs'/name) for name in epochs],
    codeSourceIdentity=[identity(ROOT/'internalFlow/cutCell'/name) for name in backend]+[identity(ROOT/'tests'/name) for name in tests],
    restrictedOverlapStructuralAndStoragePreflight=preflight,
    methodDocumentation=identity(ROOT/'docs/FULL_F_SUBSYSTEM_V83.md'), userInterface=ui, nodeDataVerification=node,
    production=production, delivery=old['delivery'],
    workAccounting=dict(exactFullFGraphExtractions=1, fullFNumericalStudyApplications=35, FStrengthNumericalStudyApplications=35,
        combinedCompleteFApplications=70, originalCDiagnosticApplications=2,
        existingLocalFApplications=32, newStrengthLocalApplications=33,
        newStrengthFactorConstructionsPersisted=1, failedStrengthFactorConstructionNotPersisted=1,
        failedShapeRunFApplications=0, overlapPreflightNewMatrixApplications=0, overlapNumericalInverseConstructions=0,
        newApproximateSchurAssemblies=0, newSelectedJApplications=0, newActualNSEvaluations=0,
        newFlowAcceptances=0, newProductionSolvedAreaM2=0),
    nextAction='Do not repeat either tested 32-step nonoverlapping F grouping or blindly add directions, merely tune the two saved candidate weights, old group-constant hierarchies, global momentum ILU or stopped scalar ILUTP. Exact full A_FF is now persisted and audited. Reuse its file and the frozen core/at-most-two-neighbor overlap members from fullFoverlappreflight01. First implement a fixed extended-local-inverse map returning only original core rows; prove core ownership/extension restriction, all-five-field dense references, linearity, original-F action and explicit rejection/fallback policies on small tests. Establish complete construction/storage budgets before a real factor: inverse numeric bytes 460382800 and factor upper bound 484807360 fit a proposed 512MiB factor cap but do not include full graph, basis, temporaries or process runtime. No overlap inverse exists yet. Preserve all physical CVs/F-C roles, original full C and frozen same-Q best NS direction. Changing F maps requires an explicit consistent Schur/block-preconditioner treatment or a properly variable outer method; never claim the old S identity. Complete selected-J and original-NS total-direction 0.1/original norm(R), temperature/state/descent/spatial/steady/flux and regional-physics/experimental gates remain required; scientific UI stays parallel.',
    goalStatus='active', wholeGoalComplete=False, recorderSHA256=sha(__file__))
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], manifest=identity(target), researchRows=16)))
