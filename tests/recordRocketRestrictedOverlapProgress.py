# -*- coding: utf-8 -*-
"""登记受限重叠真实试算、独立核对、只读输出筛查和工作台；不登记物理交付。"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
TARGET = ROOT / 'docs/validation/v83/progress_restricted_overlap.json'
assert not TARGET.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def receipt(path):
    return dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path), bytes=path.stat().st_size)


prior_file = ROOT / 'docs/validation/v83/progress_full_F.json'
assert sha(prior_file) == 'f920e244a2a36d6d36ecdbc2b88a5ae44e08532ade85eed17b9a3d9c5ac1bceb'
prior = read(prior_file)
processes = {}
for name in ('overlapFunit01', 'overlapFreal01', 'overlapFui01'):
    file = VAL / (name + '_process.json')
    p = read(file)
    assert p['status'] == 'terminal' and p['exitCode'] == 0 and all(p['sourceStableAfter'].values())
    assert p['engineProcessIDs']
    processes[name] = dict(exitCode=0, sourceStable=True,
        engineProcessIDs=p['engineProcessIDs'], receipt=receipt(file))

unit_file = VAL / 'overlapFunit01_independent_restricted_overlap.json'
real_file = VAL / 'overlapFreal01_independent_restricted_overlap.json'
ui_file = VAL / 'workbench27/independent_restricted_overlap_study_view.json'
node_file = VAL / 'workbench27/node_wall_velocity_summary.json'
checks = {}
for name, file in [('small', unit_file), ('real', real_file), ('workbench', ui_file), ('node', node_file)]:
    a = read(file)
    assert a['status'].startswith('passed-')
    for path, digest in a.get('inputSHA256', {}).items():
        path = Path(path)
        assert sha(path if path.is_absolute() else ROOT / path) == digest
    checks[name] = dict(report=a, receipt=receipt(file))

summary_file = VAL / 'overlapFreal01/restricted_overlap_summary.json'
r = read(summary_file)
assert r['newActualNSEvaluations'] == r['newFlowAcceptances'] == r['newProductionSolvedAreaM2'] == 0
assert r['innerRunComplete'] and r['parentDiagnosticComplete'] and not r['fullFInternalTargetMet']
for name, digest in r['sourceIdentity'].items():
    assert sha(ROOT / 'internalFlow/cutCell' / (name + '.m')) == digest
screen_file = VAL / 'overlap_output_diagnostic01.json'
screen = read(screen_file)
assert screen['status'].startswith('tested-') and len(screen['smallDenseReferences']) == 2
assert len(screen['realActionComparisons']) == 3 and screen['fullFActionEvaluations'] == 6
assert screen['newActualNSEvaluations'] == screen['newKrylovSolves'] == screen['newProductionSolvedAreaM2'] == 0
for path, digest in screen['inputSHA256'].items():
    assert sha(ROOT / path) == digest
assert sha(ROOT / 'tests/diagnoseRocketOverlapOutput.py') == screen['scriptSHA256']
failure_file = VAL / 'overlap_output_diagnostic01_reference_failure.json'
failure = read(failure_file)
assert failure['driverExitCode'] == 1 and not failure['realScreeningStarted']
assert sha(ROOT / failure['sourceFile']) == failure['sourceSHA256']

production = prior['production']
assert sha(ROOT / production['latestEndpoint']['path']) == production['latestEndpoint']['sha256']
assert sha(ROOT / production['originalNSSource']['path']) == production['originalNSSource']['sha256']
head = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip()
assert head == prior['delivery']['gitHEAD'] == '42601c9'
files = ['internalFlow/cutCell/' + name + '.m' for name in
    ('factorRocketRestrictedOverlap', 'solveRocketRestrictedOverlap', 'solveRocketOverlapFullF')]
files += ['tests/' + name for name in ('testRocketRestrictedOverlap.m', 'studyRocketRestrictedOverlap.m',
    'auditRocketRestrictedOverlapUnit.py', 'auditRocketRestrictedOverlapReal.py',
    'buildRocketRestrictedOverlapSummary.py', 'auditRocketRestrictedOverlapWorkbench.py',
    'diagnoseRocketOverlapOutput.py', 'recordRocketRestrictedOverlapProgress.py')]
out = dict(version=83, status='tested-restricted-overlap-real-run-unimproved-with-read-only-output-screen-and-workbench-not-delivered',
    recordedUTC=datetime.now(timezone.utc).isoformat(), processID=os.getpid(),
    priorProgress=receipt(prior_file), previousGoalTurnClassification='no-progress-status-answer-only',
    previousGoalTurnRevalidation='Read current goals, actual terminal receipts, independent audits and real overlap result; no active numerical process required waiting.',
    currentGoalTurnClassification='progress-persisted-tested-overlap-added-source-bound-workbench-and-rejected-unpromising-output-weighting',
    processes=processes, independentChecks=checks, compactNumericalSummary=r,
    compactSummaryReceipt=receipt(summary_file), outputScreening=screen, outputScreeningReceipt=receipt(screen_file),
    retainedSmallReferenceFailure=failure, retainedSmallReferenceFailureReceipt=receipt(failure_file),
    codeSourceIdentity=[receipt(ROOT / file) for file in files],
    methodDocumentation=receipt(ROOT / 'docs/RESTRICTED_OVERLAP_V83.md'),
    goalAndReadinessRecordSnapshots=[receipt(ROOT / file) for file in ('docs/PROJECT_GOALS.md', 'docs/GITHUB_STAGE_READINESS.md')],
    userInterface=dict(workbench='outputs/v83_validation/workbench27/index.html', researchRows=17,
        originalFacetFieldsAndHistoryByteIdenticalToPrevious=True, previousSixteenRowsExactlyPreserved=True,
        falseSourceAndQualificationRejectionsPassed=True, browserVisualOrWebGLRuntimeVerified=False),
    production=production, delivery=prior['delivery'],
    workAccounting=dict(restrictedOverlapFactorsConstructed=1, numericalStudyFullFApplications=35,
        numericalStudyRestrictedOverlapApplications=33, numericalStudyOriginalCDiagnosticApplications=1,
        separateReadOnlyOutputScreenFullFApplications=6, newFullFGraphExtractions=0,
        newApproximateSchurAssemblies=0, newSelectedJApplications=0, newActualNSEvaluations=0,
        newFlowAcceptances=0, newProductionSolvedAreaM2=0),
    nextAction='Do not repeat the tested 32-step local-F/RAS configurations, enlarge halo radius blindly or launch a weighted-output Krylov run based on this screening. Reuse the exact complete A_FF, saved best residual and original CV/five-field identity to diagnose dominant residual cross-group coupling and difficult directions, then justify a bounded global correction by actual complete operator images and resource evidence. Preserve the complete original C and best same-Q NS direction; any changed F mapping needs an explicitly consistent Schur/block formulation or suitable variable outer solve. Original selected-J/NS direction, positivity/temperature/descent/spatial/steady/boundary, regional physical validity and public validation gates remain necessary. UI and actual regional heat solutions stay highest priorities; no trajectory/history expansion.',
    goalStatus='active', wholeGoalComplete=False, recorderSHA256=sha(__file__))
TARGET.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], manifest=str(TARGET), sha256=sha(TARGET),
    qualifiedHeatFluxAreaPercent=production['qualifiedHeatFluxAreaPercent'], newQualifiedAreaM2=0)))
