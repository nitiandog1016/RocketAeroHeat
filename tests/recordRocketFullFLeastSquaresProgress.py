# -*- coding: utf-8 -*-
"""记录跨组诊断、独立LSQR筛查、MATLAB停止条件修正及工作台实际证据。"""
from pathlib import Path
from datetime import datetime, timezone
import json
import hashlib
import os
import subprocess

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
TARGET = ROOT/'docs/validation/v83/progress_full_F_least_squares.json'
assert not TARGET.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(8*1024**2), b''):
            h.update(part)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def receipt(path):
    return dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path), bytes=path.stat().st_size)


prior_file = ROOT/'docs/validation/v83/progress_restricted_overlap.json'
assert sha(prior_file) == 'a466b4b524a25fe6352ea617b74fe6c064ee65fd7b182a650f8a0fb2463689f7'
prior = read(prior_file)
processes = {}
for name in ('fullFlsqrunit01', 'fullFlsqrreal01', 'fullFlsqrunit02', 'fullFlsqrreal02', 'fullFlsqrui01'):
    file = VAL/(name+'_process.json')
    p = read(file)
    assert p['status'] == 'terminal' and p['exitCode'] == 0 and all(p['sourceStableAfter'].values())
    processes[name] = dict(exitCode=0, sourceStable=True, engineProcessIDs=p['engineProcessIDs'], receipt=receipt(file))
audits = []
for name in ('fullFlsqrunit01_independent_full_F_least_squares.json', 'fullFlsqrreal01_independent_full_F_least_squares.json',
    'fullFlsqrunit02_independent_full_F_least_squares.json', 'fullFlsqrreal02_independent_full_F_least_squares.json',
    'workbench28/independent_full_F_least_squares_study_view.json', 'workbench28/node_wall_velocity_summary.json'):
    file = VAL/name
    a = read(file)
    assert a['status'].startswith('passed-')
    for path, digest in a.get('inputSHA256', {}).items():
        path = Path(path)
        assert sha(path if path.is_absolute() else ROOT/path) == digest
    audits.append(dict(report=a, receipt=receipt(file)))
screening = []
for folder, summary in [('fullFcolorscreen01', 'full_F_color_screen_summary.json'), ('fullFlsqrscreen01', 'full_F_lsqr_screen_summary.json')]:
    file = VAL/folder/summary
    d = read(file)
    assert d['status'].startswith('tested-') and d['newActualNSEvaluations'] == d['newProductionSolvedAreaM2'] == 0
    for path, digest in d['inputSHA256'].items():
        assert sha(ROOT/path) == digest
    screening.append(dict(report=d, receipt=receipt(file)))
current_file = VAL/'fullFlsqrreal02/full_F_least_squares_summary.json'
r = read(current_file)
assert r['newActualNSEvaluations'] == r['newFlowAcceptances'] == r['newProductionSolvedAreaM2'] == 0
assert not r['fullFReport']['fullFInternalTargetMet']
for name, digest in r['sourceIdentity'].items():
    assert sha(ROOT/'internalFlow/cutCell'/(name+'.m')) == digest
old = read(VAL/'fullFlsqrreal01/full_F_least_squares_summary.json')
assert old['fullFReport']['nativeFlag'] == 0 and old['fullFReport']['nativeIterations'] == 1 and not old['fullFReport']['fullFInternalTargetMet']
assert old['sourceIdentity']['solveRocketFullFLeastSquares'] == sha(VAL/'source_epochs/pre_lsqr_normal_stop_solveRocketFullFLeastSquares.m')
failure_file = VAL/'fullFlsqrreal01_audit_path_failure.json'
failure = read(failure_file)
assert failure['driverExitCode'] == 1
for s in failure['sources']:
    assert sha(ROOT/s['path']) == s['sha256']
production = prior['production']
assert sha(ROOT/production['latestEndpoint']['path']) == production['latestEndpoint']['sha256']
assert sha(ROOT/production['originalNSSource']['path']) == production['originalNSSource']['sha256']
head = subprocess.check_output(['git','rev-parse','--short','HEAD'], cwd=ROOT, text=True).strip()
assert head == prior['delivery']['gitHEAD'] == '42601c9'
sources = ['internalFlow/cutCell/'+name+'.m' for name in ('applyRocketPackedGraphTranspose','solveRocketSpatialFBlocksTranspose','solveRocketFullFLeastSquares')]
sources += ['tests/'+name for name in ('screenRocketFullFColorSweep.py','screenRocketFullFLeastSquares.py','testRocketFullFLeastSquares.m',
    'studyRocketFullFLeastSquares.m','auditRocketFullFLeastSquares.py','buildRocketFullFLeastSquaresSummary.py',
    'auditRocketFullFLeastSquaresWorkbench.py','recordRocketFullFLeastSquaresProgress.py','testRocketWallVelocityWorkbench.js')]
sources += ['web/rocket-wall-velocity.js']
out = dict(version=83, status='tested-complete-F-global-LSQR-improved-but-unconverged-with-native-stop-repair-and-workbench',
    recordedUTC=datetime.now(timezone.utc).isoformat(), processID=os.getpid(), priorProgress=receipt(prior_file),
    previousGoalTurnClassification='progress-restricted-overlap-closure-and-independent-output-weight-screening',
    currentGoalTurnClassification='progress-connected-F-diagnosis-dense-adjoint-validation-real-LSQR-and-workbench',
    processes=processes, independentChecks=audits, independentPrototypeAndColorScreening=screening,
    retainedNormalEquationEarlyStop=dict(summary=old, receipt=receipt(VAL/'fullFlsqrreal01/full_F_least_squares_summary.json'),
        nativeFlagZeroNotUsedAsFullFConvergence=True, oldEvidenceOverwritten=False),
    retainedAuditPathFailure=failure, retainedAuditPathFailureReceipt=receipt(failure_file),
    currentRealSummary=r, currentRealSummaryReceipt=receipt(current_file),
    codeSourceIdentity=[receipt(ROOT/path) for path in sources],
    methodDocumentation=receipt(ROOT/'docs/FULL_F_LEAST_SQUARES_V83.md'),
    goalAndReadinessSnapshots=[receipt(ROOT/path) for path in ('docs/PROJECT_GOALS.md','docs/GITHUB_STAGE_READINESS.md')],
    userInterface=dict(workbench='outputs/v83_validation/workbench28/index.html', researchRows=18,
        oldSeventeenRowsAndSixFacetFilesExactlyPreserved=True, LSQRArnoldiMetricsExplicitlyNotApplicable=True,
        browserVisualOrWebGLRuntimeVerified=False), production=production, delivery=prior['delivery'],
    newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0,
    nextAction='Do not repeat the failed symmetric color sweep or merely add LSQR steps. Exact forward/transpose and fixed M_F transpose now exist and are independently checked. Examine the scaled complete A_FF*M_F image/column conditioning or slow global directions using actual operator images and explicit memory budgets, before selecting the next global correction. Preserve all original F/C unknowns, fixed-C diagnostic identity and original NS total-direction gate; no automatic old-Schur replacement. Scientific UI and actual region pressure/heat validation remain highest priorities, not trajectory/history expansion.',
    goalStatus='active', wholeGoalComplete=False, recorderSHA256=sha(__file__))
TARGET.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], manifest=str(TARGET), sha256=sha(TARGET))))
