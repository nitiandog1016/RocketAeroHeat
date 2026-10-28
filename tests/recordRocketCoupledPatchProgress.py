# -*- coding: utf-8 -*-
"""登记行尺度失败对照、五量局部块真实试算及十一行工作台的实际证据。"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, hashlib, subprocess

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
target = ROOT / 'docs/validation/v83/progress_coupled_patches.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def identity(path):
    path = Path(path)
    return dict(path=str(path.relative_to(ROOT)), sha256=sha(path), bytes=path.stat().st_size)


names = ['rowmomentumunit01', 'rowmomentumreal01', 'coupledpatchunit01', 'coupledpatchunit02',
         'coupledpatchunit03', 'coupledpatchreal01', 'coupledpatchui01']
processes = {}
for name in names:
    file = VAL / (name+'_process.json')
    p = read(file)
    expected = 1 if name == 'coupledpatchunit02' else 0
    assert p['status'] == 'terminal' and p['exitCode'] == expected and all(p['sourceStableAfter'].values())
    processes[name] = dict(exitCode=p['exitCode'], sourceStable=True, engineProcessIDs=p['engineProcessIDs'], receipt=identity(file))

row_unit = read(VAL / 'rowmomentumunit01_independent_row_momentum.json')
row_real = read(VAL / 'rowmomentumreal01_independent_row_momentum.json')
coupled_units = [read(VAL / (name+'_independent_coupled_patches.json')) for name in ('coupledpatchunit01', 'coupledpatchunit03')]
coupled_real = read(VAL / 'coupledpatchreal01_independent_coupled_patches.json')
assert all(r['status'].startswith('passed-') for r in [row_unit, row_real, coupled_real]+coupled_units)
assert not row_real['originalAResidualImproved'] and not row_real['originalUnitTargetMet'] and not row_real['rowScaledTargetMet']
assert not coupled_real['internalCTargetMet'] and not coupled_real['independentAllRealArnoldiColumnAuditClaimed']

row_summary_file = VAL / 'rowmomentumreal01/row_momentum_real_summary.json'
coupled_summary_file = VAL / 'coupledpatchreal01/coupled_patch_real_summary.json'
assert row_summary_file.stat().st_size < 128*1024 and coupled_summary_file.stat().st_size < 128*1024
row_summary = read(row_summary_file)
coupled_summary = read(coupled_summary_file)
assert not row_summary['physicalSolution'] and not coupled_summary['physicalSolution']
coverage = read(VAL / 'coupledpatchreal01/patch_coverage_and_remaining.json')
assert coverage['completeElevenLayerColumnsInOnePatch'] == 0 and not coverage['observedResidualReductionProvesGlobalOrUniqueFailureCause']
for audit in (row_real, coupled_real, coverage):
    for path, digest in audit['inputSHA256'].items():
        assert sha(ROOT/path) == digest

ui = read(VAL / 'workbench23/independent_coupled_patch_study_view.json')
node = read(VAL / 'workbench23/node_wall_velocity_summary.json')
assert ui['status'].startswith('passed-') and ui['researchRows'] == 11
assert node['status'].startswith('passed-') and node['CSVrows'] == 592560 and node['sectionSamples'] == 6624 and node['sameHostEvents'] == 15
assert not ui['browserVisualOrWebGLRuntimeVerified'] and not node['browserVisualOrWebGLRuntimeVerified']
for item in ui['sourceFacetPayload']:
    assert sha(ROOT/item['path']) == item['sha256']

old = read(ROOT / 'docs/validation/v83/progress_momentum_lines.json')
production = old['production']
for key in ('latestEndpoint', 'originalNSSource'):
    assert sha(ROOT/production[key]['path']) == production[key]['sha256']
assert subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, text=True).strip() == '42601c9'

backend = ['rocketMomentumRowScales.m', 'solveRocketRowScaledMomentumKrylov.m', 'compactRocketMomentumKrylovReport.m',
           'buildRocketStrongPatches.m', 'certifyRocketLocalInverse.m', 'factorRocketCoupledPatches.m',
           'solveRocketCoupledPatches.m', 'solveRocketCoupledPatchKrylov.m', 'compactRocketCoupledKrylovReport.m']
tests = ['testRocketRowScaledMomentum.m', 'studyRocketRowScaledMomentum.m', 'auditRocketRowMomentumUnit.py',
         'auditRocketRowMomentumReal.py', 'testRocketCoupledPatches.m', 'studyRocketCoupledPatches.m',
         'auditRocketCoupledPatchUnit.py', 'auditRocketCoupledPatchReal.py', 'diagnoseRocketCoupledPatchCoverage.py',
         'exportRocketCoupledPatchDiagnostics.py', 'auditRocketCoupledPatchWorkbench.py']
epochs = ['pre_coupled_patch_bucket_factorRocketCoupledPatches.m', 'pre_coupled_patch_fallback_testRocketCoupledPatches.m',
          'pre_coupled_single_group_factorRocketCoupledPatches.m']
paths = [ROOT/'internalFlow/cutCell'/name for name in backend] + [ROOT/'tests'/name for name in tests]

out = dict(status='tested-equivalent-row-scale-and-bounded-full-five-field-local-patches-with-unresolved-real-residuals',
    version=83, recordedUTC=datetime.now(timezone.utc).isoformat(), processID=os.getpid(),
    previousVisibleTurnClassification='no-progress-status-answer-to-user-version-question',
    revalidatedNextSafeAction='Checked current goals, sources and absent live jobs, then completed one bounded actual real full-C coupled-patch run.',
    currentGoalTurnClassification='progress-new-real-numerical-evidence-independent-audits-and-updated-workbench',
    processes=processes, rowScaleSmallAudit=row_unit, rowScaleRealAudit=row_real,
    coupledPatchSmallAudits=coupled_units, coupledPatchRealAudit=coupled_real,
    rowScaleSummary=row_summary, coupledPatchSummary=coupled_summary,
    compactSummaryReceipts=[identity(row_summary_file), identity(coupled_summary_file)],
    retainedFailure=dict(case='coupledpatchunit02', actualExitCode=1,
        cause='MATLAB repelem one-group row-vector output implicitly expanded against a column, causing a size mismatch.',
        repair='Convert the expanded group-start vector to an explicit column before assigning member positions.',
        repairedActualTest='coupledpatchunit03', originalSourceEpochs=[identity(VAL/'source_epochs'/name) for name in epochs]),
    patchMembershipAndRemainingDiagnostic=coverage,
    methodDocumentation=identity(ROOT/'docs/COUPLED_PATCH_PRECONDITIONER_V83.md'),
    codeSourceIdentity=[identity(path) for path in paths],
    userInterface=ui, nodeDataVerification=node, production=production, delivery=old['delivery'],
    workAccounting=dict(rowScaleSolverAApplications=33, rowScaleOriginalUnitDiagnosticAApplications=1,
        rowScaleFixedLocalMapApplications=32, rowScaleNewFactorConstructions=0, rowScaleNewSchurAssemblies=0,
        coupledPatchNumericalCApplications=35, coupledPatchNumericalFixedLocalMapApplications=33,
        coupledPatchIndependentAuditFullCApplications=2, coupledPatchIndependentAuditFixedLocalMapApplications=1,
        membershipDiagnosticNewCApplications=0, membershipDiagnosticNewLocalMapApplications=0,
        newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0),
    nextAction='Do not repeat equivalent row scaling, old global momentum ILU, stopped scalar ILUTP, old Schur assembly, unchanged normal-line scans or old 90/210/792 image projections, or blindly extend the 32-step budget. Reuse complete C and the audited fixed full-five-field local map; assess a bounded strong-connection multilevel or overlapping correction on dense small-domain references with explicit storage and restriction/prolongation checks before any real run. All physical CVs and all C couplings stay active. A successful internal solve must still be tested on the complete selected J and total original NS direction with the original norm(R), 0.1 target and all state/descent/spatial/steady/flux gates. Physical regional applicability, sensitivities, public evidence and scientific UI remain required.',
    goalStatus='active', wholeGoalComplete=False, recorderSHA256=sha(__file__))
target.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(status=out['status'], manifest=identity(target), researchRows=11)))
