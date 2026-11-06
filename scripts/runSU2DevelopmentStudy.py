"""模块用途。固定v164同包冷重启合同，先小例，随后唯一有总预算的整箭发展。"""
import argparse
import copy
import json
from pathlib import Path
import sys
import numpy as np
from auditSU2IdealAirState import audit, load_state
from prepareSU2IdealAirRestart import prepare, source_state, admissible_checkpoint_identity
from su2IdealAirInterface import ADMISSIBLE_FLAVOR, STEP_FLAVOR, settings
from su2ThermalRuntime import ROOT, read, resolve, sha
from validateIdealAir import run, save

BASE = ROOT/'outputs/v166_validation'
SOURCE = ROOT/'outputs/v164_validation'


def launch(folder, seconds, parent_seconds):
    """模块一。原生/父监督预算分别绑定，已有回执或目录禁止重复启动。"""
    names = ['scripts/runSU2DevelopmentStudy.py', 'scripts/prepareSU2IdealAirRestart.py',
        'scripts/su2IdealAirInterface.py', 'scripts/su2IdealAirCheckpoints.py',
        'scripts/auditSU2IdealAirState.py', 'scripts/runEmbeddedSU2.py',
        'scripts/su2ThermalRuntime.py', 'scripts/validateIdealAir.py', 'docs/validation/v166/PLAN.md']
    for name in names:
        target = folder/'executed_sources'/name
        target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes((ROOT/name).read_bytes())
    request = dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
        maximumSeconds=seconds, threads=2, sourceFiles=names, jobToken='v166-'+folder.name)
    if seconds == 900: request['resourceProfile'] = 'v166-single-development'
    save(folder/'engine_request.json', request)
    run(folder, 'parent', [sys.executable, str(ROOT/'scripts/runEmbeddedSU2.py'),
        '--request', str(folder/'engine_request.json')], parent_seconds)


def small():
    """模块二。真实已审核末帧接续；反例不修改历史源，原生仅一个三更新任务。"""
    folder, source = BASE/'small01', SOURCE/'small01'
    manifest = prepare(source, folder, 3, 3, 32*1024**2, checkpoint=3, runtime_flavor=ADMISSIBLE_FLAVOR)
    _, runtime, flavor, _, _ = resolve(folder, (folder/'case.cfg').read_text('utf8'))
    assert flavor == ADMISSIBLE_FLAVOR and runtime['samePhysicsAndNumerics'] is False
    init = manifest['initialization']
    assert init['samePhysicsAndNumerics'] and init['priorSuccessfulUpdates'] == 4
    assert not init['runTransition']['uninterruptedTrajectoryEquivalent']
    old_manifest, receipt = read(source/'case_manifest.json'), read(source/'engine_process.json')
    rejected = []
    for name in ('failed-source', 'running-source', 'unstable-source', 'incomplete-schedule',
                 'wrong-checkpoint', 'wrong-stage', 'changed-binary', 'changed-runtime',
                 'changed-surface', 'step-source'):
        m, r, step = copy.deepcopy(old_manifest), copy.deepcopy(receipt), 3
        if name == 'failed-source': r['exitCode'] = 1
        elif name == 'running-source': r['status'] = 'running'
        elif name == 'unstable-source': r['inputsStable'][next(iter(r['inputsStable']))] = False
        elif name == 'incomplete-schedule': r['checkpoints']['normalScheduleComplete'] = False
        elif name == 'wrong-checkpoint': step = 2
        elif name == 'wrong-stage': r['checkpoints']['lastComplete']['stage'] = 'previous-k'
        elif name == 'changed-binary': r['inputSHA256'][init['runTransition']['executablePath']] = '0'*64
        elif name == 'changed-runtime': r['inputSHA256'][init['runTransition']['runtimeManifestPath']] = '0'*64
        elif name == 'changed-surface': r['checkpoints']['lastComplete']['surfaceSHA256'] = '0'*64
        else: m['requiredRuntimeFlavor'] = r['runtimeFlavor'] = STEP_FLAVOR
        try: admissible_checkpoint_identity(source, m, r, step)
        except (AssertionError, KeyError): rejected.append(name)
    assert len(rejected) == 10
    try: source_state(ROOT/'outputs/v165_validation/small01', 3, None)
    except AssertionError as error:
        assert '不支持该保存态来源包' in str(error); rejected.append('actual-v165-source')
    assert len(rejected) == 11
    save(folder/'contract_checks.json', dict(passed=True, rejected=rejected,
        operatorUnchangedFromSource=True, runtimeHistoricallyChangedFromObserved=True,
        uninterruptedEquivalent=False, sourceSHA256=init['sourceSHA256']))
    launch(folder, 30, 120)
    receipt = read(folder/'engine_process.json')
    assert receipt['checkpoints']['normalScheduleComplete']
    assert [e['localIteration'] for e in receipt['checkpoints']['complete']] == [0, 2]
    initial = np.loadtxt(folder/'initial_flow.csv', delimiter=',', skiprows=1)
    old, _, _ = load_state(source/'restart_flow_000003.csv', True, captured=True)
    assert np.array_equal(initial, old[:, :initial.shape[1]])
    thermal0 = audit(folder, folder/'audit00', checkpoint=0)
    thermal = audit(folder)
    assert thermal0['seriesSuccessfulUpdates'] == 5 and thermal['seriesSuccessfulUpdates'] == 7
    save(BASE/'small_report.json', dict(passed=True, initialAllColumnsExact=True, rejections=rejected,
        strictThermalErrors={'0': thermal0['errors'], '2': thermal['errors']},
        priorUpdates=4, localUpdates=3, cumulativeUpdates=7, noCompilation=True,
        operatorUnchanged=True, uninterruptedTrajectoryEquivalent=False, physicsQualified=False, newQualifiedArea=0))
    print(json.dumps(dict(smallPassed=True, rejected=len(rejected), thermal=thermal['errors'])))


def real(recovery=False):
    """模块三。源累计200次与本段200上限分别记账；没有自动续算分支。"""
    assert read(BASE/'small_report.json')['passed']
    source, folder = SOURCE/'real01', BASE/('real02' if recovery else 'real01')
    if recovery:
        failed = BASE/'real01'
        assert read(failed/'parent_process.json')['status'] == 'terminal'
        assert read(failed/'parent_process.json')['exitCode'] == 1
        assert not (failed/'engine_process.json').exists() and not list(failed.glob('restart_flow*.csv'))
        assert '1 <= timeout <= 720' in (failed/'parent.log').read_text('utf8')
    m = prepare(source, folder, 200, 100, 1536*1024**2, checkpoint=49, runtime_flavor=ADMISSIBLE_FLAVOR)
    assert (m['nodes'], m['cells']) == (335726, 705490)
    assert m['initialization']['priorSuccessfulUpdates'] == 200
    assert m['initialization']['samePhysicsAndNumerics']
    assert sha(folder/'mesh.su2') == sha(source/'mesh.su2')
    old, new = [settings((p/'case.cfg').read_text('utf8')) for p in (source, folder)]
    changed = {k: {'source': old[k], 'target': v} for k, v in new.items() if old.get(k) != v}
    assert changed == {'ITER': {'source': '50', 'target': '200'}, 'OUTPUT_WRT_FREQ': {'source': '50', 'target': '100'}}
    assert new['HISTORY_OUTPUT'] == '(ITER, RMS_RES, AERO_COEFF, HEAT, LINSOL)'
    if recovery:
        assert all(sha(folder/n) == sha(BASE/'real01'/n) for n in ('case.cfg', 'mesh.su2', 'initial_flow.csv'))
        save(folder/'prelaunch_recovery.json', dict(originalDirectory=str(BASE/'real01'),
            originalParentSHA256=sha(BASE/'real01/parent_process.json'),
            originalLogSHA256=sha(BASE/'real01/parent.log'), priorNativeLaunches=0,
            inputFieldsExact=True, criteriaChanged=False, uniqueNativeDevelopmentTasks=1))
    save(folder/'development_preflight.json', dict(passed=True, planSHA256=sha(ROOT/'docs/validation/v166/PLAN.md'),
        configurationChanges=changed, expectedLocalCheckpoints=[0, 100, 199],
        sourceSuccessfulUpdates=200, maximumAdditionalUpdates=200, maximumSeriesSuccessfulUpdates=400,
        outputPlan=m['idealAir']['checkpointContract'], sameOperatorAndRuntime=True,
        uninterruptedTrajectoryEquivalent=False, engineCompilation=False, automaticAdditionalTasks=0))
    launch(folder, 900, 1020)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=('small', 'real', 'real-recovery')); args = parser.parse_args()
    small() if args.case == 'small' else real(args.case == 'real-recovery')
