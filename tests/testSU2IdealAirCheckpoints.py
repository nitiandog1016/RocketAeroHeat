"""模块用途。核对新增阶段输出、稀疏保存及冷重启动的合同，保留每项真实终态。"""
from pathlib import Path
import copy
import json
import os
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from validateIdealAir import run, save, sha, relative
from su2ThermalRuntime import resolve, read
from prepareSU2IdealAirRestart import prepare
from su2IdealAirCheckpoints import output_plan, inventory
from auditSU2IdealAirState import audit, load_state

OUT = ROOT/'outputs/v156_validation'
SOURCE = ROOT/'outputs/v155_validation/cases01/small'
SNAPSHOT = ROOT/'outputs/v155_validation/executed_sources'


def reject_contracts(folder):
    """模块一。配置、累计来源、预算和状态损坏均须在原生进程创建前拒绝。"""
    rejected = []
    for name, arguments in [('zero', (38, 6, 8, 6, 0, 32*1024**2)),
                            ('capacity', (335726, 9878, 19752, 600, 1, 32*1024**2))]:
        try:
            output_plan(*arguments)
        except AssertionError as error:
            rejected.append(dict(name=name, error=str(error)))
    base = OUT/'contracts01'
    base.mkdir()
    for name in ('wrongCumulative', 'wrongVelocity', 'changedState', 'missingStage', 'partialPreparation'):
        target = base/name
        target.mkdir()
        for path in folder.iterdir():
            if path.is_file() and path.name != 'engine_request.json':
                shutil.copyfile(str(path), str(target/path.name))
        manifest = read(target/'case_manifest.json')
        text = (target/'case.cfg').read_text()
        if name == 'wrongCumulative':
            manifest['initialization']['priorSuccessfulUpdates'] = 88
            save(target/'ideal_air_restart.json', manifest['initialization'])
            manifest['idealAir']['initializationSHA256'] = sha(target/'ideal_air_restart.json')
            prepared = read(target/'ideal_air_restart_preparation.json')
            prepared['initializationSHA256'] = manifest['idealAir']['initializationSHA256']
            save(target/'ideal_air_restart_preparation.json', prepared)
        elif name == 'wrongVelocity':
            manifest['expectedVelocity'][0] *= .5
        elif name == 'changedState':
            with (target/'initial_flow.csv').open('a') as stream:
                stream.write('\n')
        elif name == 'missingStage':
            text = text.replace(', FLOW_PHASE', '')
        else:
            save(target/'ideal_air_restart_preparation.json', dict(status='preparing', exitCode=None))
        (target/'case.cfg').write_text(text)
        manifest['idealAir']['configurationSHA256'] = sha(target/'case.cfg')
        save(target/'ideal_air_contract.json', manifest['idealAir'])
        manifest['idealAirContractSHA256'] = sha(target/'ideal_air_contract.json')
        save(target/'case_manifest.json', manifest)
        try:
            resolve(target, text)
        except AssertionError as error:
            rejected.append(dict(name=name, error=str(error)))
        assert not (target/'engine_process.json').exists()
    assert len(rejected) == 7, rejected
    # 读取历史已验证接口，不重复原生求解。
    for relative_path in ('outputs/v154_validation/cases01/laminar', 'outputs/v154_validation/cases01/sst',
                          'outputs/v145_validation/fixture02/case'):
        old = ROOT/relative_path
        resolve(old, (old/'case.cfg').read_text())
    save(base/'report.json', dict(passed=True, rejected=rejected, legacyContracts=3))


def native(folder, name, sources):
    """模块二。通过正常监督器运行，所有原生任务均只执行一次。"""
    save(folder/'engine_request.json', dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
        threads=2, maximumSeconds=60, sourceFiles=[path.relative_to(ROOT).as_posix() for path in sources]))
    run(OUT/'supervision', name, [sys.executable, '-X', 'utf8', str(ROOT/'scripts/runEmbeddedSU2.py'),
                                '--request', str(folder/'engine_request.json')], 70)
    receipt = read(folder/'engine_process.json')
    assert receipt['exitCode'] == 0 and receipt['checkpoints']['normalFinalPairPresent']
    assert not receipt['checkpoints']['rejected'] and all(receipt['inputsStable'].values())
    return receipt


def main():
    """模块三。完整与稀疏共有帧逐字段相等，独立热状态及分段来源分项报告。"""
    cases = OUT/'cases01'
    assert cases.is_dir() and (OUT/'supervision').is_dir()
    sources = [ROOT/'scripts'/name for name in ('su2IdealAirInterface.py', 'su2ThermalRuntime.py',
        'prepareSU2IdealAirRestart.py', 'su2IdealAirCheckpoints.py', 'auditSU2IdealAirState.py', 'runEmbeddedSU2.py')]
    sources += [Path(__file__).resolve(), ROOT/'docs/validation/v156/PLAN.md']
    locks = {str(path): sha(path) for path in sources}
    for path in sources:
        target = OUT/'executed_sources02'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    save(OUT/'source_lock02.json', locks)
    full, sparse, restart = cases/'full', cases/'sparse', cases/'restart'
    # 纠错后复用原控制的实际1—5帧，原首帧失败和已用预算继续保留。
    full_receipt = read(full/'engine_process.json')
    assert full_receipt['status'] == 'terminal' and full_receipt['exitCode'] == 0
    assert [row['localIteration'] for row in full_receipt['checkpoints']['complete']] == list(range(1, 6))
    full_report = audit(full)
    prepare(SOURCE, sparse, 6, 2, 32*1024**2, snapshot=SNAPSHOT)
    sparse_receipt = native(sparse, 'sparse', sources)
    assert [row['localIteration'] for row in sparse_receipt['checkpoints']['complete']] == [0, 2, 4, 5]
    reports, equality = [], {}
    for step in (0, 2, 4, 5):
        if step > 0:
            a, fields, nd = load_state(full/('restart_flow_%06d.csv' % step), True, captured=True)
            b, _, _ = load_state(sparse/('restart_flow_%06d.csv' % step), True, captured=True)
            equality[str(step)] = relative(a, b)
            assert equality[str(step)] < 2e-13
        reports.append(audit(sparse, sparse/('audit_%06d' % step), checkpoint=step))
    first, first_fields, _ = load_state(sparse/'restart_flow_000000.csv', True, captured=True)
    initial = np.genfromtxt(str(full/'initial_flow.csv'), delimiter=',', names=True)
    first_error = relative(first_fields['Flow_Primitive_K'], initial['Turb_Kin_Energy'])
    assert first_error < 2e-13
    cross_errors, separation = [], []
    for step in range(2, 6):
        _, fields, _ = load_state(full/('restart_flow_%06d.csv' % step), True, captured=True)
        _, previous, _ = load_state(full/('restart_flow_%06d.csv' % (step-1)), True, captured=True)
        cross_errors.append(relative(fields['Flow_Primitive_K'], previous['Turb_Kin_Energy']))
        separation.append(float(max(abs(fields['Flow_Primitive_K']-fields['Turb_Kin_Energy']))))
    assert max(cross_errors) < 2e-13 and max(separation) > 1e-5
    prepare(sparse, restart, 3, 2, 32*1024**2, checkpoint=2)
    restart_receipt = native(restart, 'restart', sources)
    restart_report = audit(restart)
    assert restart_report['successfulUpdates'] == 3 and restart_report['seriesSuccessfulUpdates'] == 6
    assert restart_receipt['checkpoints']['lastComplete']['localIteration'] == 2
    _, restarted, _ = load_state(restart/'restart_flow_000000.csv', True, captured=True)
    _, source_fields, _ = load_state(sparse/'restart_flow_000002.csv', True, captured=True)
    restart_k = relative(restarted['Flow_Primitive_K'], source_fields['Turb_Kin_Energy'])
    assert restart_k < 2e-13
    end, _, _ = load_state(restart/'restart_flow_000002.csv', True, captured=True)
    baseline, _, _ = load_state(full/'restart_flow_000005.csv', True, captured=True)
    stable = {name: sha(Path(name)) == digest for name, digest in locks.items()}
    assert all(stable.values())
    result = dict(passed=True, full=full_report, sparse=reports, restart=restart_report,
        commonFrameErrors=equality, firstPrimitiveKError=first_error, previousPrimitiveKErrors=cross_errors,
        stageKMaximumDifferences=separation, restartFirstPrimitiveKError=restart_k,
        coldRestartVsUninterruptedDifference=relative(end, baseline),
        coldRestartCarriesAdaptiveCFL=False, sourceSHA256=locks, inputsStable=stable,
        nativeRuns=3, realFullDomainNativeRuns=0, physicsQualified=False, newQualifiedArea=0)
    save(OUT/'native_report.json', result)
    print(json.dumps(dict(passed=True, nativeRuns=3, firstPrimitiveKError=first_error,
                         commonFrameErrors=equality, coldRestartVsUninterruptedDifference=result['coldRestartVsUninterruptedDifference'])))


if __name__ == '__main__':
    main()
