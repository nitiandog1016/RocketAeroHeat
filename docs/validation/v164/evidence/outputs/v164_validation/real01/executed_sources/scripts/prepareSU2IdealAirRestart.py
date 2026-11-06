"""同模型初场复用与冷重启动，保留原始守恒态、准确旧源和累计更新来源。"""
from pathlib import Path
import argparse
import copy
import csv
import os
import time

import numpy as np

from idealAirProperties import IdealAir
from prepareSU2SSTCase import set_value
from runEmbeddedSU2 import atomic_write_json
from su2IdealAirInterface import (MODEL, PHASE_FLAVOR, OBSERVED_FLAVOR, OBSERVATION,
    ADMISSIBLE_FLAVOR, ADMISSIBILITY, settings)
from su2IdealAirCheckpoints import output_plan, STAGE
from su2SSTInterface import mesh, state
from su2ThermalRuntime import ROOT, read, sha

KIND = 'same-mesh-ideal-air-state-v1'
OUTPUT_KEYS = {'ITER', 'OUTPUT_WRT_FREQ', 'VOLUME_OUTPUT'}


def historical_locks(identities, snapshot):
    """模块一。已变更源码仅由事前保存的准确快照承接，旧输入身份始终逐字保持。"""
    locks, rebinding = {}, []
    for name, digest in identities.items():
        path = Path(name)
        target = path
        if not path.is_file() or sha(path) != digest:
            assert snapshot is not None, '历史源码变化，须指定准确执行源快照'
            relative = path.relative_to(ROOT)
            assert relative.parts[0] in ('scripts', 'tests', 'third_party'), '只允许准确源码快照重绑定'
            target = Path(snapshot).resolve()/relative
            target.relative_to(ROOT/'outputs')
            rebinding.append(dict(original=str(path), exactSnapshot=str(target), sha256=digest))
        assert target.is_file() and sha(target) == digest, str(target)
        locks[str(target)] = digest
    return locks, rebinding


def source_state(source, checkpoint, snapshot):
    """核对来源合同，返回原守恒列及选定帧以前的实际更新次数。"""
    source = Path(source).resolve()
    source.relative_to(ROOT/'outputs')
    manifest = read(source/'case_manifest.json')
    contract = read(source/'ideal_air_contract.json')
    assert manifest['idealAir'] == contract and contract['thermodynamicModel'] == MODEL and contract['turbulent']
    assert sha(source/'ideal_air_contract.json') == manifest['idealAirContractSHA256']
    assert sha(source/'case.cfg') == contract['configurationSHA256']
    assert sha(source/'mesh.su2') == contract['meshSHA256']
    locks = {str(source/name): sha(source/name) for name in
             ('case.cfg', 'mesh.su2', 'case_manifest.json', 'ideal_air_contract.json', 'sst_contract.json')}
    assert sha(source/'sst_contract.json') == manifest['sstContractSHA256']
    assert read(source/'sst_contract.json') == manifest['sst']
    rebinding = []
    if checkpoint is None:
        assert contract['initialization'] == 'general-mesh-ideal-air-SST-wall-distance-v1'
        prepared = read(source/'ideal_air_preparation.json')
        init = read(source/'ideal_air_initialization.json')
        assert prepared['status'] == 'terminal' and prepared['exitCode'] == 0 and all(prepared['inputsStable'].values())
        assert init == manifest['initialization'] and not init['physicsQualified']
        assert sha(source/'ideal_air_initialization.json') == contract['initializationSHA256'] == prepared['initializationSHA256']
        old_locks, rebinding = historical_locks(init['sourceSHA256'], snapshot)
        locks.update(old_locks)
        for name in ('ideal_air_initialization.json', 'ideal_air_preparation.json'):
            locks[str(source/name)] = sha(source/name)
        selected = source/'initial_flow.csv'
        assert sha(selected) == init['initialStateSHA256']
        data, columns, nd, _, _, _ = state(selected, True)
        prior = 0
    else:
        from auditSU2IdealAirState import load_state
        assert manifest['requiredRuntimeFlavor'] == PHASE_FLAVOR
        receipt_path = source/'engine_process.json'
        receipt = read(receipt_path)
        assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
        assert receipt['runtimeFlavor'] == PHASE_FLAVOR
        for name in ('case.cfg', 'mesh.su2', 'case_manifest.json', 'ideal_air_contract.json'):
            assert sha(source/name) == receipt['inputSHA256'][str(source/name)]
        entry = next(row for row in receipt['checkpoints']['complete'] if row['localIteration'] == checkpoint)
        selected = source/entry['volume']
        assert selected.name == 'restart_flow_%06d.csv' % checkpoint
        assert sha(selected) == entry['volumeSHA256'] == receipt['outputSHA256'][selected.name]
        data, fields, nd = load_state(selected, True, captured=True)
        columns = list(fields)[:2*nd+5]
        data = data[:, :len(columns)]
        prior = manifest['initialization']['priorSuccessfulUpdates']+checkpoint+1
        assert entry['seriesSuccessfulUpdates'] == prior and entry['stage'] == STAGE
        locks[str(receipt_path)] = sha(receipt_path)
    locks[str(selected)] = sha(selected)
    return manifest, contract, data, columns, nd, prior, locks, rebinding, selected


def prepare(source, folder, iterations, every, maximum_bytes, checkpoint=None, snapshot=None, exported_folder=False,
            linear_history=False, runtime_flavor=PHASE_FLAVOR):
    """模块二。独占新目录发布配置，初始场和失败记录均不会覆盖历史文件。"""
    assert runtime_flavor in (PHASE_FLAVOR, OBSERVED_FLAVOR, ADMISSIBLE_FLAVOR), '未声明的三维重启运行包'
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    if exported_folder:
        # 主工作台已导出网格/冻结快照。只接受尚未发布物性、与子目录原导出逐字相等的输入。
        source = Path(source).resolve()
        assert checkpoint is None and folder.is_dir()
        assert source.parent == folder and source.name == 'ideal_air_initial_state'
        for name in ('case.cfg', 'case_manifest.json'):
            assert sha(folder/name) == sha(source/'ideal_air_original_inputs'/name)
        assert sha(folder/'mesh.su2') == sha(source/'mesh.su2')
        for name in ('ideal_air_restart_preparation.json', 'ideal_air_contract.json', 'engine_process.json'):
            assert not (folder/name).exists(), '主任务已使用或部分发布'
    else:
        folder.mkdir(parents=True, exist_ok=False)
    marker = folder/'ideal_air_restart_preparation.json'
    record = dict(status='preparing', exitCode=None)
    atomic_write_json(marker, record)
    started = time.monotonic()
    try:
        source = Path(source).resolve()
        manifest, old_contract, data, columns, nd, prior, locks, rebinding, selected = source_state(source, checkpoint, snapshot)
        assert nd == 3, '当前相容重启动合同限三维SST'
        manifest = copy.deepcopy(manifest)
        core = IdealAir()
        points, _, wall = mesh(source/'mesh.su2')
        assert len(data) == len(points) == manifest['nodes']
        assert np.allclose(data[:, 1:4], points, rtol=2e-14, atol=1e-14)
        plan = output_plan(len(data), len(np.unique([v for face in wall for v in face[1:]])),
                           len(wall), iterations, every, maximum_bytes)
        text = (source/'case.cfg').read_text(encoding='utf8')
        for key, value in dict(ITER=str(iterations), OUTPUT_WRT_FREQ=str(every),
                VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, FLOW_PHASE)').items():
            text = set_value(text, key, value)
        if linear_history:
            assert checkpoint is not None and not exported_folder
            assert settings(text)['HISTORY_OUTPUT'] == '(ITER, RMS_RES, AERO_COEFF, HEAT)'
            text = set_value(text, 'HISTORY_OUTPUT', '(ITER, RMS_RES, AERO_COEFF, HEAT, LINSOL)')
        assert settings(text)['RESTART_SOL'] == 'YES' and 'RESTART_ITER' not in settings(text)
        if not exported_folder:
            os.link(str(source/'mesh.su2'), str(folder/'mesh.su2'))
        if checkpoint is None:
            (folder/'initial_flow.csv').write_bytes(selected.read_bytes())
        else:
            np.savetxt(str(folder/'initial_flow.csv'), data, delimiter=',', fmt='%.17g', comments='',
                       header=','.join('"'+name+'"' for name in columns))
        saved, _, _, rho, _, thermal = state(folder/'initial_flow.csv', True)
        assert np.array_equal(saved, data)
        recovered = core.state('rhoe', rho, thermal)
        assert np.isfinite(recovered['T']).all()
        locks.update(core.locks)
        locks[str(Path(__file__).resolve())] = sha(Path(__file__).resolve())
        init = dict(kind=KIND, sourceDirectory=str(source), sourceLocalIteration=checkpoint,
            sourceStatePath=str(selected), sourceStateSHA256=sha(selected), priorSuccessfulUpdates=prior,
            localIterationStartsAt=0, sourceKind='prepared-initial-state' if checkpoint is None else 'native-checkpoint',
            samePhysicsAndNumerics=True, sourceRuntimeFlavor=manifest['requiredRuntimeFlavor'],
            runtimeChange='stage output and first saved iteration only',
            coldRestart='conservative rho/rhoU/rhoE and post-update k/omega; rebuild primitives and eddy viscosity; reset adaptive CFL and solver histories',
            restartTemperatureRangeK=[float(min(recovered['T'])), float(max(recovered['T']))],
            sourceSHA256=locks, historicalSourceRebinding=rebinding,
            initialStateSHA256=sha(folder/'initial_flow.csv'), physicsQualified=False, newQualifiedArea=0)
        if linear_history:
            init['historyDiagnostics'] = dict(addedGroup='LINSOL', changesNumerics=False)
        if runtime_flavor == OBSERVED_FLAVOR:
            init['runtimeChange'] = 'phase-v2 numerics; C++ termination observation and 17-digit ASCII output'
            manifest['runtimeObservation'] = dict(OBSERVATION)
        if runtime_flavor == ADMISSIBLE_FLAVOR:
            init['samePhysicsAndNumerics'] = False
            init['runtimeChange'] = 'ideal-air domain rejection restores Solution_Old with current k; invalid old state terminates'
            init['numericalChange'] = dict(ADMISSIBILITY)
            manifest['runtimeAdmissibility'] = dict(ADMISSIBILITY)
        atomic_write_json(folder/'ideal_air_restart.json', init)
        (folder/'case.cfg').write_text(text, encoding='utf8')
        contract = copy.deepcopy(old_contract)
        contract.update(schemaVersion=3, initialization=KIND, initializationSHA256=sha(folder/'ideal_air_restart.json'),
                        checkpointContract=plan, configurationSHA256=sha(folder/'case.cfg'), meshSHA256=sha(folder/'mesh.su2'))
        atomic_write_json(folder/'sst_contract.json', manifest['sst'])
        atomic_write_json(folder/'ideal_air_contract.json', contract)
        manifest.update(initialization=init, idealAir=contract, requiredRuntimeFlavor=runtime_flavor,
            idealAirContractSHA256=sha(folder/'ideal_air_contract.json'), sstContractSHA256=sha(folder/'sst_contract.json'))
        manifest['flight'].update(iterations=iterations, outputEvery=every)
        atomic_write_json(folder/'case_manifest.json', manifest)
        record.update(status='terminal', exitCode=0, elapsedSeconds=time.monotonic()-started,
                      initializationSHA256=contract['initializationSHA256'])
        atomic_write_json(marker, record)
        from su2ThermalRuntime import resolve
        resolve(folder, text)
        return manifest
    except Exception as error:
        record.update(status='terminal', exitCode=1, elapsedSeconds=time.monotonic()-started,
                      error=type(error).__name__+': '+str(error))
        atomic_write_json(marker, record)
        raise


def validate_source(folder, c, manifest, contract):
    """模块三。正常启动入口复核来源、数值设置、原节点和相容冷重启热状态。"""
    marker = folder/'ideal_air_restart_preparation.json'
    prepared = read(marker)
    assert prepared['status'] == 'terminal' and prepared['exitCode'] == 0
    init_path = folder/'ideal_air_restart.json'
    assert sha(init_path) == contract['initializationSHA256'] == prepared['initializationSHA256']
    init = read(init_path)
    assert init == manifest['initialization'] and init['kind'] == KIND
    if manifest['requiredRuntimeFlavor'] == OBSERVED_FLAVOR:
        assert manifest['runtimeObservation'] == OBSERVATION
        assert init['runtimeChange'] == 'phase-v2 numerics; C++ termination observation and 17-digit ASCII output'
    if manifest['requiredRuntimeFlavor'] == ADMISSIBLE_FLAVOR:
        assert manifest['runtimeAdmissibility'] == init['numericalChange'] == ADMISSIBILITY
        assert init['samePhysicsAndNumerics'] is False
        assert init['runtimeChange'] == 'ideal-air domain rejection restores Solution_Old with current k; invalid old state terminates'
    else:
        assert init['samePhysicsAndNumerics'] is True and 'numericalChange' not in init
    assert init['localIterationStartsAt'] == 0 and not init['physicsQualified']
    assert c['SOLUTION_FILENAME'] == 'initial_flow' and c['READ_BINARY_RESTART'] == 'NO'
    assert 'RESTART_ITER' not in c and 'continuation' not in manifest
    locks = dict(init['sourceSHA256'])
    for name, digest in locks.items():
        assert sha(Path(name)) == digest, name
    source = Path(init['sourceDirectory'])
    old = settings((source/'case.cfg').read_text(encoding='utf8'))
    expected = dict(old)
    if 'historyDiagnostics' in init:
        assert init['historyDiagnostics'] == dict(addedGroup='LINSOL', changesNumerics=False)
        assert init['sourceLocalIteration'] is not None
        assert old['HISTORY_OUTPUT'] == '(ITER, RMS_RES, AERO_COEFF, HEAT)'
        expected['HISTORY_OUTPUT'] = '(ITER, RMS_RES, AERO_COEFF, HEAT, LINSOL)'
    assert {k: v for k, v in c.items() if k not in OUTPUT_KEYS} == {k: v for k, v in expected.items() if k not in OUTPUT_KEYS}
    old_manifest = read(source/'case_manifest.json')
    changed_contract_keys = {'schemaVersion', 'initialization', 'initializationSHA256',
                             'configurationSHA256', 'checkpointContract'}
    assert {key: value for key, value in contract.items() if key not in changed_contract_keys} == {
        key: value for key, value in old_manifest['idealAir'].items() if key not in changed_contract_keys}, '物性或参考态合同变化'
    assert manifest['expectedVelocity'] == old_manifest['expectedVelocity']
    assert manifest['sst'] == old_manifest['sst'] and manifest['thermalBoundary'] == old_manifest['thermalBoundary']
    assert contract['meshSHA256'] == old_manifest['idealAir']['meshSHA256']
    assert contract['thermodynamicModel'] == old_manifest['idealAir']['thermodynamicModel']
    assert sha(folder/'initial_flow.csv') == init['initialStateSHA256']
    initial, columns, nd, rho, _, thermal = state(folder/'initial_flow.csv', True)
    points, _, wall = mesh(folder/'mesh.su2')
    assert nd == 3 and len(initial) == manifest['nodes'] == len(points)
    assert np.allclose(initial[:, 1:4], points, rtol=2e-14, atol=1e-14)
    selected = Path(init['sourceStatePath'])
    assert sha(selected) == init['sourceStateSHA256'] == locks[str(selected)]
    original = np.loadtxt(str(selected), delimiter=',', skiprows=1, usecols=range(len(columns)))
    assert np.array_equal(initial, original), '冷重启守恒态改变'
    step = init['sourceLocalIteration']
    if step is None:
        assert init['sourceKind'] == 'prepared-initial-state' and init['priorSuccessfulUpdates'] == 0
        assert selected == source/'initial_flow.csv'
    else:
        receipt = read(source/'engine_process.json')
        assert selected == source/('restart_flow_%06d.csv' % step)
        assert sha(source/'engine_process.json') == locks[str(source/'engine_process.json')]
        assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
        assert init['priorSuccessfulUpdates'] == old_manifest['initialization']['priorSuccessfulUpdates']+step+1
    plan = contract['checkpointContract']
    assert plan == output_plan(len(points), len(np.unique([v for face in wall for v in face[1:]])),
        len(wall), int(c['ITER']), int(c['OUTPUT_WRT_FREQ']), plan['maximumOutputBytes'])
    assert manifest['flight']['iterations'] == int(c['ITER']) and manifest['flight']['outputEvery'] == int(c['OUTPUT_WRT_FREQ'])
    core = IdealAir()
    core.state('rhoe', rho, thermal)
    locks.update(core.locks)
    for path in (marker, init_path, folder/'initial_flow.csv', Path(__file__).resolve(),
                 Path(__file__).with_name('su2IdealAirCheckpoints.py')):
        locks[str(path)] = sha(path)
    return locks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--folder', required=True)
    parser.add_argument('--iterations', required=True, type=int)
    parser.add_argument('--every', required=True, type=int)
    parser.add_argument('--maximum-mib', required=True, type=int)
    parser.add_argument('--checkpoint', type=int)
    parser.add_argument('--source-snapshot')
    parser.add_argument('--linear-history', action='store_true')
    parser.add_argument('--runtime-flavor', choices=(PHASE_FLAVOR, OBSERVED_FLAVOR, ADMISSIBLE_FLAVOR), default=PHASE_FLAVOR)
    args = parser.parse_args()
    prepare(args.source, args.folder, args.iterations, args.every, args.maximum_mib*1024**2,
            args.checkpoint, args.source_snapshot, linear_history=args.linear_history, runtime_flavor=args.runtime_flavor)
