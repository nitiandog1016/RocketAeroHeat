"""模块用途：具名层流空间Tw限步组合合同，拒绝旧包冒充、反馈耦合或全域保护声明。"""
import csv
from pathlib import Path
import numpy as np

from prepareNACASpatialStep import (KIND, INITIALIZATION, PLAN, PROOF_SHA256, FLAVOR, SPATIAL,
    ADAPTIVE, ACTUAL_PLAN_MARKER, INITIAL_COLUMNS, NUMERICAL_CHANGE, QUERY_CONTEXT, specification, source_state, source_files, adaptive_config)
from prepareSU2NACALaminar import capacity
from su2IdealAirInterface import MODEL, R, OBSERVATION, CONSERVATIVE_STEP
from su2NACALaminarInterface import mesh, inlet_state_check
from su2ThermalRuntime import read, sha


def validate(folder, manifest, text, runtime):
    """模块一：重验固定来源及新物理输入；源码变更只能由原精确快照承接。"""
    folder = Path(folder).resolve()
    assert manifest['caseKind'] == KIND and manifest['historicalProofSHA256'] == PROOF_SHA256
    spec = specification(manifest['specification']['variant']); assert spec == manifest['specification']
    assert manifest['requiredRuntimeFlavor'] == runtime['flavor'] == FLAVOR
    assert runtime['spatialWall'] == SPATIAL and runtime['thermodynamicModel'] == MODEL
    assert runtime['numericalParentFlavor'] == 'rah-ideal-air-spatial-wall-v1'
    assert runtime['numericalBaselineFlavor'] == FLAVOR
    assert all(runtime[k] == OBSERVATION[k] for k in ('asciiSignificantDigits','terminationObservation'))
    assert runtime['conservativeStep'] == CONSERVATIVE_STEP and runtime['queryContext'] == QUERY_CONTEXT
    assert manifest['numericalChange'] == NUMERICAL_CHANGE
    assert runtime['samePhysicsAndNumerics'] is False and runtime['samePhysicsAndNumericsWithoutSpatialTable'] is False
    assert manifest['planSHA256'] == sha(PLAN) and manifest['advancement'] == ADAPTIVE
    assert manifest['physicalQualified'] is False and 'sst' not in manifest and not (folder/'sst_contract.json').exists()
    prepared = read(folder/'step_preparation.json')
    assert prepared['status'] == 'terminal' and prepared['exitCode'] == 0 and prepared['sourceInputsStable'] is True
    checked = source_state(spec['variant']); old = checked['manifest']; source = checked['source']
    assert text == adaptive_config(checked['text'], spec), '新配置包含未登记的改变'
    assert manifest['publicBenchmark'] == old['publicBenchmark'] and manifest['expectedVelocity'] == old['expectedVelocity']
    assert manifest['thermalBoundary'] == old['thermalBoundary'] and manifest['thermalProfileSHA256'] == old['thermalProfileSHA256']
    assert manifest['nodes'] == old['nodes'] and manifest['elements'] == old['elements']
    for name in ('mesh.su2', 'mesh_arrays.npz', 'geometry.json', 'wall_temperature.csv'):
        assert sha(folder/name) == sha(source/name), '原网格、法向列或壁温表改变: ' + name
    p, cells, markers = mesh(folder/'mesh.su2')
    assert np.array_equal(p, checked['points']) and cells == checked['cells'] and markers == checked['markers']
    assert set(markers) == {'wall', 'outer', 'outlet', 'sym_y', 'sym_z'}
    with np.load(str(folder/'mesh_arrays.npz')) as arrays:
        assert np.array_equal(arrays['points'], p)
        for key, kind in (('prisms', 13), ('hexes', 12)):
            assert np.array_equal(arrays[key], np.array([cell[1:] for cell in cells if cell[0] == kind]))
    inlet_state_check(text, manifest['publicBenchmark']['conditions'], manifest['expectedVelocity'])
    contract = read(folder/'ideal_air_contract.json')
    assert contract == manifest['idealAir'] and sha(folder/'ideal_air_contract.json') == manifest['idealAirContractSHA256']
    assert contract['schemaVersion'] == 8 and contract['initialization'] == INITIALIZATION
    assert contract['thermodynamicModel'] == MODEL and contract['gasConstant'] == R
    assert contract['energyEntropyReference'] == 'OTH' and contract['temperatureRangeK'] == [60, 2000]
    assert not contract['turbulent'] and not contract['physicsQualified'] and contract['newQualifiedArea'] == 0
    for name, key in (('case.cfg', 'configurationSHA256'), ('mesh.su2', 'meshSHA256'),
                      ('geometry.json', 'geometrySHA256'), ('mesh_arrays.npz', 'meshArraysSHA256')):
        assert sha(folder/name) == contract[key]
    wallids = np.unique([i for face in markers['wall'] for i in face[1:]])
    assert contract['checkpointContract'] == capacity(len(p), len(wallids), len(markers['wall']),
        spec['updates'], spec['every'], spec['maximumOutputBytes'])
    initialization = read(folder/'step_initialization.json')
    assert initialization == manifest['initialization']
    assert sha(folder/'step_initialization.json') == contract['initializationSHA256'] == prepared['initializationSHA256']
    for key, value in checked['identity'].items():
        assert initialization[key] == value, '冷重启来源身份改变: ' + key
    assert initialization['kind'] == INITIALIZATION and initialization['localIterationStartsAt'] == 0
    assert initialization['priorSuccessfulUpdates'] == spec['sourceLocalIteration']+1
    assert initialization['samePhysicalModel'] and initialization['sameSpatialDiscretization']
    assert not initialization['samePhysicsAndNumerics'] and not initialization['uninterruptedTrajectoryEquivalent']
    assert initialization['adaptiveCFLReset'] and initialization['solverHistoryReset']
    assert initialization['numericalChange'] == NUMERICAL_CHANGE
    assert not initialization['boundaryLayerSolved'] and not initialization['physicsQualified'] and initialization['newQualifiedArea'] == 0
    assert sha(folder/'initial_flow.csv') == initialization['initialStateSHA256']
    with (folder/'initial_flow.csv').open(encoding='utf8') as stream:
        assert next(csv.reader(stream)) == INITIAL_COLUMNS
    initial = np.loadtxt(str(folder/'initial_flow.csv'), delimiter=',', skiprows=1)
    assert np.isfinite(initial).all() and np.array_equal(initial, checked['initial']), '初场守恒列或原节点改变'
    locks = dict(initialization['sourceSHA256'])
    expected_locks = dict(checked['locks']); expected_locks.update({str(path): sha(path) for path in source_files()})
    assert locks == expected_locks, '新源锁缺失、多余或发生变化'
    assert all(sha(Path(path)) == digest for path, digest in locks.items())
    for name in ('step_preparation.json', 'step_initialization.json', 'case_manifest.json', 'ideal_air_contract.json',
                 'geometry.json', 'mesh_arrays.npz', 'initial_flow.csv', 'wall_temperature.csv'):
        locks[str(folder/name)] = sha(folder/name)
    return locks


def resource_budget(request, manifest):
    """模块二：新profile独占预算；实际候选必须有事前计划登记标记。"""
    assert manifest['caseKind'] == KIND and manifest['planSHA256'] == sha(PLAN)
    spec = specification(manifest['specification']['variant']); assert spec == manifest['specification']
    assert request.get('resourceProfile') == spec['resourceProfile']
    assert request['threads'] == spec['threads'] and request['maximumSeconds'] == spec['maximumSeconds']
    assert manifest['requiredRuntimeFlavor'] == FLAVOR and manifest['advancement'] == ADAPTIVE
    assert manifest['initialization']['priorSuccessfulUpdates'] == spec['sourceLocalIteration']+1
    plan = manifest['idealAir']['checkpointContract']
    assert plan['maximumUpdates'] == spec['updates'] and plan['every'] == spec['every']
    assert plan['maximumOutputBytes'] == spec['maximumOutputBytes']
    if spec['actualValidationMesh']:
        assert ACTUAL_PLAN_MARKER in PLAN.read_text(encoding='utf8'), '实际自适应候选尚未完成事前登记'
        assert manifest['nodes'] == 190334 and manifest['elements'] == 179568
    else:
        assert manifest['nodes'] == 549 and manifest['elements'] == 384
    return float(spec['maximumSeconds']), spec['threads']


def inventory(folder, receipt):
    """模块三：共用已审核清点核心，仅增加明确的冷重启累计offset。"""
    from su2NACASpatialInterface import _inventory_for_spec
    folder = Path(folder).resolve(); manifest = read(folder/'case_manifest.json')
    assert manifest['caseKind'] == KIND
    spec = specification(manifest['specification']['variant']); assert spec == manifest['specification']
    initialization = manifest['initialization']
    assert initialization['kind'] == INITIALIZATION and initialization['priorSuccessfulUpdates'] == spec['sourceLocalIteration']+1
    assert initialization['sourceLocalIteration'] == spec['sourceLocalIteration']
    return _inventory_for_spec(folder, receipt, spec, offset=initialization['priorSuccessfulUpdates'])
