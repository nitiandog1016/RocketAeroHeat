"""模块用途：从具名已热审层流保存态冷重启；同网格/Tw，仅启用固定官方自适应推进。"""
import argparse
import copy
import csv
import os
from pathlib import Path
import time

import numpy as np

from idealAirProperties import IdealAir
from prepareNACASpatialCase import FLAVOR, SPATIAL, KIND as SOURCE_KIND
from prepareSU2NACALaminar import STAGE, capacity
from su2IdealAirInterface import MODEL, R, settings
from su2NACALaminarInterface import mesh, load_state, inlet_state_check
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save

KIND = 'naca-tn3344-adaptive-cold-restart-v1'
INITIALIZATION = 'same-mesh-spatial-wall-laminar-state-v1'
PLAN = ROOT / 'docs/validation/v172/PLAN.md'
PROOF = ROOT / 'docs/validation/v171/evidence_manifest.json'
PROOF_SHA256 = 'b20ee03659b14bbaff1e2b1db3ad9e05cb932f43cf13c45e29054d29a9155435'
VARIANTS = ('small-adaptive', 'validation-adaptive')
INITIAL_COLUMNS = ['PointID', 'x', 'y', 'z', 'Density', 'Momentum_x', 'Momentum_y', 'Momentum_z', 'Energy']
ADAPTIVE = dict(kind='su2-official-adaptive-cfl-v1', initialCFL=.4,
    parameters=[.1, 1.2, .4, 100., .001, 0.], multigridLevels=0,
    samePhysicalModel=True, sameSpatialDiscretization=True, samePhysicsAndNumerics=False)
ACTUAL_PLAN_MARKER = 'ADAPTIVE_ACTUAL_CANDIDATE_REGISTERED_V172'


def specification(variant):
    """模块一：只允许两个具名组合；实际候选仍须另行完成计划登记。"""
    assert variant in VARIANTS, '未登记的自适应变体'
    actual = variant == 'validation-adaptive'
    return dict(variant=variant, actualValidationMesh=actual, updates=1000 if actual else 8,
        every=250 if actual else 2, maximumSeconds=1800 if actual else 30, threads=2,
        maximumOutputBytes=(1024 if actual else 32)*1024**2, requiredRuntimeFlavor=FLAVOR,
        sourceDirectory='outputs/v171_validation/' + ('real01' if actual else 'spatial_profile01'),
        sourceLocalIteration=1500 if actual else 3,
        sourceAudit='outputs/v171_validation/' + ('audit_real02' if actual else 'audit_small04') + '/report.json',
        sourceVariant='validation-profile' if actual else 'spatial-profile',
        resourceProfile='v172-naca-adaptive-validation' if actual else 'v172-naca-adaptive-small')


def source_files():
    """模块二：只锁运行所需源和计划，不锁尚在编辑的评审文稿。"""
    names = ('prepareNACAAdaptiveContinuation.py', 'su2NACAAdaptiveInterface.py', 'prepareNACASpatialCase.py',
        'su2NACASpatialInterface.py', 'prepareSU2NACALaminar.py', 'su2NACALaminarInterface.py',
        'su2IdealAirInterface.py', 'su2IdealAirCheckpoints.py', 'su2ThermalRuntime.py', 'runEmbeddedSU2.py',
        'analyzeNACASpatialHeat.py', 'nacaSavedFrameStudy.py', 'nacaSpatialWallFlux.py', 'nacaWallFlux.py',
        'nacaBoundaryLayerEdge.py', 'nacaMeasurementSegments.py', 'naca3344Model.py',
        'idealAirProperties.py', 'buildIdealAirBatch.py', 'auditSU2SSTState.py')
    return [ROOT/'scripts'/name for name in names] + [PLAN, PROOF]


def _proof_rows(proof):
    return proof['files'] + proof['finalSources'] + proof['localOnly']


def _proof_identity(path, proof, locks):
    """直接证据必须匹配冻结清单；其他输入通过已冻结原回执的SHA传递绑定。"""
    path = Path(path).resolve(); relative = path.relative_to(ROOT).as_posix()
    rows = [row for row in _proof_rows(proof) if row.get('source', row.get('path')) == relative]
    assert rows and len({row['sha256'] for row in rows}) == 1, '证据身份缺失或冲突: ' + relative
    digest = rows[0]['sha256']; assert sha(path) == digest, '原证据发生变化: ' + relative
    locks[str(path)] = digest
    return digest


def _historical_inputs(identities, proof, locks, rebinding):
    """仅源码可重绑定到proof已有准确快照；物理输入、运行包与回执必须保持原件。"""
    for name, digest in identities.items():
        path = Path(name).resolve(); relative = path.relative_to(ROOT)
        if str(path) in locks and locks[str(path)] == digest:
            continue
        if path.is_file() and sha(path) == digest:
            locks[str(path)] = digest
            continue
        source = relative.parts[0] in ('scripts', 'tests', 'third_party') and path.suffix.lower() in (
            '.py', '.m', '.hpp', '.h', '.cpp', '.c', '.inc')
        assert source, '非源码输入改变或缺失，禁止快照代填: ' + name
        candidates = [ROOT/row['archive'] for row in _proof_rows(proof)
            if row.get('archive') and row['sha256'] == digest]
        assert candidates, 'proof没有该历史源码准确快照: ' + name
        target = candidates[0].resolve(); target.relative_to(ROOT/'docs/validation/v171/evidence')
        assert sha(target) == digest, '历史准确快照改变'
        locks[str(target)] = digest
        item = dict(original=str(path), exactSnapshot=str(target), sha256=digest)
        if item not in rebinding:
            rebinding.append(item)


def adaptive_config(source_text, spec):
    """模块三：配置只替换推进开关/固定参数、更新预算与保存周期。"""
    values = settings(source_text)
    assert values['CFL_NUMBER'] == '0.4' and values['CFL_ADAPT'] == 'NO' and values['MGLEVEL'] == '0'
    assert 'CFL_ADAPT_PARAM' not in values and 'RESTART_ITER' not in values
    changes = {'CFL_ADAPT': 'YES', 'ITER': str(spec['updates']), 'OUTPUT_WRT_FREQ': str(spec['every'])}
    lines = []
    for line in source_text.splitlines():
        key = line.split('=', 1)[0].strip()
        lines.append(key + '= ' + changes[key] if key in changes else line)
    lines.append('CFL_ADAPT_PARAM= (0.1, 1.2, 0.4, 100, 0.001, 0)')
    return '\n'.join(lines) + '\n'


def source_state(variant):
    """模块四：原超时/正常身份、精选完整双帧、热审核和守恒热状态逐项验证。"""
    spec = specification(variant); source = ROOT/spec['sourceDirectory']; audit_path = ROOT/spec['sourceAudit']
    assert sha(PROOF) == PROOF_SHA256, 'v171冻结证据清单身份改变'
    proof = read(PROOF); locks = {str(PROOF): PROOF_SHA256}; rebinding = []
    for name in ('engine_process.json', 'case_manifest.json', 'ideal_air_contract.json', 'case.cfg',
                 'mesh.su2', 'mesh_arrays.npz', 'geometry.json', 'wall_temperature.csv',
                 'rah_wall_temperature_applied.csv', 'history.csv'):
        _proof_identity(source/name, proof, locks)
    _proof_identity(audit_path, proof, locks)
    receipt = read(source/'engine_process.json'); manifest = read(source/'case_manifest.json')
    old_contract = read(source/'ideal_air_contract.json'); audit = read(audit_path)
    assert receipt['status'] == 'terminal' and receipt['runtimeFlavor'] == FLAVOR
    for item in (receipt, audit):
        stable = item.get('inputsStable'); inputs = item.get('inputSHA256')
        assert isinstance(stable, dict) and stable and all(value is True for value in stable.values())
        assert isinstance(inputs, dict) and set(stable) == set(inputs)
        _historical_inputs(inputs, proof, locks, rebinding)
    assert not any(receipt.get(key, False) for key in ('cancelled', 'outputBudgetExceeded', 'supervisorError'))
    if spec['actualValidationMesh']:
        assert receipt['exitCode'] == 1 and receipt['timedOut'] is True
        assert audit['savedFrameInterfacePassed'] and not audit['interfacePassed'] and not audit['normalRunCompleted']
        assert audit['auditScope'] == 'complete-saved-frames-from-abnormal-task'
    else:
        assert receipt['exitCode'] == 0 and not receipt.get('timedOut', False)
        assert receipt['checkpoints']['normalScheduleComplete'] and audit['interfacePassed']
    assert manifest['caseKind'] == SOURCE_KIND and manifest['requiredRuntimeFlavor'] == FLAVOR
    assert manifest['specification']['variant'] == spec['sourceVariant']
    assert manifest['initialization']['priorSuccessfulUpdates'] == 0
    assert old_contract == manifest['idealAir'] and sha(source/'ideal_air_contract.json') == manifest['idealAirContractSHA256']
    assert not old_contract['turbulent'] and old_contract['thermodynamicModel'] == MODEL
    assert old_contract['energyEntropyReference'] == 'OTH' and old_contract['temperatureRangeK'] == [60, 2000]
    assert old_contract['gasConstant'] == R and 'sst' not in manifest and not (source/'sst_contract.json').exists()
    assert audit['stage'] == STAGE and not audit['experimentalQualification'] and audit['newQualifiedArea_m2'] == 0
    assert audit['inputSHA256'][str(source/'engine_process.json')] == locks[str(source/'engine_process.json')]

    step = spec['sourceLocalIteration']
    entry = next(row for row in receipt['checkpoints']['complete'] if row['localIteration'] == step)
    assert entry['stage'] == STAGE and entry['successfulUpdates'] == entry['seriesSuccessfulUpdates'] == step+1
    for key, ext in (('volume', 'csv'), ('surface', 'vtk')):
        expected = ('restart_flow_' if key == 'volume' else 'surface_flow_') + ('%06d.' % step) + ext
        assert entry[key] == expected
        digest = _proof_identity(source/expected, proof, locks)
        assert digest == entry[key+'SHA256'] == receipt['outputSHA256'][expected]
    assert audit['inputSHA256'][str(source/entry['volume'])] == entry['volumeSHA256']
    frame = next(row for row in audit['frames'] if row['localIteration'] == step)
    assert frame['interfacePassed'] and frame['energyMaximumRelativeError'] < 2e-11 and frame['eosMaximumRelativeError'] < 2e-11
    assert frame['nativeHeatMaximumScaledError'] < 2e-10 and frame['powerScaledDifference'] < 2e-9
    assert frame['paperCpMaximumRelativeError'] < 2e-8 and frame['viscosityMaximumRelativeError'] < 2e-11

    # 只重验精选双帧，原周期完整性与终态由已冻结来源回执保持；不解析无关大帧。
    from su2NACASpatialInterface import inventory as old_inventory
    selected_receipt = copy.deepcopy(receipt)
    selected_receipt['outputSHA256'] = {name: digest for name, digest in receipt['outputSHA256'].items()
        if name in ('history.csv', entry['volume'], entry['surface'])}
    assert sha(source/'history.csv') == receipt['outputSHA256']['history.csv']
    selected = old_inventory(source, selected_receipt)
    assert selected['historyContiguous'] and not selected['rejected'] and len(selected['complete']) == 1
    assert selected['complete'][0] == entry
    p, cells, markers = mesh(source/'mesh.su2'); meta = read(source/'geometry.json')
    assert set(markers) == {'wall', 'outer', 'outlet', 'sym_y', 'sym_z'}
    data, fields = load_state(source/entry['volume']); initial = data[:, :9].copy()
    assert np.array_equal(initial[:, 1:4], p) and len(initial) == manifest['nodes']
    rho = initial[:, 4]; kinetic = .5*np.sum((initial[:, 5:8]/rho[:, None])**2, axis=1)
    thermal = initial[:, 8]/rho-kinetic
    core = IdealAir(); recovered = core.state('rhoe', rho, thermal); locks.update(core.locks)
    thermo = core.state('PT', fields['Pressure'], fields['Temperature'])
    energy_error = float(np.max(abs(thermal-thermo['e'])/np.maximum(1., abs(thermo['e']))))
    eos_error = float(np.max(abs(rho/thermo['rho']-1)))
    temperature_error = float(np.max(abs(recovered['T']-fields['Temperature'])/np.maximum(1., abs(fields['Temperature']))))
    pressure_error = float(np.max(abs(recovered['p']-fields['Pressure'])/np.maximum(1., abs(fields['Pressure']))))
    assert max(energy_error, eos_error) < 2e-11 and max(temperature_error, pressure_error) < 2e-10
    assert np.isfinite(recovered['T']).all() and recovered['T'].min() >= 60 and recovered['T'].max() <= 2000
    from su2NACASpatialInterface import read_profile
    from nacaSpatialWallFlux import verify_applied
    wall_ids = np.unique([node for face in markers['wall'] for node in face[1:]])
    _, tw = read_profile(source/'wall_temperature.csv', p[wall_ids, 0])
    applied = verify_applied(source/'rah_wall_temperature_applied.csv', p, {'ids': wall_ids}, tw)
    runtime_path = ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json'; runtime = read(runtime_path)
    executable = Path(runtime['executable']).resolve(); executable.relative_to(runtime_path.parent)
    for path in (runtime_path, executable):
        assert sha(path) == receipt['inputSHA256'][str(path)]
        locks[str(path)] = sha(path)
    assert runtime['flavor'] == FLAVOR and runtime['spatialWall'] == SPATIAL and runtime['thermodynamicModel'] == MODEL
    source_text = (source/'case.cfg').read_text(encoding='utf8')
    inlet_state_check(source_text, manifest['publicBenchmark']['conditions'], manifest['expectedVelocity'])
    source_identity = dict(sourceDirectory=str(source), sourceLocalIteration=step,
        sourceStatePath=str(source/entry['volume']), sourceStateSHA256=entry['volumeSHA256'],
        sourceSurfacePath=str(source/entry['surface']), sourceSurfaceSHA256=entry['surfaceSHA256'],
        sourceReceiptPath=str(source/'engine_process.json'), sourceReceiptSHA256=sha(source/'engine_process.json'),
        sourceAuditPath=str(audit_path), sourceAuditSHA256=sha(audit_path), sourceStage=STAGE,
        sourceRuntimeFlavor=FLAVOR, sourceRuntimeManifestSHA256=sha(runtime_path), sourceExecutableSHA256=sha(executable),
        sourceTermination=dict(status=receipt['status'], exitCode=receipt['exitCode'], timedOut=receipt.get('timedOut', False),
            cancelled=receipt.get('cancelled', False), outputBudgetExceeded=receipt.get('outputBudgetExceeded', False)),
        sourceLastFiniteHistoryLocalIteration=receipt['checkpoints']['lastFiniteHistoryLocalIteration'],
        sourceAbnormalInstantStateAvailable=False, priorSuccessfulUpdates=step+1,
        sourceInterfacePassed=bool(audit['interfacePassed']), sourceSavedFrameInterfacePassed=bool(frame['interfacePassed']),
        historicalSourceRebinding=rebinding,
        restartThermalCheck=dict(energyMaximumRelativeError=energy_error, eosMaximumRelativeError=eos_error,
            temperatureMaximumScaledError=temperature_error, pressureMaximumScaledError=pressure_error,
            temperatureMinimumK=float(recovered['T'].min()), temperatureMaximumK=float(recovered['T'].max())),
        sourceAppliedTemperatureCheck=applied)
    return dict(spec=spec, source=source, manifest=manifest, contract=old_contract, text=source_text,
        initial=initial, points=p, cells=cells, markers=markers, meta=meta, identity=source_identity, locks=locks)


def prepare(folder, variant):
    """模块五：独占新目录，原网格/表只读硬链接，17位保存原五个守恒量。"""
    from runEmbeddedSU2 import atomic_write_json
    folder = Path(folder).resolve(); folder.relative_to(ROOT/'outputs'); folder.mkdir(parents=True, exist_ok=False)
    marker = folder/'adaptive_preparation.json'; started = time.monotonic()
    atomic_write_json(marker, dict(status='preparing', exitCode=None))
    try:
        checked = source_state(variant); spec = checked['spec']; source = checked['source']; old = checked['manifest']
        locks = dict(checked['locks']); locks.update({str(path): sha(path) for path in source_files()})
        (folder/'executed_prepare.py').write_bytes(Path(__file__).read_bytes())
        for name in ('mesh.su2', 'mesh_arrays.npz', 'geometry.json', 'wall_temperature.csv'):
            os.link(source/name, folder/name)
        np.savetxt(str(folder/'initial_flow.csv'), checked['initial'], delimiter=',', fmt='%.17g', comments='',
            header=','.join('"'+name+'"' for name in INITIAL_COLUMNS))
        initial = np.loadtxt(str(folder/'initial_flow.csv'), delimiter=',', skiprows=1)
        assert np.array_equal(initial, checked['initial']), '17位初场没有逐元素保持原守恒态'
        (folder/'case.cfg').write_text(adaptive_config(checked['text'], spec), encoding='ascii')
        initialization = dict(checked['identity'], kind=INITIALIZATION, localIterationStartsAt=0,
            initialStateSHA256=sha(folder/'initial_flow.csv'), sourceSHA256=locks,
            samePhysicalModel=True, sameSpatialDiscretization=True, samePhysicsAndNumerics=False,
            uninterruptedTrajectoryEquivalent=False, adaptiveCFLReset=True, solverHistoryReset=True,
            coldRestart='same laminar rho/rhoU/rhoE; rebuild primitive/gradient/limiter and solver/adaptive histories',
            boundaryLayerSolved=False, physicsQualified=False, newQualifiedArea=0)
        save(folder/'adaptive_initialization.json', initialization)
        wallids = np.unique([i for f in checked['markers']['wall'] for i in f[1:]])
        contract = copy.deepcopy(checked['contract'])
        contract.update(schemaVersion=7, initialization=INITIALIZATION,
            initializationSHA256=sha(folder/'adaptive_initialization.json'), configurationSHA256=sha(folder/'case.cfg'),
            checkpointContract=capacity(len(initial), len(wallids), len(checked['markers']['wall']),
                spec['updates'], spec['every'], spec['maximumOutputBytes']))
        save(folder/'ideal_air_contract.json', contract)
        manifest = dict(caseKind=KIND, specification=spec, requiredRuntimeFlavor=FLAVOR,
            idealAir=contract, idealAirContractSHA256=sha(folder/'ideal_air_contract.json'),
            initialization=initialization, advancement=copy.deepcopy(ADAPTIVE),
            thermalBoundary=copy.deepcopy(old['thermalBoundary']), thermalProfileSHA256=old['thermalProfileSHA256'],
            publicBenchmark=copy.deepcopy(old['publicBenchmark']), nodes=old['nodes'], elements=old['elements'],
            expectedVelocity=old['expectedVelocity'], physicalQualified=False, planSHA256=sha(PLAN),
            historicalProofSHA256=PROOF_SHA256)
        save(folder/'case_manifest.json', manifest)
        request = dict(workDirectory=str(folder), configFile='case.cfg', meshFile='mesh.su2',
            maximumSeconds=spec['maximumSeconds'], threads=spec['threads'], resourceProfile=spec['resourceProfile'],
            sourceFiles=[path.relative_to(ROOT).as_posix() for path in source_files()])
        save(folder/'request.json', request)
        assert all(sha(Path(path)) == digest for path, digest in locks.items())
        atomic_write_json(marker, dict(status='terminal', exitCode=0, elapsedSeconds=time.monotonic()-started,
            initializationSHA256=contract['initializationSHA256'], sourceInputsStable=True))
        return manifest
    except Exception as exc:
        atomic_write_json(marker, dict(status='terminal', exitCode=1, elapsedSeconds=time.monotonic()-started,
            error=type(exc).__name__+': '+str(exc)))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--folder', required=True)
    parser.add_argument('--variant', choices=VARIANTS, required=True); args = parser.parse_args()
    prepare(args.folder, args.variant)
