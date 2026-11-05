"""准备显式二维轴对称变比热SST-Vm，复用公开输入和精确壁距。"""
import argparse
import os
from pathlib import Path
import time

import numpy as np

from idealAirProperties import IdealAir
from prepareSU2AxisymmetricInitialState import segment_projection
from prepareSU2IdealAirRestart import historical_locks
from runEmbeddedSU2 import atomic_write_json
from su2AxisymmetricAnalysis import mesh_data, integrate_scalar
from su2IdealAirCheckpoints import output_plan
from su2IdealAirInterface import AXISYMMETRIC_FLAVOR, MODEL, R, settings
from su2SSTInterface import PHASE_OUTPUT, state, wall_geometry
from su2ThermalRuntime import ROOT, read, sha

KIND = 'axisymmetric-ideal-air-SST-Vm-wall-distance-v1'


def fixture(path):
    """模块一。81点正半径圆柱壳，供同一轴对称入口的短步和取消检查。"""
    x = np.array([.064, .032, .016, .008, .004, .002, .001, .0005, 0])
    y = .1015+np.array([0, 2e-6, 8e-6, 32e-6, .000128, .0005, .002, .008, .032])
    xx, yy = np.meshgrid(x, y)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    ix = np.arange(len(points)).reshape(len(y), len(x))
    cells = np.stack((ix[:-1, :-1], ix[1:, :-1], ix[1:, 1:], ix[:-1, 1:]), axis=2).reshape(-1, 4)
    markers = dict(wall=np.column_stack((ix[0, 1:], ix[0, :-1])),
        inlet=np.column_stack((ix[:-1, 0], ix[1:, 0])),
        top=np.column_stack((ix[-1, :-1], ix[-1, 1:])),
        outlet=np.column_stack((ix[1:, -1], ix[:-1, -1])))
    with path.open('w', encoding='ascii', newline='\n') as stream:
        stream.write('NDIME= 2\nNELEM= %d\n' % len(cells))
        for cell in cells:
            stream.write('9 '+' '.join(map(str, cell))+'\n')
        stream.write('NPOIN= %d\n' % len(points))
        for index, point in enumerate(points):
            stream.write('%.17g %.17g %d\n' % (point[0], point[1], index))
        stream.write('NMARK= 4\n')
        for tag, edges in markers.items():
            stream.write('MARKER_TAG= %s\nMARKER_ELEMS= %d\n' % (tag, len(edges)))
            for edge in edges:
                stream.write('3 %d %d\n' % tuple(edge))


def initial_state(core, points, markers, geometry, c):
    """模块二。既有C1速度/总焓与正SST规则，以OTH焓反解热态。"""
    distance, nearest, weight, feet = segment_projection(points, points, markers['wall'])
    wall_ids = geometry['wallIDs']
    open_ids = np.setdiff1d(np.union1d(markers['inlet'], markers['top']), wall_ids)
    length = np.maximum(np.max(points[markers['inlet'], 0])-feet[:, 0], 0)
    thickness = np.zeros(len(points))
    inside = length > 0
    thickness[inside] = .37*length[inside]/(c['unitReynolds']*length[inside])**.2
    eta = np.ones(len(points))
    np.divide(distance, thickness, out=eta, where=inside)
    eta = np.clip(eta, 0, 1)
    blend = eta**2*(3-2*eta)
    blend[open_ids], blend[wall_ids] = 1, 0
    ends = core.state('PT', c['p'], [c['wallTemperature'], c['T']])
    target_h = (1-blend)*ends['h'][0]+blend*ends['h'][1]+.5*c['speed']**2*blend*(1-blend)
    thermo = core.state('Ph', c['p'], target_h)
    for key, values in thermo.items():
        values[blend == 0], values[blend == 1] = ends[key][0], ends[key][1]
    velocity = np.column_stack((-c['speed']*blend, np.zeros(len(points))))
    kinetic = c['k']*blend
    wall_lengths = np.full(len(points), np.nan)
    wall_lengths[wall_ids] = geometry['nearestWallDistance']
    segment = markers['wall'][nearest]
    local_length = (1-weight)*wall_lengths[segment[:, 0]]+weight*wall_lengths[segment[:, 1]]
    tw = c['wallTemperature']
    muwall = 1.716e-5*(tw/273.15)**1.5*383.55/(tw+110.4)
    nuwall = muwall/ends['rho'][0]
    near = 60*nuwall/(.075*(distance+local_length)**2)
    omega = np.exp((1-blend)*np.log(near)+blend*np.log(c['omega']))
    omega[open_ids] = c['omega']
    omega[wall_ids] = 60*nuwall/(.075*geometry['nearestWallDistance']**2)
    rho = thermo['rho']
    energy = rho*(thermo['e']+.5*np.sum(velocity**2, axis=1)+kinetic)
    values = np.column_stack((rho, rho[:, None]*velocity, energy, kinetic, omega))
    expected_h = (1-blend)*ends['h'][0]+blend*(ends['h'][1]+.5*c['speed']**2+c['k'])
    actual_h = (energy+c['p'])/rho
    error = float(max(abs(actual_h/expected_h-1)))
    assert error < 2e-11 and np.isfinite(values).all() and min(omega) > 0
    return values, dict(temperature=thermo['T'], maximumThicknessM=float(max(thickness)), totalEnthalpyError=error)


def prepare(folder, iterations, every, maximum_bytes, small=False):
    """模块三。先冻结数据及初场，再发布可由正式运行器验证的独立合同。"""
    started = time.monotonic()
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    assert not folder.exists()
    source = ROOT/'outputs/v158_validation/preparation01'
    prepared = read(source/'preparation.json')
    assert prepared['prepared'] and prepared['inputsStable']
    locks, rebound = historical_locks(prepared['sourceSHA256'], ROOT/'outputs/v158_validation/prepare01_source')
    for name, entry in prepared['files'].items():
        assert sha(source/name) == entry['sha256']
    locks.update({str(source/name): sha(source/name) for name in ('preparation.json', 'mesh.su2', 'profile_locations.csv')})
    core = IdealAir()
    locks.update(core.locks)
    folder.mkdir(parents=True)
    marker = folder/'axisymmetric_preparation.json'
    record = dict(status='preparing', exitCode=None)
    atomic_write_json(marker, record)
    try:
        if small:
            fixture(folder/'mesh.su2')
        else:
            os.link(str(source/'mesh.su2'), str(folder/'mesh.su2'))
        p, cells, markers = mesh_data(folder/'mesh.su2')
        if small:
            # 小壳出口位于X=0，独立单元/边界审核由测试覆盖。
            assert len(p) == 81 and len(cells) == 64 and min(p[:, 1]) > 0
        else:
            assert len(p) == 64521 and len(cells) == 64000
            assert sha(folder/'mesh.su2') == prepared['files']['mesh.su2']['sha256']
        geometry = wall_geometry(folder/'mesh.su2')
        c = prepared['conditions']
        initial, detail = initial_state(core, p, markers, geometry, c)
        headers = ['PointID', 'x', 'y', 'Density', 'Momentum_x', 'Momentum_y', 'Energy', 'Turb_Kin_Energy', 'Omega']
        data = np.column_stack((np.arange(len(p)), p, initial))
        np.savetxt(str(folder/'initial_flow.csv'), data, delimiter=',', fmt='%.17g', comments='',
                   header=','.join('"'+name+'"' for name in headers))
        loaded, _, _, rho, _, thermal = state(folder/'initial_flow.csv', True)
        assert np.array_equal(data, loaded)
        recovered = core.state('rhoe', rho, thermal)
        assert max(abs(recovered['T']/detail['temperature']-1)) < 2e-11
        capacity = output_plan(len(p), len(geometry['wallIDs']), len(markers['wall']), iterations,
                               every, maximum_bytes, extra_columns=7)
        old_config = ROOT/'outputs/v152_validation/acc_sbli01/case.cfg'
        config = settings(old_config.read_text('utf8'))
        locks[str(old_config)] = sha(old_config)
        for name in ('GAMMA_VALUE', 'RESTART_ITER'):
            config.pop(name, None)
        inlet = '%.17g, %.17g, %.17g, 0, 0' % (c['T'], c['p'], -c['speed'])
        changes = dict(FLUID_MODEL='IDEAL_AIR_VARIABLE_CP', GAS_CONSTANT=format(R, '.17g'),
            RESTART_SOL='YES', SOLUTION_FILENAME='initial_flow', READ_BINARY_RESTART='NO',
            FREESTREAM_TEMPERATURE=format(c['T'], '.17g'), FREESTREAM_PRESSURE=format(c['p'], '.17g'),
            MARKER_SUPERSONIC_INLET='(inlet, '+inlet+', top, '+inlet+')',
            MARKER_INLET_TURBULENT='(inlet, 0.00021, 0.001, top, 0.00021, 0.001)',
            CONV_NUM_METHOD_FLOW='ROE', ITER=str(iterations), OUTPUT_WRT_FREQ=str(every),
            SST_OPTIONS='(V1994m, VORTICITY)',
            WRT_RESTART_OVERWRITE='NO', WRT_SURFACE_OVERWRITE='NO',
            VOLUME_OUTPUT='(COORDINATES, SOLUTION, PRIMITIVE, FLOW_PHASE, RESIDUAL)',
            SCREEN_WRT_FREQ_INNER='100', REF_AREA=format(integrate_scalar(p, markers['wall'], np.ones(len(p))), '.17g'))
        config.update(changes)
        text = '\n'.join(k+'= '+v for k, v in config.items())+'\n'
        (folder/'case.cfg').write_text(text, encoding='ascii')
        for name in ('prepareSU2AxisymmetricAir.py', 'prepareSU2AxisymmetricInitialState.py', 'su2IdealAirCheckpoints.py',
                     'su2IdealAirInterface.py', 'su2SSTInterface.py', 'su2AxisymmetricAnalysis.py'):
            path = ROOT/'scripts'/name
            locks[str(path)] = sha(path)
        init = dict(kind=KIND, initialGuessOnly=True, priorSuccessfulUpdates=0, localIterationStartsAt=0,
            sourceLocalIteration=None, initialStateSHA256=sha(folder/'initial_flow.csv'), sourceSHA256=locks,
            historicalSourceRebinding=rebound, maximumThicknessM=detail['maximumThicknessM'],
            totalEnthalpyRelativeError=detail['totalEnthalpyError'], temperatureRangeK=[float(min(recovered['T'])), float(max(recovered['T']))],
            method='exact wall-segment distance; C1 velocity/total enthalpy; OTH h inverse; positive SST',
            boundaryPriority='wall initial state at shared inlet/outlet nodes', physicsQualified=False, newQualifiedArea=0)
        sst = dict(schemaVersion=1, model='SST-Vm', thermodynamics=MODEL, outputStateContract=PHASE_OUTPUT,
            intensity=c['turbulenceIntensity'], viscosityRatio=c['eddyToLaminarRatio'], Prt=.9,
            kInfinity=c['k'], omegaInfinity=c['omega'], provenance='NASA TMR WIND SST-Vm prescribed computational inflow',
            assumption='fully turbulent model; measured transition history unresolved', physicsQualified=False, newQualifiedArea=0)
        atomic_write_json(folder/'axisymmetric_initialization.json', init)
        contract = dict(schemaVersion=4, thermodynamicModel=MODEL, gasConstant=R, temperatureRangeK=[60, 2000],
            energyEntropyReference='OTH', initialization=KIND, initializationSHA256=sha(folder/'axisymmetric_initialization.json'),
            turbulent=True, freestreamSpeed=c['speed'], freestreamSoundSpeed=c['sound'], freestreamCp=c['cp'],
            checkpointContract=capacity, configurationSHA256=sha(folder/'case.cfg'), meshSHA256=sha(folder/'mesh.su2'),
            boundaryMethod='native prescribed supersonic exterior inlet/top and extrapolated outlet; variable-cp PT energy and EOS derivatives',
            physicsQualified=False, newQualifiedArea=0)
        atomic_write_json(folder/'sst_contract.json', sst)
        atomic_write_json(folder/'ideal_air_contract.json', contract)
        manifest = dict(requiredRuntimeFlavor=AXISYMMETRIC_FLAVOR, nodes=len(p), elements=len(cells),
            expectedVelocity=[-c['speed'], 0, 0], flight=dict(iterations=iterations, outputEvery=every),
            thermalBoundary=dict(mode='isothermal', wallTemperature=c['wallTemperature']),
            initialization=init, idealAir=contract, idealAirContractSHA256=sha(folder/'ideal_air_contract.json'),
            sst=sst, sstContractSHA256=sha(folder/'sst_contract.json'),
            publicBenchmark=dict(preparation=str(source/'preparation.json'), conditions=c, fixture=small,
                scope='attached axisymmetric cylinder-flare; zero incidence; prescribed Tw',
                directTotalTemperatureReferenceStatus='pending-reference'), physicsQualified=False)
        atomic_write_json(folder/'case_manifest.json', manifest)
        stable = {path: sha(Path(path)) == digest for path, digest in locks.items()}
        assert all(stable.values())
        record.update(status='terminal', exitCode=0, elapsedSeconds=time.monotonic()-started,
                      inputsStable=stable, initializationSHA256=contract['initializationSHA256'])
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
    """模块四。正式入口拒绝模型、维度、全状态边界、原编号或初场身份变化。"""
    record = read(folder/'axisymmetric_preparation.json')
    init_path = folder/'axisymmetric_initialization.json'
    assert record['status'] == 'terminal' and record['exitCode'] == 0 and all(record['inputsStable'].values())
    assert sha(init_path) == record['initializationSHA256'] == contract['initializationSHA256']
    init = read(init_path)
    assert init == manifest['initialization'] and init['kind'] == KIND and init['priorSuccessfulUpdates'] == 0
    assert c['AXISYMMETRIC'] == 'YES' and c['MARKER_SUPERSONIC_OUTLET'] == '(outlet)'
    assert not any(key in c for key in ('MARKER_RIEMANN', 'MARKER_GILES', 'MARKER_SYM', 'RESTART_ITER'))
    assert c['SOLUTION_FILENAME'] == 'initial_flow' and c['READ_BINARY_RESTART'] == 'NO'
    assert c['MARKER_ISOTHERMAL'] == '(wall, 311)'
    conditions = manifest['publicBenchmark']['conditions']
    assert conditions == read(Path(manifest['publicBenchmark']['preparation']))['conditions']
    assert manifest['expectedVelocity'] == [-conditions['speed'], 0, 0]
    assert manifest['sst']['kInfinity'] == conditions['k'] and manifest['sst']['omegaInfinity'] == conditions['omega']
    assert contract['freestreamSpeed'] == conditions['speed'] and contract['freestreamCp'] == conditions['cp']
    assert float(c['MACH_NUMBER']) == conditions['mach']
    assert float(c['FREESTREAM_TEMPERATURE']) == conditions['T'] and float(c['FREESTREAM_PRESSURE']) == conditions['p']
    inlet = [s.strip() for s in c['MARKER_SUPERSONIC_INLET'].strip('() ').split(',')]
    assert len(inlet) == 12 and inlet[0] == 'inlet' and inlet[6] == 'top'
    for start in (1, 7):
        assert np.array_equal(np.array(inlet[start:start+5], float), [conditions['T'], conditions['p'], -conditions['speed'], 0, 0])
    assert c['MARKER_INLET_TURBULENT'] == '(inlet, 0.00021, 0.001, top, 0.00021, 0.001)'
    locks = dict(init['sourceSHA256'])
    assert all(sha(Path(path)) == digest for path, digest in locks.items())
    assert sha(folder/'initial_flow.csv') == init['initialStateSHA256']
    values, _, nd, rho, _, thermal = state(folder/'initial_flow.csv', True)
    points, cells, markers = mesh_data(folder/'mesh.su2')
    assert nd == 2 and len(points) == manifest['nodes'] and len(cells) == manifest['elements']
    assert np.array_equal(values[:, 1:3], points) and min(points[:, 1]) > 0
    core = IdealAir()
    expected, _ = initial_state(core, points, markers, wall_geometry(folder/'mesh.su2'), conditions)
    assert np.max(abs(values[:, 3:]-expected)/np.maximum(1, abs(expected))) < 2e-12
    core.state('rhoe', rho, thermal)
    plan = contract['checkpointContract']
    assert plan == output_plan(len(points), len(np.unique(markers['wall'])), len(markers['wall']),
        int(c['ITER']), int(c['OUTPUT_WRT_FREQ']), plan['maximumOutputBytes'], extra_columns=7)
    locks.update(core.locks)
    for path in (init_path, folder/'initial_flow.csv', folder/'axisymmetric_preparation.json', Path(__file__).resolve()):
        locks[str(path)] = sha(path)
    return locks


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    parser.add_argument('--iterations', type=int, required=True)
    parser.add_argument('--every', type=int, required=True)
    parser.add_argument('--maximum-mib', type=int, required=True)
    parser.add_argument('--fixture', action='store_true')
    args = parser.parse_args()
    result = prepare(args.folder, args.iterations, args.every, args.maximum_mib*1024**2, args.fixture)
    print(result['initialization']['temperatureRangeK'])
