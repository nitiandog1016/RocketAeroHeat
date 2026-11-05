"""一般三维混合域的新参考SST初场，复用三角壁距并以真实变比热焓反解温度。"""
from pathlib import Path
import json
import re
import time

import numpy as np

from idealAirProperties import IdealAir
from prepareSU2WallDistanceState import surface_distance, boundary_faces
from su2SSTInterface import wall_geometry, mesh, state
from su2IdealAirInterface import R, FLAVOR, settings
from su2ThermalRuntime import read, sha
from runEmbeddedSU2 import ROOT, atomic_write_json

KIND = 'general-mesh-ideal-air-SST-wall-distance-v1'


def assemble(core, detail, faces, thickness, velocity_inf, pressure, temperature_inf,
             wall_temperature, kinetic_inf, omega_inf, nu_wall):
    """模块一。用同一热力学核构造线性混合总焓及正SST量，壁面/远场精确赋值。"""
    distance = detail['distance']
    assert np.isfinite(thickness) and thickness > 0
    assert np.all((distance >= 0) & (distance <= thickness))
    eta = distance/thickness
    blend = eta*eta*(3-2*eta)
    wall_ids, far_ids = detail['wallIDs'], detail['farIDs']
    assert np.all(blend[wall_ids] == 0) and np.all(blend[far_ids] == 1)
    endpoints = core.state('PT', pressure, [wall_temperature, temperature_inf])
    speed2 = float(velocity_inf @ velocity_inf)
    target = (1-blend)*endpoints['h'][0]+blend*endpoints['h'][1]+.5*speed2*blend*(1-blend)
    thermo = core.state('Ph', pressure, target)
    # 端点直接复用同一PT状态，避免反解舍入使给定边界产生微小偏移。
    for name, values in thermo.items():
        values[blend == 0] = endpoints[name][0]
        values[blend == 1] = endpoints[name][1]
    velocity = blend[:, None]*velocity_inf
    kinetic = blend*kinetic_inf
    omega = np.full(len(distance), omega_inf)
    inside = detail['nearestFace'] >= 0
    lengths = np.full(len(distance), np.nan)
    lengths[wall_ids] = detail['wallLengths']
    local_length = np.sum(lengths[faces[detail['nearestFace'][inside]]]*detail['weights'][inside], axis=1)
    assert np.all(local_length > 0) and np.isfinite(local_length).all()
    near = 60*nu_wall/(.075*(distance[inside]+local_length)**2)
    omega[inside] = np.exp((1-blend[inside])*np.log(near)+blend[inside]*np.log(omega_inf))
    omega[wall_ids] = 60*nu_wall/(.075*detail['wallLengths']**2)
    rho = thermo['rho']
    energy = rho*(thermo['e']+.5*np.sum(velocity**2, axis=1)+kinetic)
    fields = np.column_stack([rho, rho[:, None]*velocity, energy, kinetic, omega])
    assert np.isfinite(fields).all() and min(rho) > 0 and min(omega) > 0 and min(kinetic) >= 0
    total_h = thermo['h']+.5*np.sum(velocity**2, axis=1)+kinetic
    expected_h = (1-blend)*endpoints['h'][0]+blend*(endpoints['h'][1]+.5*speed2+kinetic_inf)
    error = float(np.max(abs(total_h-expected_h)/np.maximum(1., abs(expected_h))))
    assert error < 2e-11, error
    return fields, thermo, error


def wall_viscosity(c, temperature, pressure):
    """与配置中的Sutherland参数一致，返回壁面运动黏度供omega边界初猜使用。"""
    mu = float(c['MU_REF'])*(temperature/float(c['MU_T_REF']))**1.5
    mu *= (float(c['MU_T_REF'])+float(c['SUTHERLAND_CONSTANT']))/(temperature+float(c['SUTHERLAND_CONSTANT']))
    return mu/(pressure/(R*temperature))


def prepare_wall_distance(folder, spec, provenance):
    """模块二。独占准备并保留原始输入，全部状态检查通过后发布新配置。"""
    from prepareSU2IdealAirCase import configuration, publish
    from prepareSU2SSTCase import set_value
    from su2IdealAirInterface import validate
    start = time.monotonic()
    marker = folder/'ideal_air_preparation.json'
    with marker.open('x', encoding='utf8') as stream:
        json.dump(dict(status='preparing', exitCode=None), stream)
    record = dict(status='preparing', exitCode=None, physicsQualified=False)
    try:
        assert set(spec) == {'thicknessM', 'intensity', 'viscosityRatio', 'Prt'}
        thickness = float(spec['thicknessM'])
        assert np.isfinite(thickness) and thickness > 0 and provenance.strip()
        original = folder/'ideal_air_original_inputs'
        original.mkdir(exist_ok=False)
        for name in ('case.cfg', 'case_manifest.json'):
            (original/name).write_bytes((folder/name).read_bytes())
        text = (original/'case.cfg').read_text(encoding='utf-8-sig')
        c = settings(text)
        manifest = read(original/'case_manifest.json')
        assert not any(key in manifest for key in ('sst', 'initialization', 'continuation', 'thermalProfileSHA256'))
        assert manifest['thermalBoundary']['mode'] == 'isothermal'
        assert c['MARKER_FAR'] == '(farfield)' and 'MARKER_RIEMANN' not in c
        wall = re.fullmatch(r'\(\s*wall\s*,\s*([-+\deE.]+)\s*\)', c['MARKER_ISOTHERMAL'])
        assert wall
        wall_temperature = float(wall[1])
        pressure, temperature = float(c['FREESTREAM_PRESSURE']), float(c['FREESTREAM_TEMPERATURE'])
        mach = float(c['MACH_NUMBER'])
        assert np.isfinite(mach) and mach > 1
        flight = manifest['flight']
        assert flight['Mach'] == mach and flight['temperature'] == temperature and flight['pressure'] == pressure
        assert flight['wallTemperature'] == wall_temperature
        alpha, beta = np.deg2rad([float(c['AOA']), float(c['SIDESLIP_ANGLE'])])
        direction = np.array([np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)])
        project_alpha, project_beta = np.deg2rad([flight['alpha'], flight['beta']])
        expected_direction = np.array([-np.cos(project_alpha)*np.cos(project_beta),
                                       -np.sin(project_beta), np.sin(project_alpha)*np.cos(project_beta)])
        assert np.max(abs(direction-expected_direction)) < 2e-15
        old_speed = mach*np.sqrt(float(c['GAMMA_VALUE'])*float(c['GAS_CONSTANT'])*temperature)
        assert np.max(abs(np.asarray(manifest['expectedVelocity'])-old_speed*direction)) < 2e-11
        core = IdealAir()
        freestream = core.state('PT', pressure, temperature)
        core.state('PT', pressure, wall_temperature)
        text, contract, sst = configuration(text, True, provenance, freestream, spec)
        velocity = contract['freestreamSpeed']*direction
        text = re.sub(r'^\s*MARKER_FAR\s*=.*\n?', '', text, flags=re.M)
        boundary = '(farfield, STATIC_SUPERSONIC_INFLOW_PT, %.17g, %.17g, %.17g, %.17g, %.17g)'
        changes = dict(MARKER_RIEMANN=boundary % (pressure, temperature, *(mach*direction)),
                       RESTART_SOL='YES', SOLUTION_FILENAME='initial_flow', READ_BINARY_RESTART='NO')
        for key, value in changes.items():
            text = set_value(text, key, value)
        # 源锁在壁距搜索之前建立，长准备期间的修改会阻止发布。
        sources = [Path(__file__).resolve(), ROOT/'scripts/prepareSU2IdealAirCase.py',
                   ROOT/'scripts/prepareSU2WallDistanceState.py', ROOT/'scripts/su2SSTInterface.py',
                   ROOT/'scripts/prepareSU2SSTCase.py', ROOT/'scripts/su2ThermalRuntime.py',
                   ROOT/'scripts/su2IdealAirInterface.py', folder/'mesh.su2',
                   original/'case.cfg', original/'case_manifest.json']
        locks = {str(path): sha(path) for path in sources}
        locks.update(core.locks)
        g = wall_geometry(folder/'mesh.su2')
        points = g['points']
        assert points.shape == (manifest['nodes'], 3)
        boundary_faces_by_name = boundary_faces(folder/'mesh.su2')
        faces = boundary_faces_by_name['wall']
        assert len(faces) == manifest['originalSurfaceFaceCount']
        distance, nearest, weights, pairs = surface_distance(points, points, faces, thickness)
        assert np.all(distance[g['wallIDs']] < 1e-12)
        distance[g['wallIDs']] = 0.
        detail = dict(distance=distance, nearestFace=nearest, weights=weights, wallIDs=g['wallIDs'],
                      wallLengths=g['nearestWallDistance'], farIDs=np.unique(boundary_faces_by_name['farfield']))
        fields, thermo, h_error = assemble(core, detail, faces, thickness, velocity, pressure, temperature,
            wall_temperature, sst['kInfinity'], sst['omegaInfinity'], wall_viscosity(c, wall_temperature, pressure))
        pending = folder/'initial_flow.pending'
        keys = ['PointID', 'x', 'y', 'z', 'Density', 'Momentum_x', 'Momentum_y', 'Momentum_z',
                'Energy', 'Turb_Kin_Energy', 'Omega']
        np.savetxt(str(pending), np.column_stack([np.arange(len(points)), points, fields]), delimiter=',',
                   fmt='%.17g', header=','.join('"'+key+'"' for key in keys), comments='')
        _, _, _, density, _, thermal_energy = state(pending, True)
        recovered = core.state('rhoe', density, thermal_energy)
        errors = dict(temperature=float(max(abs(recovered['T']/thermo['T']-1))),
                      pressure=float(max(abs(recovered['p']/pressure-1))), totalEnthalpy=h_error)
        assert max(errors.values()) < 2e-11, errors
        stable = {path: sha(Path(path)) == digest for path, digest in locks.items()}
        assert all(stable.values())
        pending.replace(folder/'initial_flow.csv')
        geometry = folder/'ideal_air_wall_distance.npz'
        np.savez_compressed(str(geometry), **detail)
        locks[str(geometry)] = sha(geometry)
        init = dict(kind=KIND, initialGuessOnly=True, oldFlowInterpolated=False,
            localIterationStartsAt=0, sourceLocalIteration=None, thicknessM=thickness,
            nodes=len(points), originalFaces=len(faces), wallNodes=len(g['wallIDs']),
            givenWallTemperatureK=wall_temperature, freestreamConvention='preserve-Mach-static-PT-direction',
            expectedVelocity=velocity.tolist(), originalExpectedVelocity=manifest['expectedVelocity'],
            initialStateSHA256=sha(folder/'initial_flow.csv'), sourceSHA256=locks,
            errors=errors, temperatureRangeK=[float(min(thermo['T'])), float(max(thermo['T']))],
            distanceTrianglePairs=pairs, physicsQualified=False, newQualifiedArea=0,
            method='exact triangle distance; C1 velocity/total enthalpy; ideal-air h inverse; positive SST')
        atomic_write_json(folder/'ideal_air_initialization.json', init)
        manifest.update(initialization=init, expectedVelocity=velocity.tolist(),
                        thermalPropertyMode='ideal-air variable cp; OTH energy; SST-2003m; Sutherland')
        contract.update(schemaVersion=2, initialization=KIND,
                        initializationSHA256=sha(folder/'ideal_air_initialization.json'))
        publish(folder, text, manifest, contract, sst)
        record.update(status='terminal', exitCode=0, elapsedSeconds=time.monotonic()-start,
                      inputsStable=stable, initializationSHA256=contract['initializationSHA256'])
        atomic_write_json(marker, record)
        runtime = read(ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json')
        validate(folder, manifest, text, runtime)
        return contract
    except Exception as error:
        record.update(status='terminal', exitCode=1, elapsedSeconds=time.monotonic()-start,
                      error=type(error).__name__+': '+str(error))
        atomic_write_json(marker, record)
        raise


def validate_source(folder, c, manifest, contract):
    """模块三。启动前按原编号核对写出守恒态、几何和热边界，拒绝部分发布。"""
    preparation = read(folder/'ideal_air_preparation.json')
    assert preparation['status'] == 'terminal' and preparation['exitCode'] == 0
    assert all(preparation['inputsStable'].values())
    path = folder/'ideal_air_initialization.json'
    assert sha(path) == contract['initializationSHA256'] == preparation['initializationSHA256']
    init = read(path)
    assert init == manifest['initialization'] and init['kind'] == KIND
    assert init['initialGuessOnly'] and not init['oldFlowInterpolated'] and not init['physicsQualified']
    assert init['localIterationStartsAt'] == 0 and init['sourceLocalIteration'] is None
    assert init['freestreamConvention'] == 'preserve-Mach-static-PT-direction'
    for source, digest in init['sourceSHA256'].items():
        assert sha(Path(source)) == digest, source
    initial = folder/'initial_flow.csv'
    assert sha(initial) == init['initialStateSHA256']
    values, _, dimension, density, velocity, energy = state(initial, True)
    points, _, wall = mesh(folder/'mesh.su2')
    assert dimension == 3 and values.shape[0] == init['nodes'] and np.array_equal(values[:, 1:4], points)
    assert np.array_equal(manifest['expectedVelocity'], init['expectedVelocity'])
    core = IdealAir()
    for source, digest in core.locks.items():
        assert init['sourceSHA256'][source] == digest
    pressure, temperature = float(c['FREESTREAM_PRESSURE']), float(c['FREESTREAM_TEMPERATURE'])
    reference = core.state('PT', pressure, temperature)
    expected_speed = float(c['MACH_NUMBER'])*reference['sound'][0]
    u_inf = np.asarray(init['expectedVelocity'])
    assert abs(np.linalg.norm(u_inf)/expected_speed-1) < 2e-14
    assert abs(contract['freestreamSpeed']/expected_speed-1) < 2e-14
    assert abs(contract['freestreamCp']/reference['cp'][0]-1) < 2e-14
    model = manifest['sst']
    expected_k = 1.5*(expected_speed*float(c['FREESTREAM_TURBULENCEINTENSITY']))**2
    expected_omega = expected_k/(wall_viscosity(c, temperature, pressure)*float(c['FREESTREAM_TURB2LAMVISCRATIO']))
    assert abs(model['kInfinity']/expected_k-1) < 2e-14
    assert abs(model['omegaInfinity']/expected_omega-1) < 2e-14
    boundary = [v.strip() for v in c['MARKER_RIEMANN'].strip('() ').split(',')]
    assert boundary[:2] == ['farfield', 'STATIC_SUPERSONIC_INFLOW_PT'] and len(boundary) == 7
    assert float(boundary[2]) == pressure and float(boundary[3]) == temperature
    assert np.max(abs(np.array(boundary[4:], float)*reference['sound'][0]-u_inf)) < 2e-11
    wall_temperature = init['givenWallTemperatureK']
    wall_value = re.fullmatch(r'\(\s*wall\s*,\s*([-+\deE.]+)\s*\)', c['MARKER_ISOTHERMAL'])
    assert wall_value and float(wall_value[1]) == wall_temperature
    with np.load(str(folder/'ideal_air_wall_distance.npz')) as saved:
        detail = dict(saved)
    faces = np.array([row[1:] for row in wall], np.int64)
    expected, thermo, _ = assemble(core, detail, faces, init['thicknessM'], u_inf, pressure, temperature,
        wall_temperature, model['kInfinity'], model['omegaInfinity'], wall_viscosity(c, wall_temperature, pressure))
    assert np.max(abs(values[:, 4:]-expected)/np.maximum(1., abs(expected))) < 2e-12
    recovered = core.state('rhoe', density, energy)
    assert max(abs(recovered['T']/thermo['T']-1)) < 2e-11 and max(abs(recovered['p']/pressure-1)) < 2e-11
    locks = dict(init['sourceSHA256'])
    locks.update({str(path): sha(path), str(initial): sha(initial), str(folder/'ideal_air_preparation.json'): sha(folder/'ideal_air_preparation.json')})
    return locks
