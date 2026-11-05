"""准备变比热公开热验证的工况、唯一网格与测站合同，不启动求解器。"""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from idealAirFreestream import total_to_static, static_to_total
from idealAirProperties import IdealAir
from prepareSU2AxisymmetricBenchmark import DATA, ROOT, audit_geometry, sha, source_blocks, write_mesh
from su2IdealAirInterface import MODEL, R
from su2SSTInterface import wall_geometry


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def profile_locations(block, contract):
    """模块一。在同一X截面逐层求交，再定位21个法向实验采样点。"""
    rows, ni = block.shape[:2]
    x = -contract['experimental_local_reference']['location_source_x_m']
    brackets, radii = [], []
    for j in range(rows):
        order = np.argsort(block[j, :, 0])
        xs = block[j, order, 0]
        assert np.all(np.diff(xs) > 0) and xs[0] < x < xs[-1]
        high = int(np.searchsorted(xs, x))
        lo, hi = int(order[high-1]), int(order[high])
        weight = (x-block[j, lo, 0])/(block[j, hi, 0]-block[j, lo, 0])
        radii.append((1-weight)*block[j, lo, 1]+weight*block[j, hi, 1])
        brackets.append((j*ni+lo, j*ni+hi, weight))
    heights = np.asarray(radii)-radii[0]
    assert np.all(np.diff(heights) > 0)
    locations = []
    for point in contract['measurements']['upstream_profile']:
        y = point['y_from_wall_m']
        assert 0 <= y <= heights[-1]
        high = max(1, min(rows-1, int(np.searchsorted(heights, y))))
        fraction = (y-heights[high-1])/(heights[high]-heights[high-1])
        lower, upper = brackets[high-1], brackets[high]
        ids = [lower[0], lower[1], upper[0], upper[1]]
        weights = [(1-fraction)*(1-lower[2]), (1-fraction)*lower[2],
                   fraction*(1-upper[2]), fraction*upper[2]]
        rebuilt = np.dot(weights, block.reshape(-1, 2)[ids])
        assert np.max(np.abs(rebuilt-[x, radii[0]+y])) < 2e-14
        row = dict(sourceRow=point['source_row'], sourceX=-x, projectX=x,
                   distanceFromWall=y, nominalMeasurementRadius=contract['geometry']['cylinder_radius_m']+y,
                   sampleRadius=radii[0]+y, inletDevelopmentLength=block[0, 0, 0]-x,
                   measuredVelocityRatio=point['velocity_ratio'],
                   measuredTotalTemperatureRatio=point['total_temperature_ratio'])
        row.update({'node%d' % i: ids[i] for i in range(4)})
        row.update({'weight%d' % i: weights[i] for i in range(4)})
        locations.append(row)
    return locations


def near_wall_estimate(old_mesh, new_mesh, old_report):
    """模块二。逐壁点保持旧场u_tau/nu，仅用真实壁距预测新网格y+。"""
    old = wall_geometry(old_mesh)
    new = wall_geometry(new_mesh)
    wall = old_report['currentFieldAssessment']['wall']
    assert np.array_equal(old['wallIDs'], wall['ids'])
    old_s = np.asarray(wall['S'])
    x = new['points'][new['wallIDs'], 0]
    new_s = np.where(x >= 0, -x, -x/np.cos(np.deg2rad(20)))
    order = np.argsort(old_s)
    assert np.all(np.diff(old_s[order]) > 0)
    assert min(new_s) >= min(old_s)-2e-14 and max(new_s) <= max(old_s)+2e-14
    inverse_scale = np.asarray(wall['yplus'])/old['nearestWallDistance']
    predicted = np.interp(new_s, old_s[order], inverse_scale[order])*new['nearestWallDistance']
    rows = [dict(node=int(i), S=float(s), distanceM=float(d), fixedOldFrictionYPlus=float(y))
            for i, s, d, y in zip(new['wallIDs'], new_s, new['nearestWallDistance'], predicted)]
    summary = dict(sourceLocalIteration=old_report['savedLocalIteration'],
                   sourceLastHistoryIteration=old_report['lastRecordedIteration'],
                   oldYPlus99=float(np.quantile(wall['yplus'], .99)),
                   estimatedYPlus99=float(np.quantile(predicted, .99)),
                   estimatedYPlusMaximum=float(max(predicted)),
                   minimumWallDistanceM=float(min(new['nearestWallDistance'])),
                   maximumWallDistanceM=float(max(new['nearestWallDistance'])),
                   targetYPlus99=1, newFlowYPlusVerified=False,
                   assumption='hold old unconverged u_tau/nu; new flow must be assessed independently')
    assert summary['estimatedYPlus99'] <= 1, '唯一网格的事前近壁估计未通过'
    return rows, summary


def prepare(folder):
    """模块三。绑定公开来源与旧场，生成可复用网格、来流和插值合同。"""
    started = time.monotonic()
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    assert not folder.exists(), '请保留既有准备记录，使用新的输出目录'
    contract_path = DATA/'contract.json'
    source_path = DATA/'source_manifest.json'
    contract = json.loads(contract_path.read_text('utf8'))
    assert sha(source_path) == contract['source_manifest_sha256']
    sources = json.loads(source_path.read_text('utf8'))
    required = {'raw/Axiswblim7_grids.zip', 'raw/case.html', 'raw/grids.html',
                'raw/sstvm.html', 'raw/wind_sstvm.inp', 'raw/wind_freestream.html'}
    locks = {str(p): sha(p) for p in (contract_path, source_path)}
    found = set()
    for source in sources['files']:
        if source['file'] in required:
            path = DATA/source['file']
            assert sha(path) == source['sha256'], str(path)
            locks[str(path)] = source['sha256']
            found.add(source['file'])
    assert found == required
    old_folder = ROOT/'outputs/v152_validation/acc_sbli01'
    old_report_path = ROOT/'outputs/v152_validation/review01/report.json'
    for path in (old_folder/'mesh.su2', old_report_path):
        locks[str(path)] = sha(path)
    old_report = json.loads(old_report_path.read_text('utf8'))
    assert old_report['savedLocalIteration'] == 5000 and old_report['physicsQualified'] is False
    model = IdealAir()
    locks.update(model.locks)
    for name in ('prepareSU2PublicValidation.py', 'idealAirFreestream.py',
                 'prepareSU2AxisymmetricBenchmark.py', 'su2SSTInterface.py', 'su2IdealAirInterface.py'):
        path = ROOT/'scripts'/name
        locks[str(path)] = sha(path)

    c = contract['computational_truncated_inflow']
    value = total_to_static(model, c['Mach'], c['total_temperature_K_WIND_input'],
                            c['total_pressure_Pa_WIND_input'])
    mu = 1.716e-5*(value['T']/273.15)**1.5*383.55/(value['T']+110.4)
    k = 1.5*(value['speed']*.00021)**2
    value.update(gasConstant=R, wallTemperature=311, turbulenceIntensity=.00021,
                 eddyToLaminarRatio=.001, viscosity=mu, unitReynolds=value['rho']*value['speed']/mu,
                 k=k, omega=value['rho']*k/(mu*.001), conservativeSpecificEnergy=value['e']+.5*value['speed']**2+k,
                 velocity=[-value['speed'], 0, 0],
                 totalStateDefinition='mean-flow thermodynamic total state; SST k added separately')
    reference = contract['experimental_local_reference']
    inferred = static_to_total(model, reference['pressure_Pa'], reference['temperature_K'], reference['velocity_m_s'])
    cp_old = 1.4*287.058/.4
    normalization = dict(experimentalTotalTemperatureK=reference['total_temperature_K'],
        derivedVariableCpTotalTemperatureK=inferred['totalTemperature'],
        derivedConstantCpTotalTemperatureK=reference['temperature_K']+reference['velocity_m_s']**2/(2*cp_old),
        referenceStatus='inferred from independently rounded local T and U; measured total-temperature reference unresolved',
        directExperimentalTotalTemperatureGate='pending-reference',
        velocityReferenceMS=reference['velocity_m_s'], pressureReferencePa=reference['pressure_Pa'],
        heatingReferenceWM2=reference['reported_wall_heating_magnitude_W_m2'])
    assert normalization['experimentalTotalTemperatureK'] is None

    folder.mkdir(parents=True)
    source = write_mesh(folder, level=1)
    geometry = audit_geometry(folder/'mesh.su2', expected_nodes=64521, expected_cells=64000)
    blocks, _ = source_blocks(level=1)
    old_blocks, _ = source_blocks()
    for current, previous in zip(blocks, old_blocks):
        assert np.array_equal(current[::2, ::2], previous), '公开嵌套网格身份不符'
    locations = profile_locations(blocks[0], contract)
    write_csv(folder/'profile_locations.csv', locations)
    wall_rows, near_wall = near_wall_estimate(old_folder/'mesh.su2', folder/'mesh.su2', old_report)
    write_csv(folder/'near_wall_estimates.csv', wall_rows)
    original = folder/'default_compatibility'
    original.mkdir()
    write_mesh(original)
    assert sha(original/'mesh.su2') == sha(old_folder/'mesh.su2'), '旧默认网格输出逐字回归失败'
    old_conditions = old_report['conditions']
    differences = {key: value[key]/old_conditions[old_key]-1 for key, old_key in
                   [('T', 'temperatureK'), ('p', 'pressurePa'), ('speed', 'speedMS'),
                    ('rho', 'densityKgM3'), ('unitReynolds', 'unitReynoldsPerM')]}
    source_stable = all(sha(Path(path)) == expected for path, expected in locks.items())
    assert source_stable
    output = dict(schema='RocketAeroHeat.public-validation-preparation.v1',
        prepared=True, physicsQualified=False, newQualifiedArea=0, solverRuns=0,
        thermodynamicModel=MODEL, energyEntropyReference='OTH', conditions=value,
        inputChoice='WIND Mach and total temperature/pressure; mean-flow total enthalpy and entropy',
        relativeChangeFromOldConstantCp=differences,
        unitReynoldsRelativeToFigure=value['unitReynolds']/c['unit_Re_m_inv_figure']-1,
        normalization=normalization, grid=source, geometry=geometry, nearWall=near_wall,
        profile=dict(axialStationCount=1, normalSampleCount=len(locations), projectX=.06,
                     nominalInletX=c['project_X_m'], actualInletX=float(blocks[0][0, 0, 0]),
                     actualDevelopmentLengthM=float(blocks[0][0, 0, 0]-.06),
                     nominalRadiusM=contract['geometry']['cylinder_radius_m'],
                     actualWallRadiusM=float(blocks[0][0, 0, 1])),
        nextModelContract=dict(model='SST-Vm', options=['V1994m', 'VORTICITY'],
            implemented=False, dimension=2, axisymmetric=True, intensity=.00021, viscosityRatio=.001,
            reason='existing variable-cp phase-v2 contract validates only 3D SST-2003m'),
        defaultMeshByteIdentical=True, sourceSHA256=locks, inputsStable=source_stable,
        files={p.name: dict(bytes=p.stat().st_size, sha256=sha(p))
               for p in folder.iterdir() if p.is_file()}, elapsedSeconds=time.monotonic()-started)
    (folder/'preparation.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf8')
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True)
    result = prepare(parser.parse_args().out)
    print(json.dumps({key: result[key] for key in ('prepared', 'elapsedSeconds', 'geometry', 'nearWall')}))
