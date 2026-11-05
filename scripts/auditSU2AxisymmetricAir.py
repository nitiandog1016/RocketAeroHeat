"""公开变比热轴对称保存场的独立热状态、旋转守恒和26项验收。"""
import argparse
import csv
from pathlib import Path

import numpy as np

from auditSU2IdealAirState import load_state
from idealAirProperties import IdealAir
from runEmbeddedSU2 import atomic_write_json
from su2AxisymmetricAnalysis import (mesh_data, boundary_balance, residual_scales,
    wall_values, sample_stations, integrate_scalar)
from su2BenchmarkAcceptance import acceptance, stability
from su2SSTInterface import wall_geometry
from su2ThermalRuntime import ROOT, read, sha


def serial(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: serial(item) for key, item in value.items()}
    if isinstance(value, list):
        return [serial(item) for item in value]
    return value


def assess(folder, entry, model, conditions, geometry, points, cells, markers, contract, locations):
    """模块一。每帧使用直接捕获的k核对流动阶段，所有积分保持实际OTH参考。"""
    path = folder/entry['volume']
    assert sha(path) == entry['volumeSHA256']
    data, fields, nd = load_state(path, True, True, True)
    assert nd == 2 and np.allclose(data[:, 1:3], points, rtol=2e-14, atol=1e-14)
    f = fields
    rho = f['Density']
    u = np.column_stack((f['Momentum_x'], f['Momentum_y']))/rho[:, None]
    f['velocity'] = u
    thermo = model.state('PT', f['Pressure'], f['Temperature'])
    eos_error = float(max(abs(rho/thermo['rho']-1)))
    thermal = f['Energy']/rho-.5*np.sum(u**2, axis=1)-f['Flow_Primitive_K']
    phase_error = float(max(abs(thermal-thermo['e'])/np.maximum(1, abs(thermo['e']))))
    ids = geometry['wallIDs']
    heat = f['Laminar_Viscosity'][ids]*thermo['cp'][ids]/.72
    heat *= (f['Temperature'][geometry['normalNeighbor']]-conditions['wallTemperature'])/geometry['thermalDistance']
    heat_error = float(max(abs(heat-f['Heat_Flux'][ids])/np.maximum(abs(heat), 1)))
    assert eos_error < 2e-12 and phase_error < 2e-10 and heat_error < 2e-6
    wall = wall_values(points, markers, f)
    stations = sample_stations(wall, contract)
    indices = np.array([[int(row['node%d' % i]) for i in range(4)] for row in locations])
    weights = np.array([[float(row['weight%d' % i]) for i in range(4)] for row in locations])
    sampled_u = np.sum(weights*(-u[indices, 0]), axis=1)/1274.
    mean_total_h = thermo['h']+.5*np.sum(u**2, axis=1)
    # 总温反解保留温区审核。实验分母缺失，推导参考只用于诊断。
    total_t = model.state('Ph', f['Pressure'], mean_total_h)['T']
    inferred_reference = read(Path(read(folder/'case_manifest.json')['publicBenchmark']['preparation']))['normalization']
    inferred_t0 = inferred_reference['derivedVariableCpTotalTemperatureK']
    sampled_t = np.sum(weights*total_t[indices], axis=1)/inferred_t0
    measured_u = np.array([float(row['measuredVelocityRatio']) for row in locations])
    measured_t = np.array([float(row['measuredTotalTemperatureRatio']) for row in locations])
    profile = dict(y=np.array([float(row['distanceFromWall']) for row in locations]),
        predictedU=sampled_u, measuredU=measured_u, inferredReferenceT0K=inferred_t0,
        predictedTotalTemperatureRatioWithInferredReference=sampled_t, measuredTotalTemperatureRatio=measured_t,
        velocityRatioRMSError=float(np.sqrt(np.mean((sampled_u-measured_u)**2))),
        inferredTotalTemperatureRatioRMSError=float(np.sqrt(np.mean((sampled_t-measured_t)**2))),
        totalTemperatureRatioRMSError=None, totalTemperatureReferenceStatus='pending experimental denominator')
    c = conditions
    scales = dict(speedMS=c['speed'], densityKgM3=c['rho'], temperatureK=c['T'], turbulenceIntensity=c['turbulenceIntensity'])
    balance = boundary_balance(points, cells, markers, f, scales, specific_heat=thermo['cp'])
    residual = residual_scales(points, cells, f, scales, total_enthalpy=c['totalEnthalpy']+c['k'])
    residual['reference'] = 'OTH mean total enthalpy plus SST k; physical control-area RHS scaling'
    power = integrate_scalar(points, markers['wall'], f['Heat_Flux'])
    peak = stations['heat']
    ip, im = int(np.argmax(peak['predicted'])), int(np.argmax(peak['measured']))
    return dict(EOSRelativeError=eos_error, phaseEnergyRelativeError=phase_error, heatGradientRelativeError=heat_error,
        wall=wall, stations=stations, profile=profile, boundaryBalance=balance, residuals=residual,
        revolvedHeatPowerW=power, yPlus99=float(np.quantile(wall['yplus'], .99)),
        maximumWallFluidTemperatureDifferenceK=float(max(abs(wall['fluidT']-c['wallTemperature']))),
        savedWallFluidTemperatureRangeK=[float(min(wall['fluidT'])), float(max(wall['fluidT']))],
        localReference=dict(pressurePa=float(np.interp(-.06, wall['S'], wall['p'])),
            heatFluxWm2=float(np.interp(-.06, wall['S'], wall['q']))),
        measuredStationHeatPeak=dict(relativeMagnitudeError=float(abs(peak['predicted'][ip]/peak['measured'][im]-1)),
            positionErrorM=float(abs(peak['S'][ip]-peak['S'][im])), predictedWm2=float(peak['predicted'][ip]),
            measuredWm2=float(peak['measured'][im])))


def audit(folder, output=None):
    """模块二。仅分析已完整保存帧，原任务终态、缺失帧和资格分别记录。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    out = Path(output).resolve() if output else folder/'public_review'
    out.relative_to(folder)
    out.mkdir(exist_ok=False)
    manifest = read(folder/'case_manifest.json')
    receipt = read(folder/'engine_process.json')
    assert not manifest['publicBenchmark']['fixture']
    assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
    for name in ('mesh.su2', 'case.cfg', 'case_manifest.json', 'ideal_air_contract.json'):
        assert sha(folder/name) == receipt['inputSHA256'][str(folder/name)]
    source = Path(manifest['publicBenchmark']['preparation']).parent
    conditions = manifest['publicBenchmark']['conditions']
    contract = read(ROOT/'tests/data/tmr_mach7/contract.json')
    with (source/'profile_locations.csv').open(encoding='utf8') as stream:
        locations = list(csv.DictReader(stream))
    points, cells, markers = mesh_data(folder/'mesh.su2')
    geometry = wall_geometry(folder/'mesh.su2')
    model = IdealAir()
    dependencies = [Path(__file__).resolve(), ROOT/'scripts/su2AxisymmetricAnalysis.py',
        ROOT/'scripts/su2BenchmarkAcceptance.py', ROOT/'scripts/auditSU2IdealAirState.py',
        source/'profile_locations.csv', folder/'engine_process.json']
    locks = {str(p): sha(p) for p in dependencies}
    locks.update(model.locks)
    frames = []
    current = None
    for entry in receipt['checkpoints']['complete']:
        current = assess(folder, entry, model, conditions, geometry, points, cells, markers, contract, locations)
        step = entry['localIteration']
        atomic_write_json(out/('frame_%06d.json' % step), serial(current))
        frames.append(dict(iteration=step, p=current['wall']['p'], q=current['wall']['q'], power=current['revolvedHeatPowerW']))
        locks[str(folder/entry['volume'])] = entry['volumeSHA256']
    with (folder/'history.csv').open(encoding='utf8') as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    keys = [key.strip() for key in rows[0]]
    finite = [list(map(float, row)) for row in rows[1:] if len(row) == len(keys)]
    finite = [row for row in finite if np.isfinite(row).all()]
    last = finite[-1] if finite else None
    last_iteration = int(last[keys.index('Inner_Iter')]) if last else None
    residuals = {name: float(last[keys.index(name)]) if last else None for name in ('rms[Rho]', 'rms[RhoE]', 'rms[k]', 'rms[w]')}
    steady = stability(frames, last_iteration, integrate_scalar(points, markers['wall'], np.ones(len(points))))
    normal = receipt['exitCode'] == 0 and not any(receipt.get(k, False) for k in
        ('timedOut', 'cancelled', 'supervisorError', 'outputBudgetExceeded'))
    same_final = bool(frames and frames[-1]['iteration'] == last_iteration)
    accepted = acceptance(current, residuals, steady, normal, same_final,
        missing_reasons={'upstream total temperature profile': 'experimental total-temperature denominator unresolved'})
    stable = {path: sha(Path(path)) == digest for path, digest in locks.items()}
    assert all(stable.values())
    report = dict(schema='RocketAeroHeat.axisymmetric-ideal-air-review.v1', nativeExitCode=receipt['exitCode'],
        nativeNormalTerminal=normal, elapsedSeconds=receipt['elapsedSeconds'], lastRecordedIteration=last_iteration,
        savedLocalIteration=frames[-1]['iteration'] if frames else None, completedUpdates=len(finite),
        currentFieldAssessment=current, finalRecordedResiduals=residuals, acceptance=accepted, temporalStability=steady,
        sourceSHA256=locks, inputsStable=stable, model=manifest['sst']['model'], conditions=conditions,
        abnormalInstantStateAvailable=False, physicsQualified=False, newQualifiedArea=0)
    atomic_write_json(out/'report.json', serial(report))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output')
    args = parser.parse_args()
    report = audit(args.folder, args.output)
    print(report['acceptance'])
