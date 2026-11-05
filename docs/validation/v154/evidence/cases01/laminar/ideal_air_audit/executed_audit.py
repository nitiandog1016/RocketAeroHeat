"""用固定CoolProp独立理想项核对新引擎的非紧凑保存态和非零热梯度。"""
from pathlib import Path
import argparse
import csv
import json
import re

import numpy as np
from scipy.spatial import cKDTree

from validateIdealAir import ROOT, run, save, sha, write_rows, read_rows, relative
from su2IdealAirInterface import FLAVOR, R, settings
from su2ThermalRuntime import read
from su2SSTInterface import wall_geometry, phase_columns
from auditSU2SSTState import vtk_tokens


def load_state(path, turbulent):
    """按实际具名列读取，绝对内能按OTH参考解释。"""
    with path.open(encoding='utf-8-sig') as stream:
        columns = next(csv.reader(stream))
    nd = 3 if 'z' in columns else 2
    axes = ['x', 'y', 'z'][:nd]
    expected = ['PointID']+axes+['Density']+['Momentum_'+a for a in axes]+['Energy']
    primitive = phase_columns(nd)
    if turbulent:
        expected += ['Turb_Kin_Energy', 'Omega']
    else:
        primitive = primitive[:-1]
    assert columns == expected+primitive, columns
    data = np.loadtxt(str(path), delimiter=',', skiprows=1)
    assert np.isfinite(data).all() and np.array_equal(data[:, 0], np.arange(len(data)))
    values = {name: data[:, index] for index, name in enumerate(columns)}
    assert min(values['Density']) > 0 and min(values['Pressure']) > 0
    assert np.all((values['Temperature'] >= 60) & (values['Temperature'] <= 2000))
    if turbulent:
        assert min(values['Turb_Kin_Energy']) >= 0 and min(values['Omega']) > 0
    return data, values, nd


def audit(folder):
    """分开验证流动阶段、后更新k、给定Tw和高精度热流，保留真实资格边界。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    receipt = read(folder/'engine_process.json')
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0
    assert not receipt['timedOut'] and not receipt['cancelled'] and all(receipt['inputsStable'].values())
    assert receipt['runtimeFlavor'] == FLAVOR
    contract = receipt['thermodynamicModel']
    turbulent = contract['turbulent']
    config = settings((folder/'case.cfg').read_text(encoding='utf8'))
    for name in ('case.cfg', 'mesh.su2', 'ideal_air_contract.json', 'case_manifest.json'):
        assert sha(folder/name) == receipt['inputSHA256'][str(folder/name)]
    for name in ('restart_flow.csv', 'surface_flow.vtk', 'history.csv'):
        assert sha(folder/name) == receipt['outputSHA256'][name]
    with (folder/'history.csv').open() as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    keys = [key.strip() for key in rows[0]]
    history = np.array(rows[1:], float)
    steps = history[:, keys.index('Inner_Iter')]
    assert np.array_equal(steps, np.arange(len(history))) and len(history) >= 2
    step = int(steps[-1])
    for stem, extension in (('restart_flow', 'csv'), ('surface_flow', 'vtk')):
        named = '%s_%06d.%s' % (stem, step, extension)
        assert sha(folder/named) == receipt['outputSHA256'][named] == sha(folder/(stem+'.'+extension))
    data, values, nd = load_state(folder/'restart_flow.csv', turbulent)
    axes = ['x', 'y', 'z'][:nd]
    rho = values['Density']
    velocity = np.column_stack([values['Momentum_'+axis]/rho for axis in axes])
    flow_velocity = np.column_stack([values['Velocity_'+axis] for axis in axes])
    kinetic = values['Turb_Kin_Energy'] if turbulent else np.zeros(len(data))
    previous_path = folder/('restart_flow_%06d.csv' % (step-1))
    assert sha(previous_path) == receipt['outputSHA256'][previous_path.name]
    previous, previous_values, previous_nd = load_state(previous_path, turbulent)
    assert previous_nd == nd and np.array_equal(previous[:, :nd+1], data[:, :nd+1])
    old_k = previous_values['Turb_Kin_Energy'] if turbulent else kinetic
    energy_without_k = values['Energy']/rho-.5*np.sum(velocity*velocity, axis=1)
    out = folder/'ideal_air_audit'
    out.mkdir(exist_ok=False)
    # 独立程序使用固定CoolProp的umass_idealgas/cp0mass，核对完整保存阶段。
    helper = ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    upstream = read(ROOT/'outputs/v153_validation/native03/runtime.json')
    assert sha(helper) == upstream['sha256']
    write_rows(out/'oracle_input.csv', [('rhoT', float(r), float(t)) for r, t in zip(rho, values['Temperature'])])
    write_rows(out/'current_input.csv', [('rhoe', float(r), float(e)) for r, e in zip(rho, energy_without_k-kinetic)])
    inputs = [Path(__file__).resolve(), ROOT/'scripts/su2IdealAirInterface.py', ROOT/'scripts/su2SSTInterface.py',
              ROOT/'scripts/auditSU2SSTState.py', ROOT/'scripts/validateIdealAir.py', helper, previous_path,
              folder/'restart_flow.csv', folder/'surface_flow.vtk', folder/'engine_process.json',
              out/'oracle_input.csv', out/'current_input.csv']
    lock = {str(path): sha(path) for path in inputs}
    save(out/'input_locks.json', lock)
    (out/'executed_audit.py').write_bytes(Path(__file__).read_bytes())
    for mode, option in (('oracle', '--oracle'), ('current', '--core')):
        run(out, mode, [str(helper), option, str(out/(mode+'_input.csv')), str(out/(mode+'.csv'))], 30)
    oracle = read_rows(out/'oracle.csv')
    current = read_rows(out/'current.csv')
    assert np.all(oracle['accepted'] == 1) and np.all(current['accepted'] == 1)
    errors = dict(phaseEnergy=relative(energy_without_k-old_k, oracle['e']),
                  EOS=relative(values['Pressure'], oracle['p']), velocity=relative(velocity, flow_velocity))
    assert errors['phaseEnergy'] < 2e-10 and errors['EOS'] < 2e-11 and errors['velocity'] < 2e-11, errors
    g = wall_geometry(folder/'mesh.su2')
    points = np.column_stack([data[:, 1:nd+1], np.zeros((len(data), 3-nd))])
    xyz, surface = vtk_tokens(folder/'surface_flow.vtk')
    distances, ids = cKDTree(points).query(xyz)
    assert max(distances) < 1e-6 and set(ids) == set(g['wallIDs'])
    order = {int(v): i for i, v in enumerate(g['wallIDs'])}
    ix = np.array([order[int(j)] for j in ids])
    temperature = values['Temperature']
    mu = float(config['MU_REF'])*(temperature/float(config['MU_T_REF']))**1.5
    mu *= (float(config['MU_T_REF'])+float(config['SUTHERLAND_CONSTANT']))/(temperature+float(config['SUTHERLAND_CONSTANT']))
    wall_temperature = float(re.fullmatch(r'\(\s*wall\s*,\s*([-+\deE.]+)\s*\)', config['MARKER_ISOTHERMAL'])[1])
    conductivity = mu[ids]*oracle['cp'][ids]/float(config['PRANDTL_LAM'])
    q = conductivity*(temperature[g['normalNeighbor'][ix]]-wall_temperature)/g['thermalDistance'][ix]
    errors.update(viscosity=relative(values['Laminar_Viscosity'], mu),
                  heatGradient=relative(values['Heat_Flux'][ids], q),
                  surfacePressure=relative(surface['Pressure'][:, 0], values['Pressure'][ids]),
                  surfaceTemperature=relative(surface['Temperature'][:, 0], temperature[ids]),
                  surfaceHeat=relative(surface['Heat_Flux'][:, 0], q))
    assert errors['viscosity'] < 2e-11 and errors['heatGradient'] < 2e-10, errors
    assert max(errors[name] for name in ('surfacePressure', 'surfaceTemperature', 'surfaceHeat')) < 2e-6
    assert max(abs(q)) > 1
    eddy_fraction = 0.
    if turbulent:
        eddy_fraction = float(max(abs(values['Eddy_Viscosity'][ids])/mu[ids]))
        assert eddy_fraction < 1e-12 and np.all((kinetic[ids] == 0) | (kinetic[ids] == 1e-10))
    density_inf = float(config['FREESTREAM_PRESSURE'])/(R*float(config['FREESTREAM_TEMPERATURE']))
    friction = np.column_stack([values['Skin_Friction_Coefficient_'+axis][ids] for axis in axes])
    shear = np.linalg.norm(friction, axis=1)*.5*density_inf*contract['freestreamSpeed']**2
    yplus = g['nearestWallDistance'][ix]*np.sqrt(rho[ids]*shear)/mu[ids]
    errors['yPlus'] = relative(values['Y_Plus'][ids], yplus)
    assert errors['yPlus'] < 2e-10
    power = float(g['nodeArea'][ix]@q)
    errors['power'] = relative(np.array([history[-1, keys.index('HF')]]), np.array([power]))
    assert errors['power'] < 2e-6
    residual_names = ['rms[Rho]', 'rms[RhoE]']+(['rms[k]', 'rms[w]'] if turbulent else [])
    residuals = [float(history[-1, keys.index(name)]) for name in residual_names]
    stable = {path: sha(Path(path)) == expected for path, expected in lock.items()}
    assert all(stable.values())
    report = dict(interfacePassed=True, thermodynamicModel=contract, dimension=nd, nodes=len(data),
        wallNodes=len(ids), checkpointLocalIteration=step, successfulUpdates=len(history),
        elapsedSeconds=receipt['elapsedSeconds'], errors=errors, inputSHA256=lock, inputsStable=stable,
        nativeHeatPower=power, nativeHeatPowerUnit='W/m' if nd == 2 else 'W',
        wallHeatRangeWm2=[float(min(q)), float(max(q))], givenWallTemperatureK=wall_temperature,
        wallFluidTemperatureRangeK=[float(min(temperature[ids])), float(max(temperature[ids]))],
        temperatureRangeK=[float(min(temperature)), float(max(temperature))],
        currentKTemperatureRangeK=[float(min(current['T'])), float(max(current['T']))],
        currentKvsFlowPhaseTemperatureDifferenceK=float(max(abs(current['T']-temperature))),
        flowPrimitivePhase='flow update with previous saved k', turbulencePhase='subsequent turbulence update',
        previousKCheckpoint=step-1, wallEddyToMolecularMaximum=eddy_fraction,
        finalResiduals=residuals, originalResidualThresholdPassed=step >= 10 and all(v < -8 for v in residuals),
        yPlusP99=float(np.percentile(yplus, 99)), physicsQualified=False, newQualifiedArea=0)
    table = np.column_stack([ids, xyz, values['Pressure'][ids], q, temperature[ids],
                             np.full(len(ids), wall_temperature), yplus])
    np.savetxt(str(out/'wall.csv'), table, delimiter=',', fmt='%.17g', comments='',
        header='nodeID,x_m,y_m,z_m,p_Pa,q_W_m2,fluid_T_K,given_Tw_K,yPlus')
    report['wallCSV_SHA256'] = sha(out/'wall.csv')
    save(out/'report.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    args = parser.parse_args()
    result = audit(args.folder)
    print(json.dumps({k: v for k, v in result.items() if k not in ('inputSHA256', 'inputsStable')}, indent=2))
