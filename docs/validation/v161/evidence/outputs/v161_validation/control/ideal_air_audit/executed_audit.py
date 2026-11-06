"""用固定CoolProp独立理想项核对新引擎的非紧凑保存态和非零热梯度。"""
from pathlib import Path
import argparse
import csv
import json
import re

import numpy as np
from scipy.spatial import cKDTree

from validateIdealAir import ROOT, run, save, sha, write_rows, read_rows, relative
from su2IdealAirInterface import FLAVOR, PHASE_FLAVORS, R, settings
from su2ThermalRuntime import read
from su2SSTInterface import wall_geometry, phase_columns, spatial_columns
from auditSU2SSTState import vtk_tokens


def load_state(path, turbulent, captured=False, spatial=False):
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
    assert not spatial or (turbulent and captured and nd == 2)
    assert columns == expected+primitive+(spatial_columns(nd)[:nd+4] if spatial else [])+(
        ['Flow_Primitive_K'] if captured else []), columns
    data = np.loadtxt(str(path), delimiter=',', skiprows=1)
    assert np.isfinite(data).all() and np.array_equal(data[:, 0], np.arange(len(data)))
    values = {name: data[:, index] for index, name in enumerate(columns)}
    assert min(values['Density']) > 0 and min(values['Pressure']) > 0
    assert np.all((values['Temperature'] >= 60) & (values['Temperature'] <= 2000))
    if turbulent:
        assert min(values['Turb_Kin_Energy']) >= 0 and min(values['Omega']) > 0
    if captured:
        assert min(values['Flow_Primitive_K']) >= 0
    return data, values, nd


def publish_review(folder, audit_directory, checkpoint=None):
    """将已完成的审核绑定到工作台；不重新计算，也不覆盖原审核身份。"""
    folder = Path(folder).resolve()
    directory = Path(audit_directory).resolve()
    directory.relative_to(folder)
    report = read(directory/'report.json')
    assert report['interfacePassed'] and not report['physicsQualified'] and all(report['inputsStable'].values())
    assert sha(folder/'engine_process.json') == report['inputSHA256'][str(folder/'engine_process.json')]
    assert sha(directory/'wall.csv') == report['wallCSV_SHA256']
    explicit = {}
    if checkpoint is not None:
        receipt = read(folder/'engine_process.json')
        assert not report['nativeNormalTerminal'] and receipt['status'] == 'terminal'
        assert receipt['runtimeFlavor'] == 'rah-ideal-air-phase-v2'
        entry = receipt['checkpoints']['lastComplete']
        assert checkpoint == report['checkpointLocalIteration'] == entry['localIteration']
        for name, suffix in (('restart_flow', 'csv'), ('surface_flow', 'vtk')):
            base = name+'.'+suffix
            selected = '%s_%06d.%s' % (name, checkpoint, suffix)
            assert sha(folder/base) == sha(folder/selected) == receipt['outputSHA256'][selected] == receipt['outputSHA256'][base]
        explicit = dict(explicitAbnormalCheckpoint=checkpoint, abnormalInstantStateAvailable=False)
    target = folder/'ideal_air_review.json'
    assert not target.exists()
    save(target, dict(auditPath=(directory/'report.json').relative_to(folder).as_posix(),
        auditSHA256=sha(directory/'report.json'), wallPath=(directory/'wall.csv').relative_to(folder).as_posix(),
        wallSHA256=sha(directory/'wall.csv'), thermoPath=(directory/'oracle.csv').relative_to(folder).as_posix(),
        thermoSHA256=sha(directory/'oracle.csv'), receiptSHA256=sha(folder/'engine_process.json'), **explicit))


def audit(folder, output=None, reuse_properties=None, checkpoint=None, recheck_checkpoints=False):
    """分开验证流动阶段、后更新k、给定Tw和高精度热流，保留真实资格边界。"""
    folder = Path(folder).resolve()
    folder.relative_to(ROOT/'outputs')
    receipt = read(folder/'engine_process.json')
    assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
    captured = receipt['runtimeFlavor'] in PHASE_FLAVORS
    assert receipt['runtimeFlavor'] in (FLAVOR, *PHASE_FLAVORS)
    normal = receipt['exitCode'] == 0 and not any(receipt.get(key, False)
        for key in ('timedOut', 'cancelled', 'outputBudgetExceeded', 'supervisorError'))
    assert normal or (captured and checkpoint is not None), '异常终态须显式选择已完整保存的检查点'
    contract = receipt['thermodynamicModel']
    turbulent = contract['turbulent']
    config = settings((folder/'case.cfg').read_text(encoding='utf8'))
    for name in ('case.cfg', 'mesh.su2', 'ideal_air_contract.json', 'case_manifest.json'):
        assert sha(folder/name) == receipt['inputSHA256'][str(folder/name)]
    for name in (('history.csv',) if captured else ('restart_flow.csv', 'surface_flow.vtk', 'history.csv')):
        assert sha(folder/name) == receipt['outputSHA256'][name]
    with (folder/'history.csv').open() as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    keys = [key.strip() for key in rows[0]]
    history = np.array(rows[1:], float)
    steps = history[:, keys.index('Inner_Iter')]
    assert np.array_equal(steps, np.arange(len(history))) and len(history) >= (1 if captured else 2)
    step = int(steps[-1]) if checkpoint is None else checkpoint
    assert isinstance(step, int) and 0 <= step < len(history)
    assert np.isfinite(history[step]).all()
    checked_inventory = None
    if recheck_checkpoints:
        assert captured and checkpoint is not None, '只读复核必须显式指定具名帧'
        from su2IdealAirCheckpoints import inventory
        checked_inventory = inventory(folder, receipt)
    if captured:
        selected_inventory = checked_inventory if recheck_checkpoints else receipt['checkpoints']
        entry = next(row for row in selected_inventory['complete'] if row['localIteration'] == step)
        volume_path, surface_path = folder/entry['volume'], folder/entry['surface']
        assert sha(volume_path) == entry['volumeSHA256'] == receipt['outputSHA256'][volume_path.name]
        assert sha(surface_path) == entry['surfaceSHA256'] == receipt['outputSHA256'][surface_path.name]
    else:
        for stem, extension in (('restart_flow', 'csv'), ('surface_flow', 'vtk')):
            named = '%s_%06d.%s' % (stem, step, extension)
            assert sha(folder/named) == receipt['outputSHA256'][named] == sha(folder/(stem+'.'+extension))
        volume_path, surface_path = folder/'restart_flow.csv', folder/'surface_flow.vtk'
    spatial = receipt['runtimeFlavor'] == 'rah-ideal-air-axisymmetric-v1'
    data, values, nd = load_state(volume_path, turbulent, captured, spatial)
    axes = ['x', 'y', 'z'][:nd]
    rho = values['Density']
    velocity = np.column_stack([values['Momentum_'+axis]/rho for axis in axes])
    flow_velocity = np.column_stack([values['Velocity_'+axis] for axis in axes])
    kinetic = values['Turb_Kin_Energy'] if turbulent else np.zeros(len(data))
    if captured:
        previous_path = volume_path
        old_k = values['Flow_Primitive_K']
    else:
        previous_path = folder/('restart_flow_%06d.csv' % (step-1))
        assert sha(previous_path) == receipt['outputSHA256'][previous_path.name]
        previous, previous_values, previous_nd = load_state(previous_path, turbulent)
        assert previous_nd == nd and np.array_equal(previous[:, :nd+1], data[:, :nd+1])
        old_k = previous_values['Turb_Kin_Energy'] if turbulent else kinetic
    energy_without_k = values['Energy']/rho-.5*np.sum(velocity*velocity, axis=1)
    out = Path(output).resolve() if output else folder/'ideal_air_audit'
    out.relative_to(folder)
    out.mkdir(exist_ok=False)
    if checked_inventory is not None:
        save(out/'rechecked_inventory.json', dict(originalReceiptSHA256=sha(folder/'engine_process.json'),
            originalSupervisorError=receipt.get('supervisorError'), checkpoints=checked_inventory))
    # 独立程序使用固定CoolProp的umass_idealgas/cp0mass，核对完整保存阶段。
    helper = ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    upstream = read(ROOT/'outputs/v153_validation/native03/runtime.json')
    assert sha(helper) == upstream['sha256']
    write_rows(out/'oracle_input.csv', [('rhoT', float(r), float(t)) for r, t in zip(rho, values['Temperature'])])
    write_rows(out/'current_input.csv', [('rhoe', float(r), float(e)) for r, e in zip(rho, energy_without_k-kinetic)])
    inputs = [Path(__file__).resolve(), ROOT/'scripts/su2IdealAirInterface.py', ROOT/'scripts/su2SSTInterface.py',
              ROOT/'scripts/auditSU2SSTState.py', ROOT/'scripts/validateIdealAir.py', helper, previous_path,
              volume_path, surface_path, folder/'engine_process.json',
              out/'oracle_input.csv', out/'current_input.csv']
    lock = {str(path): sha(path) for path in inputs}
    save(out/'input_locks.json', lock)
    (out/'executed_audit.py').write_bytes(Path(__file__).read_bytes())
    for mode, option in (('oracle', '--oracle'), ('current', '--core')):
        if reuse_properties:
            old = Path(reuse_properties).resolve()
            old.relative_to(folder)
            old_locks = read(old/'input_locks.json')
            for path in (helper, volume_path, previous_path):
                assert sha(path) == old_locks[str(path)]
            assert sha(out/(mode+'_input.csv')) == sha(old/(mode+'_input.csv'))
            process = read(old/(mode+'_process.json'))
            assert process['status'] == 'terminal' and process['exitCode'] == 0
            (out/(mode+'.csv')).write_bytes((old/(mode+'.csv')).read_bytes())
            lock[str(old/(mode+'.csv'))] = sha(old/(mode+'.csv'))
        else:
            run(out, mode, [str(helper), option, str(out/(mode+'_input.csv')), str(out/(mode+'.csv'))], 30)
    oracle = read_rows(out/'oracle.csv')
    current = read_rows(out/'current.csv')
    assert np.all(oracle['accepted'] == 1) and np.all(current['accepted'] == 1)
    errors = dict(phaseEnergy=relative(energy_without_k-old_k, oracle['e']),
                  EOS=relative(values['Pressure'], oracle['p']), velocity=relative(velocity, flow_velocity))
    assert errors['phaseEnergy'] < 2e-10 and errors['EOS'] < 2e-11 and errors['velocity'] < 2e-11, errors
    g = wall_geometry(folder/'mesh.su2')
    points = np.column_stack([data[:, 1:nd+1], np.zeros((len(data), 3-nd))])
    xyz, surface = vtk_tokens(surface_path)
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
    errors['yPlus'] = relative(values['Y_Plus'][ids], yplus) if turbulent else None
    if turbulent:
        assert errors['yPlus'] < 2e-10
    else:
        # 原生层流路径未构造湍流壁距，Y_Plus保存为0；另列独立剪切尺度。
        assert np.all(values['Y_Plus'][ids] == 0)
    power = float(g['nodeArea'][ix]@q)
    errors['power'] = relative(np.array([history[step, keys.index('HF')]]), np.array([power]))
    assert errors['power'] < 2e-6
    residual_names = ['rms[Rho]', 'rms[RhoE]']+(['rms[k]', 'rms[w]'] if turbulent else [])
    residuals = [float(history[step, keys.index(name)]) for name in residual_names]
    stable = {path: sha(Path(path)) == expected for path, expected in lock.items()}
    assert all(stable.values())
    report = dict(interfacePassed=True, thermodynamicModel=contract, dimension=nd, nodes=len(data),
        wallNodes=len(ids), checkpointLocalIteration=step, successfulUpdates=step+1,
        seriesSuccessfulUpdates=(read(folder/'case_manifest.json')['initialization']['priorSuccessfulUpdates']+step+1) if captured else step+1,
        nativeNormalTerminal=normal,
        checkpointInventoryRechecked=recheck_checkpoints,
        elapsedSeconds=receipt['elapsedSeconds'], errors=errors, inputSHA256=lock, inputsStable=stable,
        nativeHeatPower=power, nativeHeatPowerUnit='W/m' if nd == 2 else 'W',
        wallHeatRangeWm2=[float(min(q)), float(max(q))], givenWallTemperatureK=wall_temperature,
        wallFluidTemperatureRangeK=[float(min(temperature[ids])), float(max(temperature[ids]))],
        temperatureRangeK=[float(min(temperature)), float(max(temperature))],
        currentKTemperatureRangeK=[float(min(current['T'])), float(max(current['T']))],
        currentKvsFlowPhaseTemperatureDifferenceK=float(max(abs(current['T']-temperature))),
        flowPrimitivePhase='captured SetPrimVar k input' if captured else 'flow update with previous saved k',
        turbulencePhase='subsequent turbulence update',
        previousKCheckpoint=None if captured else step-1, wallEddyToMolecularMaximum=eddy_fraction,
        finalResiduals=residuals, originalResidualThresholdPassed=step >= 10 and all(v < -8 for v in residuals),
        yPlusP99=float(np.percentile(yplus, 99)), nativeYPlusAvailable=bool(turbulent),
        reusedPropertiesFrom=str(reuse_properties) if reuse_properties else None,
        physicsQualified=False, newQualifiedArea=0)
    table = np.column_stack([ids, xyz, values['Pressure'][ids], q, temperature[ids],
                             np.full(len(ids), wall_temperature), yplus])
    np.savetxt(str(out/'wall.csv'), table, delimiter=',', fmt='%.17g', comments='',
        header='nodeID,x_m,y_m,z_m,p_Pa,q_W_m2,fluid_T_K,given_Tw_K,yPlus')
    report['wallCSV_SHA256'] = sha(out/'wall.csv')
    save(out/'report.json', report)
    if normal and checkpoint is None:
        publish_review(folder, out)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', required=True)
    parser.add_argument('--output')
    parser.add_argument('--reuse-properties')
    parser.add_argument('--checkpoint', type=int)
    parser.add_argument('--recheck-checkpoints', action='store_true')
    args = parser.parse_args()
    result = audit(args.folder, args.output, args.reuse_properties, args.checkpoint, args.recheck_checkpoints)
    print(json.dumps({k: v for k, v in result.items() if k not in ('inputSHA256', 'inputsStable')}, indent=2))
