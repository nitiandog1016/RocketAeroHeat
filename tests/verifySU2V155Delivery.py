"""只读复核v155初场和唯一小型终态，明确保留缺失前帧k的审核失败。"""
from pathlib import Path
import csv
import json
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from validateIdealAir import run, save, sha, write_rows, read_rows, relative
from su2ThermalRuntime import read
from su2IdealAirInterface import settings
from su2SSTInterface import wall_geometry
from auditSU2IdealAirState import load_state
from auditSU2SSTState import vtk_tokens


def main():
    output = ROOT/'outputs/v155_validation'
    review = output/'review01'
    review.mkdir(exist_ok=False)
    folder = output/'cases01/small'
    receipt = read(folder/'engine_process.json')
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0 and not receipt['timedOut'] and not receipt['cancelled']
    assert all(receipt['inputsStable'].values())
    for name in ('restart_flow.csv', 'surface_flow.vtk', 'history.csv'):
        assert sha(folder/name) == receipt['outputSHA256'][name]
    assert not (folder/'restart_flow_000000.csv').exists()
    assert not (folder/'ideal_air_review.json').exists()
    failure = read(output/'supervision/small_audit_process.json')
    assert failure['exitCode'] == 1
    assert 'restart_flow_000000.csv' in (output/'supervision/small_audit.log').read_text()
    data, values, dimension = load_state(folder/'restart_flow.csv', True)
    assert dimension == 3 and len(data) == 38
    helper = ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    assert sha(helper) == read(helper.parent/'runtime.json')['sha256']
    write_rows(review/'oracle_input.csv', [('rhoT', float(r), float(t))
        for r, t in zip(values['Density'], values['Temperature'])])
    run(review, 'oracle', [str(helper), '--oracle', str(review/'oracle_input.csv'), str(review/'oracle.csv')], 30)
    thermo = read_rows(review/'oracle.csv')
    assert np.all(thermo['accepted'] == 1)
    c = settings((folder/'case.cfg').read_text())
    geometry = wall_geometry(folder/'mesh.su2')
    ids = geometry['wallIDs']
    wall_temperature = float(re.fullmatch(r'\(\s*wall\s*,\s*([-+\deE.]+)\s*\)', c['MARKER_ISOTHERMAL'])[1])
    temperature = values['Temperature']
    mu = float(c['MU_REF'])*(temperature/float(c['MU_T_REF']))**1.5
    mu *= (float(c['MU_T_REF'])+float(c['SUTHERLAND_CONSTANT']))/(temperature+float(c['SUTHERLAND_CONSTANT']))
    q = mu[ids]*thermo['cp'][ids]/float(c['PRANDTL_LAM'])*(temperature[geometry['normalNeighbor']]-wall_temperature)/geometry['thermalDistance']
    surface_xyz, surface = vtk_tokens(folder/'surface_flow.vtk')
    distance, matched = cKDTree(data[:, 1:4]).query(surface_xyz)
    assert max(distance) < 1e-6 and set(matched) == set(ids)
    native_power = float(geometry['nodeArea'] @ q)
    with (folder/'history.csv').open() as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    keys = [key.strip() for key in rows[0]]
    history = np.asarray(rows[1:], float)
    assert np.array_equal(history[:, keys.index('Inner_Iter')], [0, 1])
    errors = dict(EOS=relative(values['Pressure'], thermo['p']),
        molecularViscosity=relative(values['Laminar_Viscosity'], mu),
        heatGradient=relative(values['Heat_Flux'][ids], q),
        surfaceHeat=relative(surface['Heat_Flux'][:, 0], values['Heat_Flux'][matched]),
        power=relative(np.array([native_power]), np.array([history[-1, keys.index('HF')]])))
    assert errors['EOS'] < 2e-11 and errors['heatGradient'] < 2e-10
    assert errors['molecularViscosity'] < 2e-11 and errors['surfaceHeat'] < 2e-6 and errors['power'] < 2e-6
    residuals = {name: float(history[-1, keys.index(name)]) for name in ('rms[Rho]', 'rms[RhoE]', 'rms[k]', 'rms[w]')}
    native_report = dict(successfulUpdates=2, localIteration=1, elapsedSeconds=receipt['elapsedSeconds'],
        exitCode=receipt['exitCode'], errors=errors, finalResiduals=residuals,
        originalResidualThresholdPassed=all(value < -8 for value in residuals.values()),
        givenWallTemperatureK=wall_temperature,
        wallFluidTemperatureRangeK=[float(min(temperature[ids])), float(max(temperature[ids]))],
        heatFluxRangeWm2=[float(min(q)), float(max(q))], heatPowerW=native_power,
        phaseEnergyIndependentCheck='pending: native steady output omits iteration zero',
        interfacePassed=False, physicsQualified=False, newQualifiedArea=0)
    initial = {name: read(output/'cases01'/name/'initial_independent_audit/report.json') for name in ('small', 'real')}
    assert all(item['passed'] for item in initial.values())
    matlab = read(output/'matlab/export01_process.json')
    assert matlab['status'] == 'terminal' and matlab['exitCode'] == 0 and all(matlab['sourceStableAfter'].values())
    old_locks = read(output/'source_lock.json')
    stable = {name: sha(Path(name)) == digest for name, digest in old_locks.items()}
    assert all(stable.values())
    for name in ('small', 'real'):
        record = read(output/'cases01'/name/'ideal_air_preparation.json')
        assert record['exitCode'] == 0 and all(record['inputsStable'].values())
        assert sha(output/'cases01'/name/'initial_flow.csv') == initial[name]['initialStateSHA256']
    initial_errors = [item['errors']['thermalEnergy'] for item in initial.values()]
    report = dict(initializationPassed=True, nativeFullInterfacePassed=False, initialStates=initial,
        maximumIndependentInitialEnergyError=max(initial_errors), nativeSmall=native_report,
        fullDomainSolverRuns=0, matlabSourceCount=len(matlab['sourceSHA256']),
        inputsStable=stable, sourceSHA256={str(Path(__file__).resolve()): sha(Path(__file__).resolve())},
        physicsQualified=False, newQualifiedArea=0)
    save(review/'report.json', report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('inputsStable', 'sourceSHA256')}))


if __name__ == '__main__':
    main()
