"""模块用途。按事前窗口复核唯一完整域发展，不用退出码或云图替代物理资格。"""
import csv
import json
from pathlib import Path
import time
import numpy as np
from scipy.io import loadmat
from auditSU2IdealAirState import audit, load_state
from idealAirProperties import IdealAir
from prepareSU2WallDistanceState import boundary_faces
from reviewSU2SSTLocalDiffusion import neighbors_of
from reviewSU2AdmissibleRun import history
from su2SSTInterface import wall_geometry
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save, run, write_rows, read_rows

BASE = ROOT/'outputs/v166_validation'
CASE = BASE/'real02'
RESIDUALS = ('rms[Rho]', 'rms[RhoU]', 'rms[RhoV]', 'rms[RhoW]', 'rms[RhoE]', 'rms[k]', 'rms[w]')


def main():
    """模块一。完整三帧、严格热审核及输入身份；审核失败不放宽门槛。"""
    out = BASE/'review01'; out.mkdir(parents=True, exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    started = time.monotonic()
    receipt, manifest = read(CASE/'engine_process.json'), read(CASE/'case_manifest.json')
    assert receipt['status'] == 'terminal' and all(receipt['inputsStable'].values())
    normal = receipt['exitCode'] == 0 and receipt['checkpoints']['normalScheduleComplete']
    entries = receipt['checkpoints']['complete']
    steps = [x['localIteration'] for x in entries]
    assert steps == ([0, 100, 199] if normal else [0, 100])
    thermal = {}
    for step in steps:
        thermal[str(step)] = audit(CASE, CASE/('audit%03d'%step), checkpoint=None if normal and step == 199 else step)
    volume_path = ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'
    original_path = ROOT/'outputs/v56_validation/state_validity/state_validity_actual_rocket.mat'
    vol = loadmat(str(volume_path), squeeze_me=True)
    original = loadmat(str(original_path), squeeze_me=True, struct_as_record=False)['result']
    V, F, P = original.V, original.F.astype(int)-1, vol['vertices']
    assert np.array_equal(V, P[:len(V)]) and np.array_equal(F, vol['sourceFaces'].astype(int)-1)
    area = np.linalg.norm(np.cross(V[F[:, 1]]-V[F[:, 0]], V[F[:, 2]]-V[F[:, 0]]), axis=1)/2
    center = V[F].mean(axis=1)
    base = original.localGeometry.isAftTipAdjacent.astype(bool)
    nose = original.heating.noseFace.astype(bool) & ~base
    front = (center[:, 0] >= -.6) & ~base & ~nose
    groups = dict(nose=nose, forward_side=front, rear_side=~base & ~nose & ~front, base=base)
    assert np.all(sum(mask.astype(int) for mask in groups.values()) == 1)
    assert len(F) == 19752 and len(P) == 335726

    # 模块二。独立重建物理外边界梯度/通量，保留与原生离散边界通量的区别。
    markers = boundary_faces(CASE/'mesh.su2'); far, wall = markers['farfield'], markers['wall']
    cells = [vol['wedges'].astype(int)-1, vol['tetrahedra'].astype(int)-1]
    assert sum(len(g) for g in cells) == manifest['cells'] == 705490
    operators = {}
    for i, adjacent in neighbors_of(cells, np.unique(far), len(P)).items():
        js = np.array(sorted(adjacent), int); delta = P[js]-P[i]
        weighted = delta/np.sum(delta**2, axis=1)[:, None]
        operators[i] = (js, np.linalg.solve(weighted.T@delta, weighted.T))
    bq = P[np.concatenate((far, wall))]
    bav = -np.cross(bq[:, 1]-bq[:, 0], bq[:, 2]-bq[:, 0])/2
    volume_error = abs(np.einsum('ij,ij->', bav, bq.mean(axis=1))/(3*manifest['fluidVolume'])-1)
    assert volume_error < 2e-10
    fq = P[far]; far_av = -np.cross(fq[:, 1]-fq[:, 0], fq[:, 2]-fq[:, 0])/2
    geometry = wall_geometry(CASE/'mesh.su2'); wallids = geometry['wallIDs']
    core = IdealAir()
    paths = [Path(__file__).resolve(), volume_path, original_path, CASE/'engine_process.json',
             CASE/'case_manifest.json', CASE/'history.csv', CASE/'mesh.su2', ROOT/'docs/validation/v166/PLAN.md']
    locks = {str(p): sha(p) for p in paths}; locks.update(core.locks)
    frames, surfaces = [], []
    for entry in entries:
        path = CASE/entry['volume']; assert sha(path) == entry['volumeSHA256']; locks[str(path)] = sha(path)
        _, f, nd = load_state(path, True, captured=True)
        assert nd == 3
        rho, temperature, pressure = f['Density'], f['Temperature'], f['Pressure']
        velocity = np.column_stack([f['Momentum_'+a] for a in 'xyz'])/rho[:, None]
        thermo = core.state('PT', pressure, temperature)
        surface = {a: f[b][F].mean(axis=1) for a,b in
                   [('p','Pressure'), ('q','Heat_Flux'), ('T','Temperature'), ('yplus','Y_Plus')]}
        surfaces.append(surface)
        regions = {}
        for name, mask in groups.items():
            regions[name] = dict(faces=int(sum(mask)), areaM2=float(sum(area[mask])),
                pressureRangePa=[float(min(surface['p'][mask])), float(max(surface['p'][mask]))],
                heatRangeWm2=[float(min(surface['q'][mask])), float(max(surface['q'][mask]))],
                meanWallFluidTemperatureRangeK=[float(min(surface['T'][mask])), float(max(surface['T'][mask]))],
                physicalHeatPowerW=float(area[mask]@surface['q'][mask]), physicsQualified=False)
        values = np.column_stack((velocity, temperature)); grad = np.zeros((len(P), 4, 3))
        for i, (adjacent, op) in operators.items(): grad[i] = (op@(values[adjacent]-values[i])).T
        mu = f['Laminar_Viscosity']+f['Eddy_Viscosity']
        stress = mu[:, None, None]*(grad[:, :3]+grad[:, :3].transpose(0, 2, 1))
        stress -= (2*mu*np.trace(grad[:, :3], axis1=1, axis2=2)/3)[:, None, None]*np.eye(3)
        conductivity = thermo['cp']*(f['Laminar_Viscosity']/.72+f['Eddy_Viscosity']/.9)
        mass_flux = rho[:, None]*velocity
        energy_flux = (f['Energy']+pressure)[:, None]*velocity-np.einsum('nij,nj->ni', stress, velocity)-conductivity[:, None]*grad[:, 3]
        mass_faces = np.einsum('ij,ij->i', far_av, mass_flux[far].mean(axis=1))
        energy_faces = np.einsum('ij,ij->i', far_av, energy_flux[far].mean(axis=1))
        power = float(area@surface['q'])
        mass_in, energy_in = float(sum(np.maximum(-mass_faces, 0))), float(sum(np.maximum(-energy_faces, 0)))
        assert mass_in > 0 and energy_in > 0
        balance = dict(massNetKgS=float(sum(mass_faces)), massIncomingKgS=mass_in,
            massRelativeImbalance=float(abs(sum(mass_faces))/mass_in),
            energyNetW=float(sum(energy_faces)+power), energyIncomingW=energy_in,
            energyRelativeImbalance=float(abs(sum(energy_faces)+power)/energy_in),
            method='v160 physical flux diagnostic: cell-edge WLS, effective Newtonian stress and Pr conduction, nodal flux triangle mean',
            nativeDiscreteBoundaryFluxAvailable=False, nativeConservationProven=False,
            normalizationReference='OTH total energy including captured primitive SST k')
        error = temperature[F]-300
        # 三角形上线性温差的平方精确面积平均，避免仅取面均值掩盖点偏差。
        square_mean = (np.sum(error, axis=1)**2+np.sum(error**2, axis=1))/12
        frame = dict(localIteration=entry['localIteration'], seriesSuccessfulUpdates=entry['seriesSuccessfulUpdates'],
            givenWallTemperatureK=300, wallFluidTemperatureRangeK=[float(min(temperature[wallids])), float(max(temperature[wallids]))],
            maximumWallTemperatureErrorK=float(max(abs(temperature[wallids]-300))),
            areaRMSWallTemperatureErrorK=float(np.sqrt(area@square_mean/sum(area))),
            domainTemperatureRangeK=[float(min(temperature)), float(max(temperature))],
            yPlusP99=float(np.quantile(f['Y_Plus'][wallids], .99)), physicalSurfaceHeatPowerW=power,
            nativeDualHeatPowerW=float(geometry['nodeArea']@f['Heat_Flux'][wallids]),
            boundaryBalanceDiagnostic=balance, regions=regions, physicsQualified=False, newQualifiedArea=0)
        frames.append(frame); save(out/('frame_%06d.json'%entry['localIteration']), frame)

    # 模块三。冷重启跳变与本段发展分开；所有阈值取自启动前计划。
    oldpath = ROOT/'outputs/v164_validation/real01/restart_flow_000049.csv'
    locks[str(oldpath)] = sha(oldpath)
    _, old, _ = load_state(oldpath, True, captured=True)
    jump = {a: float(area@abs(surfaces[0][a]-old[b][F].mean(axis=1))/(area@abs(old[b][F].mean(axis=1))))
            for a,b in [('p','Pressure'),('q','Heat_Flux')]}
    def change(first, last, field, updates):
        return float(area@abs(surfaces[last][field]-surfaces[first][field])/(area@abs(surfaces[first][field]))*100/updates)
    changes = {field: dict(first100=change(0, 1, field, 100),
        last99Per100=change(1, 2, field, 99) if normal else None) for field in ('p', 'q')}
    keys, h = history(CASE/'history.csv'); assert len(h) == (200 if normal else 180) and np.isfinite(h).all()
    mean = {window: {name: float(np.mean(rows[:, keys.index(name)])) for name in RESIDUALS}
            for window,rows in [('first20',h[:20]),('previous20',h[-40:-20]),('last20',h[-20:])]}
    delta = {name: mean['last20'][name]-mean['first20'][name] for name in RESIDUALS}
    tail = {name: mean['last20'][name]-mean['previous20'][name] for name in RESIDUALS}
    power = h[:, keys.index('HF')]
    power_range = dict(first20=float(np.ptp(power[:20])/np.mean(abs(power[:20]))),
                       last20=float(np.ptp(power[-20:])/np.mean(abs(power[-20:]))))
    tracepath = CASE/'admissibility.csv'; events = []
    if tracepath.exists():
        assert sha(tracepath) == receipt['outputSHA256'][tracepath.name]; locks[str(tracepath)] = sha(tracepath)
        with tracepath.open() as stream: events = list(csv.DictReader(stream))
        assert len(events) <= 4096
        assert [int(e['event']) for e in events] == list(range(1,len(events)+1))
    # 模块四。失败时只复核已保存状态和实际事件；不补造199帧或异常瞬间的全场。
    def event_array(key): return np.array([float(e[key]) for e in events])
    def event_state(prefix): return np.column_stack([event_array(prefix+'_'+n) for n in ('rho','rhou','rhov','rhow','rhoE')])
    trial, old_state, k = event_state('trial'), event_state('old'), event_array('k')
    def e_of(u): return u[:,4]/u[:,0]-.5*np.sum((u[:,1:4]/u[:,0,None])**2,axis=1)-k
    old_e, trial_e = e_of(old_state), e_of(trial)
    replay = max(float(max(abs(old_e-event_array('old_e')))), float(max(abs(trial_e-event_array('trial_e')))))
    assert replay < 1e-8
    recovered = event_array('recovered').astype(bool)
    bounds = core.state('PT', 101325, [60, 2000])['e']
    assert np.all((old_e[recovered] >= bounds[0]) & (old_e[recovered] <= bounds[1]))
    recovered_thermo = core.state('rhoe', old_state[recovered,0], old_e[recovered])
    helper = ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    assert sha(helper) == read(ROOT/'outputs/v153_validation/native03/runtime.json')['sha256']
    oracle_rows = [('rhoT', float(r), float(t)) for r,t in zip(old_state[recovered,0], recovered_thermo['T'])]
    oracle_rows.append(('rhoT', float(old_state[-1,0]), 60.))
    write_rows(out/'event_oracle_input.csv', oracle_rows)
    run(out, 'event_oracle', [str(helper), '--oracle', str(out/'event_oracle_input.csv'), str(out/'event_oracle.csv')], 30)
    oracle = read_rows(out/'event_oracle.csv'); assert np.all(oracle['accepted'] == 1)
    oracle_error = float(max(abs(oracle['e'][:-1]-old_e[recovered])/np.maximum(1,abs(old_e[recovered]))))
    assert oracle_error < 2e-11 and abs(oracle['e'][-1]-bounds[0]) < 1e-7
    failure = None
    if not normal:
        assert len(events) == 338 and sum(~recovered) == 1 and not recovered[-1]
        point = events[-1]['local_point']
        before = max(i for i,event in enumerate(events[:-1]) if event['local_point'] == point)
        assert np.array_equal(old_state[-1], trial[-1]) and np.array_equal(old_state[before], old_state[-1])
        assert old_e[before] >= bounds[0] > old_e[-1]
        coupling_error = abs((old_e[-1]-old_e[before])+(k[-1]-k[before]))
        assert coupling_error < 1e-8
        failure = dict(recoveredEvent=int(events[before]['event']), failedEvent=int(events[-1]['event']),
            localRCMPoint=int(point), conservedComponentsExact=True, globalLocationKnown=False,
            oldK=float(k[before]), newK=float(k[-1]), kIncreaseJkg=float(k[-1]-k[before]),
            minimumThermalEnergyJkg=float(bounds[0]), priorThermalEnergyJkg=float(old_e[before]),
            failedThermalEnergyJkg=float(old_e[-1]), deficitJkg=float(bounds[0]-old_e[-1]),
            pairedEnergyChangeErrorJkg=float(coupling_error),
            actualSSTIncrementAndRHSAvailable=False,
            scope='same five conserved values and recorded current k explain this EOS rejection; not the cause of whole-domain nonconvergence')
    locks[str(helper)] = sha(helper)
    first, last = frames[0], frames[-1]
    gases = ('rms[Rho]', 'rms[RhoE]'); turbulence = ('rms[k]', 'rms[w]')
    gates = dict(normal200Updates=bool(normal and len(h)==200), strictStateAndSource=True, originalGeometry=True,
        densityEnergyDevelopment=all(delta[n] <= -.30 for n in gases),
        turbulenceDevelopment=all(delta[n] <= .10 for n in turbulence),
        flowLateWindow=all(tail[n] <= .05 for n in gases) and any(tail[n] <= -.05 for n in gases),
        turbulenceLateWindow=all(tail[n] <= .10 for n in turbulence),
        pressureStability=changes['p']['last99Per100'] <= .9*changes['p']['first100'] if normal else None,
        heatStability=changes['q']['last99Per100'] <= .9*changes['q']['first100'] if normal else None,
        powerStability=power_range['last20'] <= .9*power_range['first20'],
        maximumWallTemperature=last['maximumWallTemperatureErrorK'] <= first['maximumWallTemperatureErrorK']+1,
        rmsWallTemperature=last['areaRMSWallTemperatureErrorK'] <= .95*first['areaRMSWallTemperatureErrorK'],
        yPlusNoRegression=last['yPlusP99'] <= 1.05*first['yPlusP99'])
    for field in ('mass', 'energy'):
        key = field+'RelativeImbalance'
        gates[field+'BalanceNoRegression'] = last['boundaryBalanceDiagnostic'][key] <= 1.1*first['boundaryBalanceDiagnostic'][key]+1e-5
    if not normal:
        # 事前收益窗要求200次完整场，提前失败时这些门槛待验，禁止换成100或179次授予通过。
        gates = {key: value if key in ('normal200Updates','strictStateAndSource','originalGeometry') else None for key,value in gates.items()}
    assert all(sha(Path(p)) == digest for p,digest in locks.items())
    report = dict(passedChecks=True, developmentBenefitPassed=all(gates.values()), gates=gates,
        thermalErrors={step:r['errors'] for step,r in thermal.items()}, frames=frames,
        coldRestartSurfaceJump=jump, surfaceChangePer100Updates=changes, nativePowerRelativeRanges=power_range,
        meanLogResiduals=mean, lastMinusFirst20Log=delta, lastMinusPrevious20Log=tail,
        finalResiduals={name: float(h[-1,keys.index(name)]) for name in RESIDUALS},
        originalResidualThresholdPassed=thermal[str(steps[-1])]['originalResidualThresholdPassed'],
        originalYPlusThresholdPassed=last['yPlusP99'] <= 1, regionDefinitions='v56 nose and base masks; side split at X=-0.6 m; statistics only',
        admissibilityEvents=len(events), recoveredEvents=int(sum(recovered)), failedOldStates=int(sum(~recovered)),
        energyReplayErrorJkg=replay, eventOracleEnergyError=oracle_error, failurePair=failure,
        eventGlobalLocationsKnown=False, eventVolumeCorrectionsAvailable=False,
        normalTerminal=bool(normal), lastSavedLocalIteration=steps[-1], lastFiniteHistoryLocalIteration=len(h)-1,
        historyWindowIndices=dict(first20=[0,19], previous20=[len(h)-40,len(h)-21], last20=[len(h)-20,len(h)-1]),
        requestedUpdates=200, finiteHistoryUpdates=len(h), priorUpdates=200, seriesFiniteHistoryUpdates=200+len(h),
        failureInstantFullStateAvailable=False, plannedFinalWindowAvailable=bool(normal),
        originalFaces=len(F), originalAreaM2=float(sum(area)), boundaryVolumeRelativeError=volume_error,
        inputSHA256=locks, inputsStable=True, elapsedSeconds=time.monotonic()-started,
        automaticAdditionalTasks=0, engineCompilation=0, physicsQualified=False, newQualifiedArea=0)
    save(out/'report.json', report)
    print(json.dumps({k:report[k] for k in ('developmentBenefitPassed','gates','lastMinusFirst20Log',
        'lastMinusPrevious20Log','surfaceChangePer100Updates','nativePowerRelativeRanges','finalResiduals','admissibilityEvents','elapsedSeconds')}, indent=2))


if __name__ == '__main__': main()
