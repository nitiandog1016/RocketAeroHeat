"""真实整箭三帧的原面分区、阶段状态及独立物理边界积分。"""
import csv
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from auditSU2IdealAirState import load_state
from idealAirProperties import IdealAir
from prepareSU2WallDistanceState import boundary_faces
from reviewSU2SSTLocalDiffusion import neighbors_of
from su2SSTInterface import wall_geometry
from su2ThermalRuntime import read, sha

BASE = ROOT/'outputs/v160_validation'
CASE = BASE/'real01'
OUT = BASE/'field_review02'
OUT.mkdir(exist_ok=False)
started = time.monotonic()
receipt = read(CASE/'engine_process.json')
manifest = read(CASE/'case_manifest.json')
assert receipt['exitCode'] == 0 and all(receipt['inputsStable'].values())
volume_path = ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'
original_path = ROOT/'outputs/v56_validation/state_validity/state_validity_actual_rocket.mat'
vol = loadmat(str(volume_path), squeeze_me=True)
original = loadmat(str(original_path), squeeze_me=True, struct_as_record=False)['result']
V, F = original.V, original.F.astype(int)-1
P = vol['vertices']
assert np.array_equal(V, P[:len(V)]) and np.array_equal(F, vol['sourceFaces'].astype(int)-1)
av = np.cross(V[F[:, 1]]-V[F[:, 0]], V[F[:, 2]]-V[F[:, 0]])/2
area = np.linalg.norm(av, axis=1)
center = V[F].mean(axis=1)
base = original.localGeometry.isAftTipAdjacent.astype(bool)
nose = original.heating.noseFace.astype(bool) & ~base
front_side = (center[:, 0] >= -.6) & ~base & ~nose
rear_side = ~base & ~nose & ~front_side
groups = dict(original_nose_heating_region=nose, forward_side_to_x_minus_0_6=front_side,
              rear_side=rear_side, original_aft_fan=base)
assert np.all(sum(mask.astype(int) for mask in groups.values()) == 1)
markers = boundary_faces(CASE/'mesh.su2')
far = markers['farfield']
wall = markers['wall']
far_ids = np.unique(far)
groups_of_cells = [vol['wedges'].astype(int)-1, vol['tetrahedra'].astype(int)-1]
assert sum(len(g) for g in groups_of_cells) == manifest['cells']
neighbors = neighbors_of(groups_of_cells, far_ids, len(P))
operators = {}
for i, adjacent in neighbors.items():
    js = np.array(sorted(adjacent), int)
    delta = P[js]-P[i]
    weighted = delta/np.sum(delta**2, axis=1)[:, None]
    operators[i] = (js, np.linalg.solve(weighted.T@delta, weighted.T))
all_boundary = np.concatenate((far, wall))
bq = P[all_boundary]
bav = -np.cross(bq[:, 1]-bq[:, 0], bq[:, 2]-bq[:, 0])/2
volume_by_boundary = np.sum(np.einsum('ij,ij->i', bav, bq.mean(axis=1)))/3
assert abs(volume_by_boundary/manifest['fluidVolume']-1) < 2e-10
far_q = P[far]
far_av = -np.cross(far_q[:, 1]-far_q[:, 0], far_q[:, 2]-far_q[:, 0])/2
geometry = wall_geometry(CASE/'mesh.su2')
core = IdealAir()
inputs = [Path(__file__).resolve(), volume_path, original_path, CASE/'engine_process.json', CASE/'case_manifest.json']
locks = {str(p): sha(p) for p in inputs}
locks.update(core.locks)
frames = []
face_fields = []
for entry in receipt['checkpoints']['complete']:
    path = CASE/entry['volume']
    assert sha(path) == entry['volumeSHA256']
    locks[str(path)] = entry['volumeSHA256']
    raw, f, nd = load_state(path, True, True)
    assert nd == 3 and np.allclose(raw[:, 1:4], P, rtol=2e-14, atol=1e-14)
    rho, T, p = f['Density'], f['Temperature'], f['Pressure']
    U = np.column_stack([f['Momentum_'+axis] for axis in 'xyz'])/rho[:, None]
    thermo = core.state('PT', p, T)
    internal = f['Energy']/rho-np.sum(U**2, axis=1)/2-f['Flow_Primitive_K']
    eos = float(max(abs(rho/thermo['rho']-1)))
    energy_error = float(max(abs(internal-thermo['e'])/np.maximum(1, abs(thermo['e']))))
    ids = geometry['wallIDs']
    heat = f['Laminar_Viscosity'][ids]*thermo['cp'][ids]/.72
    heat *= (T[geometry['normalNeighbor']]-300)/geometry['thermalDistance']
    heat_error = float(max(abs(heat-f['Heat_Flux'][ids])/np.maximum(1, abs(heat))))
    assert eos < 2e-12 and energy_error < 2e-10 and heat_error < 2e-6
    fields = dict(p=p[F].mean(axis=1), q=f['Heat_Flux'][F].mean(axis=1),
                  temperature=T[F].mean(axis=1), yplus=f['Y_Plus'][F].mean(axis=1))
    face_fields.append(fields)
    table = {}
    for name, mask in groups.items():
        table[name] = dict(faces=int(sum(mask)), areaM2=float(sum(area[mask])),
            pressureRangePa=[float(min(fields['p'][mask])), float(max(fields['p'][mask]))],
            heatRangeWm2=[float(min(fields['q'][mask])), float(max(fields['q'][mask]))],
            physicalSurfaceHeatPowerW=float(area[mask]@fields['q'][mask]),
            maximumMeanWallFluidTemperatureErrorK=float(max(abs(fields['temperature'][mask]-300))))
    # 独立WLS边界梯度及牛顿流体物理通量，保留与原生离散通量的证据区别。
    values = np.column_stack((U, T))
    grad = np.zeros((len(P), 4, 3))
    for i, (adjacent, op) in operators.items():
        grad[i] = (op@(values[adjacent]-values[i])).T
    mu = f['Laminar_Viscosity']+f['Eddy_Viscosity']
    stress = mu[:, None, None]*(grad[:, :3]+grad[:, :3].transpose(0, 2, 1))
    stress -= (2*mu*np.trace(grad[:, :3], axis1=1, axis2=2)/3)[:, None, None]*np.eye(3)
    kappa = thermo['cp']*(f['Laminar_Viscosity']/.72+f['Eddy_Viscosity']/.9)
    mass_flux = rho[:, None]*U
    energy_flux = (f['Energy']+p)[:, None]*U-np.einsum('nij,nj->ni', stress, U)-kappa[:, None]*grad[:, 3]
    mass_faces = np.einsum('ij,ij->i', far_av, mass_flux[far].mean(axis=1))
    energy_faces = np.einsum('ij,ij->i', far_av, energy_flux[far].mean(axis=1))
    wall_power = float(area@fields['q'])
    balance = dict(massNetKgS=float(sum(mass_faces)),
        massRelativeImbalance=float(abs(sum(mass_faces))/sum(np.maximum(-mass_faces, 0))),
        energyNetW=float(sum(energy_faces)+wall_power),
        energyRelativeImbalance=float(abs(sum(energy_faces)+wall_power)/sum(np.maximum(-energy_faces, 0))),
        method='independent physical flux with cell-edge WLS; linear nodal flux triangle integration; face-net incoming normalization',
        nativeDiscreteBoundaryFluxAvailable=False, normalizationReference='OTH energy including captured primitive SST k')
    h0 = (f['Energy']+p)/rho
    hot = int(np.argmax(T[ids]))
    frame = dict(iteration=entry['localIteration'], updates=entry['successfulUpdates'],
        EOSRelativeError=eos, phaseEnergyRelativeError=energy_error, heatGradientRelativeError=heat_error,
        temperatureRangeK=[float(min(T)), float(max(T))],
        wallFluidTemperatureRangeK=[float(min(T[ids])), float(max(T[ids]))],
        hottestWallNode=dict(globalNodeID=int(ids[hot]), positionM=P[ids[hot]].tolist(), temperatureK=float(T[ids[hot]])),
        yPlus99=float(np.quantile(f['Y_Plus'][ids], .99)),
        physicalSurfaceHeatPowerW=wall_power, nativeDualHeatPowerW=float(geometry['nodeArea']@f['Heat_Flux'][ids]),
        regions=table, boundaryBalanceDiagnostic=balance, totalEnthalpyMaximumJkg=float(max(h0)),
        physicsQualified=False, newQualifiedArea=0)
    frames.append(frame)
    (OUT/('frame_%06d.json' % entry['localIteration'])).write_text(json.dumps(frame, indent=2)+'\n', encoding='utf8')
last, previous = face_fields[-1], face_fields[-2]
stability = {field+'AreaL1RelativeChange':float(area@abs(last[field]-previous[field])/(area@abs(last[field]))) for field in ('p', 'q')}
stability.update(firstIteration=frames[-2]['iteration'], lastIteration=frames[-1]['iteration'],
    physicalPowerRelativeChange=abs(frames[-1]['physicalSurfaceHeatPowerW']/frames[-2]['physicalSurfaceHeatPowerW']-1),
    minimum1000UpdateWindowAvailable=False)
assert all(sha(Path(path)) == digest for path, digest in locks.items())
report = dict(passedStateAndMappingChecks=True, frames=frames, shortWindow=stability,
    originalFaces=len(F), originalAreaM2=float(sum(area)),
    regionDefinition='original v56 nose-heating mask and aft-tip fan; remaining sidewall split at X=-0.6 m for statistics',
    boundaryVolumeRelativeError=abs(volume_by_boundary/manifest['fluidVolume']-1),
    sourceSHA256=locks, inputsStable=True, elapsedSeconds=time.monotonic()-started,
    physicsQualified=False, newQualifiedArea=0)
(OUT/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
print(json.dumps(dict(elapsedSeconds=report['elapsedSeconds'], frames=len(frames), shortWindow=stability, final=frames[-1]), indent=2))
