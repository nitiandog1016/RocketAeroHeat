# -*- coding: utf-8 -*-
"""由真实原面通量分解收缩冷区的对流、Fourier和黏性内能收支；不改变方程。"""
from pathlib import Path
import hashlib, json, os, sys
import numpy as np
from scipy.io import loadmat
root = Path(__file__).resolve().parents[1]
stage = sys.argv[2] if len(sys.argv)>2 else 'v82_validation'
assert stage in ('v82_validation','v83_validation')
val = root / 'outputs' / stage
case = sys.argv[1] if len(sys.argv) > 1 else 'exact_krylov01'
assert case.replace('_', '').isalnum()
target = val / (case + '_energy_diagnostic.json'); assert not target.exists()
source = val / case / ('wall_velocity_krylov_evidence.mat' if stage=='v83_validation' else 'steady_krylov_evidence.mat')
domain_file = root / 'outputs/v80_validation/domain02/shell_domain_evidence.mat'
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()
bindings = {str(p.relative_to(root)): sha(p) for p in [source, domain_file]}
s = loadmat(source, simplify_cells=True)
d = loadmat(domain_file, simplify_cells=True, variable_names=['domain'])['domain']
r, gas = s['result'], s['gas']; info = r['info']; P = info['primitive']
f, w = d['openFaces'], d['walls']; volume = d['cells']['volume']
m, p, owner = f['cellMinus'].astype(int) - 1, f['cellPlus'].astype(int) - 1, w['ownerCell'].astype(int) - 1
cv = gas['R'] / (gas['gamma'] - 1)
rows = []
selected=[14280,16339,int(np.argmin(P[:,4])),int(np.argmax(abs(r['rate'][:,4]/r['residualScale'][4])))]
if stage=='v83_validation':selected += [16396,16398,16400,36152]
for cell in sorted(set(selected)):
    ids = np.flatnonzero((m == cell) | (p == cell)); wi = np.flatnonzero(owner == cell)
    assert np.all((m[ids] >= 0) & (p[ids] >= 0)), '当前诊断选取远离开放边界的实际近壁CV'
    sign = np.where(m[ids] == cell, 1., -1.)
    area = f['area'][ids]
    convection = -np.sum((sign * area)[:, None] * info['convectiveFaceFlux'][ids], axis=0)
    wall_conv = np.zeros((len(wi), 5)); wall_conv[:, 1:4] = info['wallPressure'][wi, None] * w['normal'][wi]
    convection -= np.sum(w['area'][wi, None] * wall_conv, axis=0)
    diffusion = np.sum((sign * area)[:, None] * info['viscousFaceFlux'][ids], axis=0)
    wall_visc = wall_conv - info['wallFlux'][wi]
    diffusion += np.sum(w['area'][wi, None] * wall_visc, axis=0)
    assert np.max(abs((convection + diffusion) / volume[cell] - r['rate'][cell])) / max(1, np.max(abs(r['rate'][cell]))) < 3e-13
    rho, u, T = P[cell, 0], P[cell, 1:4], P[cell, 4]
    def slope(integrated_rate):
        v = integrated_rate / volume[cell]
        return float((v[4] - np.dot(u, v[1:4]) + (.5 * np.dot(u, u) - cv * T) * v[0]) / (rho * cv))
    left, right = m[ids], p[ids]
    face_velocity = .5 * (P[left, 1:4] + P[right, 1:4])
    traction = info['viscousFaceFlux'][ids, 1:4]
    work_flux = np.sum(face_velocity * traction, axis=1)
    fourier = info['viscousFaceFlux'][ids, 4] - work_flux
    face_heat = sign * area * fourier / (volume[cell] * rho * cv)
    face_mechanical = sign * area * np.sum((face_velocity - u) * traction, axis=1) / (volume[cell] * rho * cv)
    wall_heat = float(np.sum(w['area'][wi] * wall_visc[:, 4]) / (volume[cell] * rho * cv))
    wall_mechanical = float(-np.sum(w['area'][wi] * (wall_visc[:, 1:4] @ u)) / (volume[cell] * rho * cv))
    conduction = float(face_heat.sum()) + wall_heat
    mechanical = float(face_mechanical.sum()) + wall_mechanical
    assert abs(conduction + mechanical - slope(diffusion)) / max(1, abs(slope(diffusion))) < 3e-13
    contributors = []
    for j in np.argsort(face_mechanical)[:5]:
        contributors.append(dict(faceID=int(ids[j]) + 1, neighbourCV=int(right[j] if left[j] == cell else left[j]) + 1,
                                 centroidM=f['centroid'][ids[j]].tolist(), mechanicalTemperatureRateKs=float(face_mechanical[j]),
                                 fourierTemperatureRateKs=float(face_heat[j])))
    rows.append(dict(CV=cell + 1, layer=cell // 19752 + 1, centroidM=d['cells']['centroid'][cell].tolist(), temperatureK=float(T),
                     convectiveTemperatureRateKs=slope(convection), diffusiveTemperatureRateKs=slope(diffusion),
                     fourierTemperatureRateKs=conduction, mechanicalInternalTemperatureRateKs=mechanical,
                     wallFourierTemperatureRateKs=wall_heat, wallMechanicalTemperatureRateKs=wall_mechanical,
                     totalTemperatureRateKs=slope(convection + diffusion), strongestNegativeMechanicalContributors=contributors))
assert all(sha(root / p) == h for p, h in bindings.items())
report = dict(status='observed-real-exact-NS-contraction-cooling-energy-decomposition', case=case, processID=os.getpid(),
              rows=rows, newEquationSolved=False, physicalSolution=False, newProductionSolvedAreaM2=0,
              inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf8'); print(json.dumps(report))
