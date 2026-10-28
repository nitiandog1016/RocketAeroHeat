# -*- coding: utf-8 -*-
"""保留真实三维新方程每个方向的温域限制位置和状态，不更改候选或门槛。"""
from pathlib import Path
import hashlib, json, os, sys
import numpy as np
from scipy.io import loadmat
root = Path(__file__).resolve().parents[1]
val = root / 'outputs/v82_validation'
case = sys.argv[1] if len(sys.argv) > 1 else 'exact_krylov01'
assert case.replace('_', '').isalnum()
target = val / (case + '_candidate_diagnostic.json')
assert not target.exists()
source = val / case / 'steady_krylov_evidence.mat'
domain_file = root / 'outputs/v80_validation/domain02/shell_domain_evidence.mat'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
bindings = {str(p.relative_to(root)): sha(p) for p in [source, domain_file]}
s = loadmat(source, simplify_cells=True)
d = loadmat(domain_file, simplify_cells=True, variable_names=['domain'])['domain']
r, gas = s['result'], s['gas']
cv = gas['R'] / (gas['gamma'] - 1)
trace = r['trace']; trace = [trace] if isinstance(trace, dict) else list(trace)
rows = []
for item in trace:
    initial = item['initialQ']
    direction = item['linearUpdateU'] * r['residualScale']
    for alpha in sorted(set([1., .5, .25, float(item['alpha'])]), reverse=True):
        q = initial + alpha * direction
        rho = q[:, 0]
        u = q[:, 1:4] / rho[:, None]
        T = (q[:, 4] - .5 * rho * np.sum(u ** 2, axis=1)) / (rho * cv)
        lo, hi = int(np.argmin(T)), int(np.argmax(T))
        rows.append(dict(iteration=int(item['iteration']), alpha=alpha,
                         isActualAcceptedAlpha=alpha == item['alpha'], minDensity=float(rho.min()),
                         temperatureRangeK=[float(T.min()), float(T.max())],
                         below80KCount=int(np.count_nonzero(T < gas['temperatureRange'][0])),
                         above2000KCount=int(np.count_nonzero(T > gas['temperatureRange'][1])),
                         minimumTemperatureCV=lo + 1, maximumTemperatureCV=hi + 1,
                         minimumTemperatureCentroidM=d['cells']['centroid'][lo].tolist(),
                         maximumTemperatureCentroidM=d['cells']['centroid'][hi].tolist(),
                         minimumTemperatureLayer=lo // 19752 + 1, maximumTemperatureLayer=hi // 19752 + 1))
z = abs(r['rate'] / r['residualScale'])
index = np.unravel_index(np.argmax(z), z.shape)
dominant = dict(CV=int(index[0]) + 1, equation=int(index[1]) + 1,
                scaledResidual=float(z[index]), centroidM=d['cells']['centroid'][index[0]].tolist(),
                layer=int(index[0]) // 19752 + 1)
assert all(sha(root / p) == h for p, h in bindings.items())
report = dict(status='observed-actual-full-NS-candidate-temperature-range-without-clipping',
              case=case, processID=os.getpid(), directions=rows, finalDominantResidual=dominant,
              modelTemperatureRangeUnchanged=True, physicalSolution=False, newProductionSolvedAreaM2=0,
              inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
print(json.dumps(report))
