"""只读两个原生检查点，量化真实域启动变化及总温超调，不再运行求解器。"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/v118_validation'
ENGINE = OUT / 'actual01/budget_run'
sys.path.insert(0, str(ROOT / 'tests'))
from auditSU2V110 import vtk_tokens


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    target = OUT / 'checkpoint_progress.json'
    assert not target.exists()
    r = loadmat(str(OUT / 'review01/result.mat'), squeeze_me=True, struct_as_record=False)['r']
    V = r.V
    F = r.F.astype(int) - 1
    area = np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1) / 2
    centers = V[F].mean(axis=1)
    fields = []
    sources = [Path(__file__).resolve(), ROOT / 'tests/auditSU2V110.py', OUT / 'review01/result.mat']
    for step in [150, 299]:
        path = ENGINE / ('surface_flow_%06d.vtk' % step)
        points, data = vtk_tokens(path)
        dist, ids = cKDTree(V).query(points)
        assert max(dist) < 1e-6 and len(np.unique(ids)) == len(V)
        values = np.zeros((len(V), 3))
        values[ids] = np.column_stack([data[key][:, 0] for key in ['Pressure', 'Heat_Flux', 'Temperature']])
        fields.append(values[F].mean(axis=1))
        sources.append(path)
    old, new = fields

    def metrics(mask):
        weight = area[mask]
        a, b = old[mask], new[mask]
        # L1分母采用绝对值积分，避免正负热流相消制造小变化。
        return dict(faces=int(sum(mask)), areaM2=float(sum(weight)),
                    pressureRelativeL1Change=float(weight @ abs(b[:, 0]-a[:, 0]) / (weight @ abs(a[:, 0]))),
                    heatFluxRelativeL1Change=float(weight @ abs(b[:, 1]-a[:, 1]) / (weight @ abs(a[:, 1]))),
                    signedPowerW=[float(weight @ a[:, 1]), float(weight @ b[:, 1])])

    bins = np.linspace(min(V[:, 0]), max(V[:, 0]), 13)
    rows = []
    for i in range(12):
        mask = (centers[:, 0] >= bins[i]) & (centers[:, 0] <= bins[i+1] if i == 11 else centers[:, 0] < bins[i+1])
        if mask.any():
            rows.append(dict(xRangeM=bins[i:i+2].tolist(), **metrics(mask)))
    path = ENGINE / 'restart_flow_000299.csv'
    state = np.loadtxt(str(path), skiprows=1, delimiter=',')
    rho = state[:, 4]
    speed2 = np.sum((state[:, 5:8] / rho[:, None])**2, axis=1)
    T = (state[:, 8] / rho - speed2 / 2) * .4 / 287.05287
    cp = 1.4 * 287.05287 / .4
    T0 = T + speed2 / (2 * cp)
    nominalT0 = r.flight.temperature * (1 + .2 * r.flight.Mach**2)
    history = np.genfromtxt(str(ENGINE / 'history.csv'), names=True, delimiter=',')
    tail = history['HF'][-50:]
    sources += [path, ENGINE / 'history.csv']
    result = dict(physicsQualified=False, newQualifiedArea=0, noSolverRun=True,
                  checkpointSteps=[150, 299], entireSurface=metrics(np.ones(len(F), dtype=bool)),
                  axialBinsAreGeometricNotPhysicalRegions=True, axialBins=rows,
                  nominalFreestreamTotalTemperatureK=float(nominalT0),
                  staticTemperatureRangeK=[float(min(T)), float(max(T))],
                  totalTemperatureRangeK=[float(min(T0)), float(max(T0))],
                  nodesWithStaticTAboveFreestreamTotalT=int(sum(T > nominalT0)),
                  nodesWithTotalTAboveFreestreamTotalT1Percent=int(sum(T0 > nominalT0 * 1.01)),
                  final50NativeHeatPowerEndpointChange=float((tail[-1]-tail[0])/max(abs(tail[0]), 1)),
                  final50NativeHeatPowerPeakToPeakOverMean=float(np.ptp(tail)/max(abs(np.mean(tail)), 1)),
                  interpretation='startup evolution and overshoot diagnostics, not accuracy or physical convergence',
                  inputSHA256={str(p): sha(p) for p in sources})
    target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['axialBins','inputSHA256']}, indent=2))


if __name__ == '__main__':
    main()
