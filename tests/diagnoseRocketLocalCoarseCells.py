# -*- coding: utf-8 -*-
"""定位保存原NS方向剩余的原控制体和正面积邻接；不求解或修改原面/物理方程。"""
from pathlib import Path
import hashlib, json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
case, label = sys.argv[1:3]
assert all(x.replace('_', '').isalnum() for x in (case, label))
folder = VAL / label
assert not folder.exists()
folder.mkdir()

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

source = VAL / case / 'wall_velocity_krylov_evidence.mat'
r = loadmat(source, simplify_cells=True, variable_names=['result'])['result']
logs = r['linearLog']
log = logs if isinstance(logs, dict) else list(logs)[-1]
remaining = r['initialRate'] / r['residualScale'] + log['steadyDirectionEvidence']['value']
coarse_file = Path(r['options']['coarseOperatorFile'])
if not coarse_file.is_absolute():
    coarse_file = ROOT / coarse_file
coarse = loadmat(coarse_file, simplify_cells=True, variable_names=['coarse'])['coarse']
cut = coarse['space']['geometry']
assert np.array_equal(r['Q'], coarse['baseQ'])
n = len(r['Q']); nf = len(cut['inputGeometry']['F']); nl = n // nf
assert n == nf * nl
energy = np.sum(remaining ** 2, axis=1)
order = np.lexsort((np.arange(n), -energy))
faces = cut['openFaces']
m = faces['cellMinus'].astype(int) - 1
p = faces['cellPlus'].astype(int) - 1
inside = (faces['area'] > 0) & (m >= 0) & (p >= 0)
graph = coo_matrix((np.ones(2 * inside.sum()),
    (np.r_[m[inside], p[inside]], np.r_[p[inside], m[inside]])), shape=(n, n)).tocsr()
norm2 = float(energy.sum())
fractions = {}
for count in (8, 16, 32, 64, 128, 256, 512, 1024):
    selected = order[:min(n, count)]
    near = np.unique(np.r_[selected, graph[selected].nonzero()[1]])
    fractions[str(count)] = dict(selectedCells=len(selected),
        selectedRemainingSquaredNormFraction=float(energy[selected].sum() / norm2),
        selectedAndPositiveFaceNeighbours=len(near),
        neighbourhoodRemainingSquaredNormFraction=float(energy[near].sum() / norm2))

def entry(i):
    adjacent = graph.indices[graph.indptr[i]:graph.indptr[i + 1]]
    return dict(originalCellID=int(i + 1), originalFaceID=int(i % nf + 1),
        layerID=int(i // nf + 1), centroidM=cut['cells']['centroid'][i].tolist(),
        volumeM3=float(cut['cells']['volume'][i]),
        remainingSquaredNormFraction=float(energy[i] / norm2),
        positiveAreaNeighbourCellIDs=(adjacent + 1).tolist(),
        positiveAreaNeighbourOriginalFaceIDs=(adjacent % nf + 1).tolist())

report = dict(status='mapped-original-CV-residual-and-positive-area-neighbourhood',
    sourceCase=case, originalCellCount=n, originalFaceCount=nf, layers=nl,
    actualNSRelativeDefect=float(np.linalg.norm(remaining) / np.linalg.norm(r['initialRate'] / r['residualScale'])),
    dominantCellCoverage=fractions, dominantOriginalCells=[entry(i) for i in order[:64]],
    newModelApplications=0, newActualNSEvaluations=0, physicalSolution=False,
    newProductionSolvedAreaM2=0, processID=os.getpid(),
    inputSHA256={str(f.relative_to(ROOT)): sha(f) for f in (source, coarse_file)},
    scriptSHA256=sha(__file__))
assert all(sha(ROOT / path) == digest for path, digest in report['inputSHA256'].items())
(folder / 'local_cells_diagnostic.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
savemat(folder / 'local_cells_selection.mat', dict(cellOrder=order[:min(n, 1024)] + 1,
    remaining=remaining, sourceFile=str(source), coarseFile=str(coarse_file),
    originalFaceCount=nf, originalCellCount=n, noPhysicalMeshChange=True))
print(json.dumps(dict(status=report['status'], topCells=report['dominantOriginalCells'][:3], coverage=fractions)))
