# -*- coding: utf-8 -*-
"""独立检查每个局部像的正面积两环支持及不触边界列的五量积分抵消。"""
from pathlib import Path
import hashlib, json, os, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix
ROOT = Path(__file__).resolve().parents[1]; VAL = ROOT / 'outputs/v83_validation'
case = sys.argv[1]; assert case.replace('_', '').isalnum()
source = VAL / case / 'local_coarse_evidence.mat'; target = VAL / (case + '_local_support.json')
assert not target.exists()
c = loadmat(source, simplify_cells=True, variable_names=['coarse'])['coarse']; cut = c['space']['geometry']
image = c['candidateImage'].tocsc(); cells = np.atleast_1d(c['candidateCellIDs']).astype(int) - 1
n = len(c['baseQ']); f = cut['openFaces']; m = f['cellMinus'].astype(int) - 1; p = f['cellPlus'].astype(int) - 1
inside = (f['area'] > 0) & (m >= 0) & (p >= 0)
graph = coo_matrix((np.ones(2 * inside.sum()), (np.r_[m[inside], p[inside]], np.r_[p[inside], m[inside]])), shape=(n, n)).tocsr()
boundary = (f['area'] > 0) & ((m >= 0) != (p >= 0))
boundary_cells = set(np.maximum(m[boundary], p[boundary]))
boundary_cells.update(np.atleast_1d(cut['walls']['ownerCell']).astype(int) - 1)
boundary_cells.discard(-1)
volume = cut['cells']['volume']; cache = {}; cancellation = []; nonzero_counts = []
for j, cell in enumerate(cells):
    if cell not in cache:
        first = np.unique(np.r_[cell, graph.indices[graph.indptr[cell]:graph.indptr[cell + 1]]])
        second = np.unique(np.r_[first, graph[first].nonzero()[1]])
        cache[cell] = set(second)
    support = cache[cell]; values = image[:, j]; row = values.indices
    active_cells = row % n
    assert set(active_cells) <= support, (j, cell, set(active_cells) - support)
    nonzero_counts.append(len(row))
    if not support & boundary_cells:
        flux = values.data * volume[active_cells]
        sums = np.bincount(row // n, weights=flux, minlength=5)
        magnitude = np.bincount(row // n, weights=np.abs(flux), minlength=5)
        relative = np.abs(sums) / np.maximum(magnitude, 1e-30)
        assert max(relative) < 1e-10
        cancellation.append(float(max(relative)))

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

out = dict(status='passed-independent-positive-area-two-ring-support-and-interior-five-equation-cancellation',
    case=case, originalCells=n, candidateColumns=len(cells), allCandidateColumnsChecked=True,
    minImageColumnNonzeros=min(nonzero_counts), maxImageColumnNonzeros=max(nonzero_counts),
    interiorColumnsWithNoBoundaryInTwoRings=len(cancellation),
    maximumInteriorFiveEquationRelativeCancellation=max(cancellation, default=None),
    actualOriginalNSEvaluations=0, modelApplications=0, physicalSolution=False, newProductionSolvedAreaM2=0,
    processID=os.getpid(), inputSHA256=sha(source), scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2) + '\n', encoding='utf8'); print(json.dumps(out))
