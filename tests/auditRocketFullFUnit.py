# -*- coding: utf-8 -*-
"""独立核对全部F系数与小域所有Arnoldi方向，固定C回代用原完整矩阵验算。"""
from pathlib import Path
import os, json, sys, hashlib
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
name = sys.argv[1]
assert name.replace('_', '').isalnum()
folder = VAL/name
target = VAL/(name+'_independent_full_F.json')
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


def rel(a, b):
    return float(np.linalg.norm(a-b)/max(1, np.linalg.norm(b)))


def dense(packed):
    n = int(packed['cellCount'])
    C = np.zeros((5*n, 5*n))
    ptr = np.atleast_1d(packed['rowPointer']).astype(int)-1
    col = np.atleast_1d(packed['column']).astype(int)-1
    order = np.atleast_1d(packed['order']).astype(int)-1
    values = np.asarray(packed['values']).reshape((5, 5, len(col)), order='F')
    for row in range(n):
        for edge in range(ptr[row], ptr[row+1]):
            rr, cc = order[row]+np.arange(5)*n, order[col[edge]]+np.arange(5)*n
            C[np.ix_(rr, cc)] = values[:, :, edge]
    return C


file = folder/'full_F_unit_evidence.mat'
s = loadmat(file, simplify_cells=True)
results = []
for case in rows(s['cases']):
    C = dense(case['packed'])
    split = case['split']
    n = int(split['cellCount'])
    f, c = split['fineOriginalCVs'].astype(int)-1, split['coarseOriginalCVs'].astype(int)-1
    ids = (f[:, None]+np.arange(5)[None, :]*n).ravel(order='F')
    A = dense(case['graph'])
    np.testing.assert_array_equal(A, C[np.ix_(ids, ids)])
    np.testing.assert_array_equal(A, case['A'])
    np.testing.assert_array_equal(case['graph']['originalPhysicalCVs'].astype(int)-1, f)
    rhs, x = case['rhs'].ravel(order='F'), case['value'].ravel(order='F')
    residual = float(np.linalg.norm(rhs-A@x)/np.linalg.norm(rhs))
    assert residual < 1.1e-10 and rel(x, np.linalg.solve(A, rhs)) < 1e-8
    d = case['detail']
    cycle = rows(d['krylov']['cycles'])[0]
    V, Z, H = cycle['V'], cycle['Z'], cycle['smallH']
    relation = rel(A@Z, V@H)
    mapping = rel((case['P']@V[:, :Z.shape[1]])/cycle['imageNormalizationScales'][None, :], Z)
    assert relation < 1e-8 and mapping < 1e-8
    full = case['complete'].ravel(order='F')
    parent = case['parentRhs'].ravel(order='F')
    image = C@full
    assert rel(image, case['parentImage'].ravel(order='F')) < 1e-8
    parent_residual = float(np.linalg.norm(parent-image)/np.linalg.norm(parent))
    assert parent_residual < 1e-8 and abs(parent_residual-case['parentResidual']) < 1e-8
    results.append(dict(originalCVs=n, FCVs=len(f), allFFCoefficientsExactlyPreserved=True,
        fullFResidual=residual, fixedCParentCResidual=parent_residual, allPhysicalArnoldiRelation=relation,
        allPhysicalLocalFMapDifference=mapping, arnoldiSteps=int(d['krylov']['arnoldiSteps'])))
out = dict(status='passed-independent-all-small-FF-coefficients-full-F-Arnoldi-and-fixed-C-parent-reference',
    processID=os.getpid(), cases=results, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):sha(file)}, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
