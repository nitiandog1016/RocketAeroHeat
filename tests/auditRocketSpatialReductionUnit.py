# -*- coding: utf-8 -*-
"""独立稠密核对空间C/F全部系数、实际局部F映射、近似Schur和完整Arnoldi。"""
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
target = VAL/(name+'_independent_spatial_reduction.json')
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def vec(value, dtype=float):
    return np.atleast_1d(value).astype(dtype)


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


def rel(a, b):
    return float(np.linalg.norm(a-b)/max(1, np.linalg.norm(b)))


def dense_packed(packed):
    n = int(packed['cellCount'])
    matrix = np.zeros((5*n, 5*n))
    ptr, col, order = vec(packed['rowPointer'], int)-1, vec(packed['column'], int)-1, vec(packed['order'], int)-1
    values = np.asarray(packed['values']).reshape((5, 5, len(col)), order='F')
    for row in range(n):
        for edge in range(ptr[row], ptr[row+1]):
            rr, cc = order[row]+np.arange(5)*n, order[col[edge]]+np.arange(5)*n
            matrix[np.ix_(rr, cc)] = values[:, :, edge]
    return matrix


def dense_rect(coupling):
    nr, nc = int(coupling['rowCount']), int(coupling['columnCount'])
    owner, column = vec(coupling['owner'], int)-1, vec(coupling['column'], int)-1
    values = np.asarray(coupling['values']).reshape((5, 5, len(column)), order='F')
    matrix = np.zeros((5*nr, 5*nc))
    for edge in range(len(column)):
        rr, cc = owner[edge]+np.arange(5)*nr, column[edge]+np.arange(5)*nc
        matrix[np.ix_(rr, cc)] = values[:, :, edge]
    return matrix


def local_inverse(factor, original, key):
    n = int(factor['cellCount'])
    matrix = np.zeros((5*n, 5*n))
    maximum = 0.
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(vec(bucket[key]))
        if not count:
            continue
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(int)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F')
        fallback = vec(bucket['fallback'], bool)
        for page in range(count):
            ids = (members[:, page][None, :]+np.arange(5)[:, None]*n).ravel(order='F')
            block = original[np.ix_(ids, ids)]
            if fallback[page]:
                diagonal = np.zeros_like(block)
                for cell in range(size):
                    ii = slice(5*cell, 5*cell+5)
                    diagonal[ii, ii] = block[ii, ii]
                block = diagonal
            expected = np.linalg.solve(block, np.eye(5*size))
            maximum = max(maximum, rel(expected, inverse[:, :, page]))
            assert maximum < 1e-8
            matrix[np.ix_(ids, ids)] = expected
    return matrix, maximum


file = folder/'spatial_reduction_unit_evidence.mat'
s = loadmat(file, simplify_cells=True)
results = []
for case in rows(s['cases']):
    C = dense_packed(case['packed'])
    reduction = case['reduction']
    split = reduction['split']
    n = int(split['cellCount'])
    f, c = vec(split['fineOriginalCVs'], int)-1, vec(split['coarseOriginalCVs'], int)-1
    np.testing.assert_array_equal(np.sort(np.r_[f, c]), np.arange(n))
    fi = (f[:, None]+np.arange(5)[None, :]*n).ravel(order='F')
    ci = (c[:, None]+np.arange(5)[None, :]*n).ravel(order='F')
    A, B, E, D = C[np.ix_(fi, fi)], C[np.ix_(fi, ci)], C[np.ix_(ci, fi)], C[np.ix_(ci, ci)]
    for key, expected in [('B', B), ('E', E), ('D', D)]:
        np.testing.assert_array_equal(dense_rect(reduction[key]), expected)
    MF, f_error = local_inverse(reduction['Ffactor'], A, 'groupIDs')
    assert rel(MF, case['MF']) < 1e-8
    S = D-E@MF@B
    schur_error = rel(dense_packed(reduction['schur']), S)
    assert schur_error < 1e-12
    MS, coarse_error = local_inverse(case['coarseFactor'], S, 'patchIDs')
    M = np.zeros_like(C)
    M[np.ix_(fi, fi)] = MF+MF@B@MS@E@MF
    M[np.ix_(fi, ci)] = -MF@B@MS
    M[np.ix_(ci, fi)] = -MS@E@MF
    M[np.ix_(ci, ci)] = MS
    mapping = rel(M, case['M'])
    assert mapping < 1e-9 and rel(M@case['rhs'].ravel(order='F'), case['single'].ravel(order='F')) < 1e-9
    b, x = case['rhs'].ravel(order='F'), case['value'].ravel(order='F')
    residual = float(np.linalg.norm(b-C@x)/np.linalg.norm(b))
    assert residual < 1.1e-10 and rel(x, np.linalg.solve(C, b)) < 1e-8
    d = case['detail']
    cycle = rows(d['krylov']['cycles'])[0]
    Z, V, H = cycle['Z'], cycle['V'], cycle['smallH']
    relation = rel(C@Z, V@H)
    maps = rel((M@V[:, :Z.shape[1]])/cycle['imageNormalizationScales'][None, :], Z)
    assert relation < 1e-8 and maps < 1e-8
    steps = int(d['krylov']['arnoldiSteps'])
    assert d['FLocalApplications'] == 2*steps and d['coarseLocalApplications'] == steps and d['rectangularCouplingApplications'] == 2*steps
    # 局部F映射与全F逆的差保持可见，不把近似分块当精确Schur。
    full_F_inverse_difference = rel(MF, np.linalg.solve(A, np.eye(len(fi))))
    assert not reduction['constructionReport']['fullFFInverseOrExactSchurClaimed']
    results.append(dict(originalCVs=n, coarseCVs=len(c), FCVs=len(f), maximumFInverseReferenceDifference=f_error,
        maximumCoarseInverseReferenceDifference=coarse_error, approximateSchurAllCoefficientDifference=schur_error,
        fullBlockMapDifference=mapping, localFMapToCompleteFInverseDifference=full_F_inverse_difference,
        fullCResidual=residual, allPhysicalArnoldiRelationDifference=relation, allPhysicalSpatialMapDifference=maps, arnoldiSteps=steps))
fallback = s['fallback']
np.testing.assert_array_equal(dense_packed(fallback['schur']), -np.eye(5))
assert fallback['Ffactor']['factorizationReport']['fallbackFBlockCount'] == 1
out = dict(status='passed-independent-spatial-five-field-split-all-coefficients-block-map-and-all-small-Arnoldi-directions',
    processID=os.getpid(), cases=results, explicitSingularLocalFFallbackVerified=True, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):sha(file)}, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
