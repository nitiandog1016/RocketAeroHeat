# -*- coding: utf-8 -*-
"""独立稠密核对扩展全部逆、仅核心行映射、零重叠及全小域F方向。"""
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
target = VAL/(name+'_independent_restricted_overlap.json')
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


def vec(value, dtype=float):
    return np.atleast_1d(value).astype(dtype)


def rel(a, b):
    return float(np.linalg.norm(a-b)/max(1, np.linalg.norm(b)))


def dense(graph):
    n = int(graph['cellCount'])
    C = np.zeros((5*n, 5*n))
    ptr, col, order = vec(graph['rowPointer'], int)-1, vec(graph['column'], int)-1, vec(graph['order'], int)-1
    values = np.asarray(graph['values']).reshape((5, 5, len(col)), order='F')
    for row in range(n):
        for edge in range(ptr[row], ptr[row+1]):
            rr, cc = order[row]+np.arange(5)*n, order[col[edge]]+np.arange(5)*n
            C[np.ix_(rr, cc)] = values[:, :, edge]
    return C


def map_reference(factor, A):
    n = int(factor['cellCount'])
    P = np.zeros((5*n, 5*n))
    core_coverage = []
    maximum = 0.
    for bucket in rows(factor['buckets']):
        count = len(vec(bucket['groupIDs']))
        if not count:
            continue
        m, k = int(bucket['coreSize']), int(bucket['overlapSize'])
        members = np.asarray(bucket['members']).reshape((k, count), order='F').astype(int)-1
        core = np.asarray(bucket['coreMembers']).reshape((m, count), order='F').astype(int)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*k, 5*k, count), order='F')
        fallback = vec(bucket['fallback'], bool)
        for page in range(count):
            np.testing.assert_array_equal(members[:m, page], core[:, page])
            assert len(np.unique(members[:, page])) == k
            ids = (members[:, page][None, :]+np.arange(5)[:, None]*n).ravel(order='F')
            rr = (core[:, page][None, :]+np.arange(5)[:, None]*n).ravel(order='F')
            block = A[np.ix_(ids, ids)]
            if fallback[page]:
                diagonal = np.zeros_like(block)
                for cell in range(k):
                    ii = slice(5*cell, 5*cell+5)
                    diagonal[ii, ii] = block[ii, ii]
                block = diagonal
            expected = np.linalg.solve(block, np.eye(5*k))
            maximum = max(maximum, rel(expected, inverse[:, :, page]))
            assert maximum < 1e-8
            P[np.ix_(rr, ids)] = expected[:5*m]
            core_coverage.extend(core[:, page].tolist())
    np.testing.assert_array_equal(np.sort(core_coverage), np.arange(n))
    return P, maximum


file = folder/'restricted_overlap_unit_evidence.mat'
s = loadmat(file, simplify_cells=True)
results = []
for case in rows(s['cases']):
    A = dense(case['graph'])
    np.testing.assert_array_equal(A, case['A'])
    P, inverse_error = map_reference(case['factor'], A)
    plain, plain_error = map_reference(case['plainFactor'], A)
    assert rel(P, case['P']) < 1e-8 and rel(plain@case['rhs'].ravel(order='F'), case['plainValue'].ravel(order='F')) < 1e-8
    b, x = case['rhs'].ravel(order='F'), case['value'].ravel(order='F')
    assert rel(P@b, case['single'].ravel(order='F')) < 1e-8
    residual = float(np.linalg.norm(b-A@x)/np.linalg.norm(b))
    assert residual < 1.1e-10 and rel(x, np.linalg.solve(A, b)) < 1e-8
    d = case['detail']
    cycle = rows(d['krylov']['cycles'])[0]
    V, Z, H = cycle['V'], cycle['Z'], cycle['smallH']
    relation = rel(A@Z, V@H)
    mapping = rel((P@V[:, :Z.shape[1]])/cycle['imageNormalizationScales'][None, :], Z)
    assert relation < 1e-8 and mapping < 1e-8
    C = dense(case['packed'])
    parent, full = case['parentRhs'].ravel(order='F'), case['complete'].ravel(order='F')
    parent_residual = float(np.linalg.norm(parent-C@full)/np.linalg.norm(parent))
    assert parent_residual < 1e-8 and rel(C@full, case['parentImage'].ravel(order='F')) < 1e-8
    results.append(dict(FCVs=int(case['graph']['cellCount']), maximumFullExtendedInverseReferenceDifference=inverse_error,
        maximumZeroOverlapInverseReferenceDifference=plain_error, allCoreOutputsExactlyOnce=True,
        fullFResidual=residual, fixedCParentCResidual=parent_residual, allPhysicalArnoldiRelation=relation,
        allPhysicalRestrictedCoreMapDifference=mapping, arnoldiSteps=int(d['krylov']['arnoldiSteps'])))
toy = dense(s['toy'])
toy_map, _ = map_reference(s['fallback'], toy)
np.testing.assert_array_equal(toy_map, np.eye(10))
np.testing.assert_array_equal(s['toyValue'], s['toyRhs'])
np.testing.assert_array_equal((toy@s['toyValue'].ravel(order='F')).reshape((2, 5), order='F'), s['toyImage'])
assert s['fallback']['factorizationReport']['fallbackOverlapBlocks'] == 2
out = dict(status='passed-independent-all-small-extended-inverses-core-restriction-and-complete-F-directions', processID=os.getpid(),
    cases=results, explicitSingularExtendedFallbackKeepsOriginalF=True, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):sha(file)}, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
