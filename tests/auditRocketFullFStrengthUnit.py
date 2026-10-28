# -*- coding: utf-8 -*-
"""独立稠密核对F内强连接缓存、完整F迭代方向及固定C回代。"""
from pathlib import Path
import os, json, hashlib, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
name = sys.argv[1]
assert name.replace('_', '').isalnum()
folder = VAL/name
target = VAL/(name+'_independent_full_F_strength.json')
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
    ptr, col, order = packed['rowPointer'].astype(int)-1, packed['column'].astype(int)-1, packed['order'].astype(int)-1
    for i in range(n):
        for edge in range(ptr[i], ptr[i+1]):
            rr, cc = order[i]+np.arange(5)*n, order[col[edge]]+np.arange(5)*n
            C[np.ix_(rr, cc)] = packed['values'][:, :, edge]
    return C


file = folder/'full_F_strength_unit_evidence.mat'
s = loadmat(file, simplify_cells=True)
results = []
for case in rows(s['cases']):
    A = dense(case['graph'])
    np.testing.assert_array_equal(A, case['A'])
    factor = case['factor']
    nf = int(factor['cellCount'])
    P = np.zeros_like(A)
    maximum = 0.
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(np.atleast_1d(bucket['patchIDs']))
        if not count:
            continue
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(int)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F')
        fallback = np.atleast_1d(bucket['fallback']).astype(bool)
        for page in range(count):
            ids = (members[:, page][None, :]+np.arange(5)[:, None]*nf).ravel(order='F')
            block = A[np.ix_(ids, ids)]
            if fallback[page]:
                diagonal = np.zeros_like(block)
                for cell in range(size):
                    ii = slice(5*cell, 5*cell+5)
                    diagonal[ii, ii] = block[ii, ii]
                block = diagonal
            expected = np.linalg.solve(block, np.eye(5*size))
            maximum = max(maximum, rel(expected, inverse[:, :, page]))
            assert maximum < 1e-8
            P[np.ix_(ids, ids)] = expected
    assert rel(P, case['P']) < 1e-8
    if 'coverage' in case:
        graph = case['graph']
        old_group = case['split']['fineGroupID'].astype(int)-1
        new_group = factor['partition']['groupID'].astype(int)-1
        fallback_group = np.zeros(int(factor['partition']['patchCount']), bool)
        for bucket in rows(factor['buckets']):
            fallback_group[np.atleast_1d(bucket['patchIDs']).astype(int)-1] = np.atleast_1d(bucket['fallback']).astype(bool)
        original_cross = included_cross = effective_count = inside_count = 0
        ptr = graph['rowPointer'].astype(int)-1
        column = graph['column'].astype(int)-1
        for row in range(nf):
            for edge in range(ptr[row], ptr[row+1]):
                q = column[edge]
                cross = old_group[row] != old_group[q]
                inside = new_group[row] == new_group[q]
                effective = inside and (not fallback_group[new_group[row]] or row == q)
                original_cross += cross
                inside_count += inside
                included_cross += cross and effective
                effective_count += effective
        coverage = case['coverage']
        assert original_cross == coverage['oldCrossPatchBlocks'] and included_cross == coverage['newlyIncludedOldCrossPatchBlocks']
        assert effective_count == coverage['newEffectiveLocalBlocksAfterFallback'] and inside_count == coverage['newWithinStrengthPatchBlocks']
    b, x = case['rhs'].ravel(order='F'), case['value'].ravel(order='F')
    assert rel(P@b, case['single'].ravel(order='F')) < 1e-8
    residual = float(np.linalg.norm(b-A@x)/np.linalg.norm(b))
    assert residual < 1.1e-10
    d = case['detail']
    cycle = rows(d['krylov']['cycles'])[0]
    V, Z, H = cycle['V'], cycle['Z'], cycle['smallH']
    relation = rel(A@Z, V@H)
    mapping = rel((P@V[:, :Z.shape[1]])/cycle['imageNormalizationScales'][None, :], Z)
    assert relation < 1e-8 and mapping < 1e-8
    C = dense(case['packed'])
    full, parent = case['complete'].ravel(order='F'), case['parentRhs'].ravel(order='F')
    parent_residual = float(np.linalg.norm(parent-C@full)/np.linalg.norm(parent))
    assert parent_residual < 1e-8
    results.append(dict(FCVs=nf, maximumLocalInverseReferenceDifference=maximum, fullFResidual=residual,
        allPhysicalFArnoldiRelation=relation, allPhysicalStrengthMapDifference=mapping,
        fixedCParentCResidual=parent_residual, arnoldiSteps=int(d['krylov']['arnoldiSteps'])))
out = dict(status='passed-independent-F-only-strength-inverses-all-small-F-directions-and-parent-reference', processID=os.getpid(),
    cases=results, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):sha(file)}, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
