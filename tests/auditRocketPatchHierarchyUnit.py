# -*- coding: utf-8 -*-
"""独立稠密核对两排序的完整RCP、每层五量逆、V循环及全部Arnoldi方向。"""
from pathlib import Path
import os, json, sys, hashlib
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
name = sys.argv[1]
assert name.replace('_', '').isalnum()
folder = VAL / name
target = VAL / (name+'_independent_patch_hierarchy.json')
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(8*1024**2), b''):
            h.update(part)
    return h.hexdigest()


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


def rel(value, reference):
    return float(np.linalg.norm(value-reference)/max(1, np.linalg.norm(reference)))


def dense_graph(packed):
    n = int(packed['cellCount'])
    C = np.zeros((5*n, 5*n))
    order = np.atleast_1d(packed['order']).astype(int)-1
    ptr = np.atleast_1d(packed['rowPointer']).astype(int)-1
    col = np.atleast_1d(packed['column']).astype(int)-1
    values = np.asarray(packed['values']).reshape((5, 5, len(col)), order='F')
    for i in range(n):
        rr = order[i]+np.arange(5)*n
        for edge in range(ptr[i], ptr[i+1]):
            cc = order[col[edge]]+np.arange(5)*n
            C[np.ix_(rr, cc)] = values[:, :, edge]
    return C


file = folder / 'patch_hierarchy_unit_evidence.mat'
s = loadmat(file, simplify_cells=True)
results = []
for case in rows(s['cases']):
    levels = rows(case['hierarchy']['levels'])
    matrices = [dense_graph(level['packed']) for level in levels]
    maximum_coarse = 0.
    maximum_local_inverse = 0.
    current = None
    for j in range(len(levels)-1, -1, -1):
        level = levels[j]
        C = matrices[j]
        n = C.shape[0]//5
        if level['terminal']:
            current = level['inverse']
            defect = rel(C@current, np.eye(5*n))
            assert defect < 1e-8 and rel(current, np.linalg.solve(C, np.eye(5*n))) < 1e-8
            continue
        factor = level['factor']
        part = factor['partition']
        ng = int(part['patchCount'])
        group = part['groupID'].astype(int)-1
        P = np.eye(ng)[group]
        R = P.T/np.atleast_1d(part['sizes'])[:, None]
        P5 = np.kron(np.eye(5), P)
        R5 = np.kron(np.eye(5), R)
        assert rel(R5@P5, np.eye(5*ng)) < 1e-13
        coarse_error = rel(R5@C@P5, matrices[j+1])
        assert coarse_error < 1e-12
        maximum_coarse = max(maximum_coarse, coarse_error)
        local = np.zeros_like(C)
        for size, bucket in enumerate(rows(factor['buckets']), 1):
            count = len(np.atleast_1d(bucket['patchIDs']))
            if not count:
                continue
            members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(int)-1
            inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F')
            fallback = np.atleast_1d(bucket['fallback']).astype(bool)
            for page in range(count):
                ids = (members[:, page][None, :]+np.arange(5)[:, None]*n).ravel(order='F')
                block = C[np.ix_(ids, ids)]
                if fallback[page]:
                    diagonal = np.zeros_like(block)
                    for cell in range(size):
                        ii = slice(5*cell, 5*cell+5)
                        diagonal[ii, ii] = block[ii, ii]
                    block = diagonal
                expected = np.linalg.solve(block, np.eye(5*size))
                maximum_local_inverse = max(maximum_local_inverse, rel(inverse[:, :, page], expected))
                assert maximum_local_inverse < 1e-8
                local[np.ix_(ids, ids)] = expected
        current = local+P5@current@R5@(np.eye(5*n)-C@local)
    C = matrices[0]
    b = case['rhs'].ravel(order='F')
    mapping = rel(current, case['M'])
    assert mapping < 1e-9 and rel(current@b, case['single'].ravel(order='F')) < 1e-9
    candidate = case['value'].ravel(order='F')
    residual = float(np.linalg.norm(b-C@candidate)/np.linalg.norm(b))
    assert residual < 1.1e-10 and rel(candidate, np.linalg.solve(C, b)) < 1e-8
    detail = case['detail']
    cycle = rows(detail['krylov']['cycles'])[0]
    Z, V, H = cycle['Z'], cycle['V'], cycle['smallH']
    relation = rel(C@Z, V@H)
    map_error = rel((current@V[:, :Z.shape[1]])/cycle['imageNormalizationScales'][None, :], Z)
    assert relation < 1e-8 and map_error < 1e-8
    assert detail['totalRootCApplications'] == detail['solverCApplications']+detail['VcycleApplications']
    results.append(dict(levels=len(levels), maximumAllCoefficientRCPDifference=maximum_coarse,
        maximumLocalInverseReferenceDifference=maximum_local_inverse, fixedDenseVcycleMapDifference=mapping,
        fullCResidual=residual, allPhysicalArnoldiRelationDifference=relation,
        allPhysicalVcycleMapDifference=map_error, arnoldiSteps=int(detail['krylov']['arnoldiSteps'])))
    if 'deepHierarchy' in case:
        deep_levels = rows(case['deepHierarchy']['levels'])
        assert len(deep_levels) >= 3
        current = None
        maximum_deep_coarse = 0.
        for j in range(len(deep_levels)-1, -1, -1):
            level = deep_levels[j]
            C = dense_graph(level['packed'])
            n = C.shape[0]//5
            if level['terminal']:
                current = np.linalg.solve(C, np.eye(5*n))
                assert rel(current, level['inverse']) < 1e-8
                continue
            part = level['factor']['partition']
            P = np.eye(int(part['patchCount']))[part['groupID'].astype(int)-1]
            R = P.T/np.atleast_1d(part['sizes'])[:, None]
            P5, R5 = np.kron(np.eye(5), P), np.kron(np.eye(5), R)
            coarse_error = rel(R5@C@P5, dense_graph(deep_levels[j+1]['packed']))
            assert coarse_error < 1e-12
            maximum_deep_coarse = max(maximum_deep_coarse, coarse_error)
            local = np.zeros_like(C)
            for size, bucket in enumerate(rows(level['factor']['buckets']), 1):
                count = len(np.atleast_1d(bucket['patchIDs']))
                if not count:
                    continue
                members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(int)-1
                inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F')
                assert not np.any(bucket['fallback'])
                for page in range(count):
                    ids = (members[:, page][None, :]+np.arange(5)[:, None]*n).ravel(order='F')
                    block = C[np.ix_(ids, ids)]
                    expected = np.linalg.solve(block, np.eye(5*size))
                    assert rel(expected, inverse[:, :, page]) < 1e-8
                    local[np.ix_(ids, ids)] = expected
            current = local+P5@current@R5@(np.eye(5*n)-C@local)
        deep_map = rel(current, case['deepMatrix'])
        assert deep_map < 1e-9 and rel(current@b, case['deepSingle'].ravel(order='F')) < 1e-9
        work = case['deepWork']
        np.testing.assert_array_equal(work['CApplicationsByLevel'], [1]*(len(deep_levels)-1)+[0])
        results[-1].update(deepRecursiveLevels=len(deep_levels), maximumDeepRCPDifference=maximum_deep_coarse,
            deepFixedDenseVcycleMapDifference=deep_map, deepLevelWorkAccountingVerified=True)
if s['report'].get('singularTerminalRejected'):
    rejected = loadmat(folder/'singular_terminal_evidence.mat', simplify_cells=True)
    assert not rejected['certificate']['usable'] and np.linalg.matrix_rank(rejected['dense']) < rejected['dense'].shape[0]
out = dict(status='passed-independent-all-small-RCP-levels-fixed-V-cycle-and-all-physical-Arnoldi-columns',
    processID=os.getpid(), cases=results, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):sha(file)}, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
