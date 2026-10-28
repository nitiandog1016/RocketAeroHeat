# -*- coding: utf-8 -*-
"""独立核对MATLAB完整F/局部逆的全部小域转置列及真实F、父C最终像。"""
from pathlib import Path
import json
import hashlib
import os
import sys
import gc
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
mode, name = sys.argv[1:3]
assert mode in ('unit', 'real') and name.replace('_', '').isalnum()
folder = VAL / name
target = VAL / (name + '_independent_full_F_least_squares.json')
assert not target.exists()
bindings = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(part)
    return h.hexdigest()


def read(path, variables=None):
    path = Path(path)
    if not path.is_absolute():
        path = ROOT/path
    bindings[path.relative_to(ROOT).as_posix()] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def norm(value):
    return float(np.linalg.norm(value))


def relative(actual, expected):
    return norm(actual-expected)/max(1., norm(expected))


def dense(graph):
    n = int(graph['cellCount'])
    ptr = np.atleast_1d(graph['rowPointer']).astype(np.int64)-1
    col = np.atleast_1d(graph['column']).astype(np.int64)-1
    order = np.atleast_1d(graph['order']).astype(np.int64)-1
    values = np.asarray(graph['values']).reshape(5, 5, len(col), order='F')
    A = np.zeros((5*n, 5*n))
    for row in range(n):
        for edge in range(ptr[row], ptr[row+1]):
            rr, cc = order[row]+np.arange(5)*n, order[col[edge]]+np.arange(5)*n
            A[np.ix_(rr, cc)] = values[:, :, edge]
    return A


def local_dense(factor):
    n = int(factor['cellCount'])
    M = np.zeros((5*n, 5*n))
    buckets = [factor['buckets']] if isinstance(factor['buckets'], dict) else list(np.atleast_1d(factor['buckets']).ravel())
    for bucket in buckets:
        count = np.size(bucket['groupIDs'])
        if not count:
            continue
        k = np.size(bucket['members']) // count
        members = np.asarray(bucket['members']).reshape(k, count, order='F').astype(np.int64)-1
        inverse = np.asarray(bucket['inverse']).reshape(5*k, 5*k, count, order='F')
        for page in range(count):
            ids = (members[:, page, None] + np.arange(5)[None, :]*n).ravel()
            M[np.ix_(ids, ids)] = inverse[:, :, page]
    return M


def action(graph, value):
    n = int(graph['cellCount'])
    ptr = np.atleast_1d(graph['rowPointer']).astype(np.int64)-1
    col = np.atleast_1d(graph['column']).astype(np.int64)-1
    order = np.atleast_1d(graph['order']).astype(np.int64)-1
    owner = np.repeat(np.arange(n), np.diff(ptr))
    blocks = np.asarray(graph['values']).reshape(5, 5, len(col), order='F').transpose(2, 0, 1)
    ordered, output = value[order], np.zeros((n, 5))
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(output, owner[ii], np.einsum('kab,kb->ka', blocks[ii], ordered[col[ii]]))
    final = np.empty_like(output)
    final[order] = output
    return final


if mode == 'unit':
    s = read(folder / 'full_F_least_squares_unit_evidence.mat')
    results = []
    for case in s['cases']:
        A, M, PA = dense(case['graph']), local_dense(case['factor']), dense(case['permutedGraph'])
        errors = dict(F=relative(case['A'], A), FT=relative(case['AT'], A.T),
            M=relative(case['M'], M), MT=relative(case['MT'], M.T),
            permutedF=relative(case['permutedA'], PA), permutedFT=relative(case['permutedAT'], PA.T))
        assert max(errors.values()) < 1e-12
        direct = np.linalg.solve(A, case['rhs'].ravel(order='F'))
        error = relative(case['value'].ravel(order='F'), direct)
        residual = norm(case['rhs'].ravel(order='F')-A @ case['value'].ravel(order='F'))/norm(case['rhs'])
        assert error < 1e-8 and residual < 1e-10 and case['detail']['fullFInternalTargetMet']
        assert case['detail']['adjointIdentityRelativeDifference'] < 1e-12
        assert case['zeroDetail']['nativeIterations'] == case['earlyDetail']['nativeIterations'] == 0
        results.append(dict(allDenseMatrixRelativeDifferences=errors, fullFRelativeResidual=residual,
            solutionRelativeDifference=error, actualNativeIterations=int(case['detail']['nativeIterations'])))
    out = dict(status='passed-independent-all-small-F-local-inverse-transpose-columns-permutations-and-dense-LSQR',
        processID=os.getpid(), cases=results, newActualNSEvaluations=0, newProductionSolvedAreaM2=0)
else:
    s = read(folder / 'full_F_least_squares_real_evidence.mat')
    report = json.loads((folder / 'full_F_least_squares_summary.json').read_text(encoding='utf8'))
    for name, digest in s['sourceIdentity'].items():
        assert sha(ROOT / 'internalFlow/cutCell' / (name+'.m')) == digest
    assert sha(s['sourceBinding']['sourceFile']) == s['sourceBinding']['sourceFileSHA256']
    graph_file = Path(s['graphFile'])
    if not graph_file.is_absolute():
        graph_file = ROOT/graph_file
    prior = json.loads((VAL/'fullFreal01_independent_full_F.json').read_text())
    assert sha(graph_file) == prior['inputSHA256'][str(graph_file.relative_to(ROOT))]
    graph = read(graph_file, ['graph'])['graph']
    base = read(Path(s['baseFile']), ['valueF', 'value', 'rhsF', 'Fimage'])
    np.testing.assert_array_equal(s['initial'], base['valueF'])
    np.testing.assert_array_equal(s['rhsF'], base['rhsF'])
    image = action(graph, s['valueF'])
    error = relative(image, s['Fimage'])
    residual = norm(s['rhsF']-image)/norm(s['rhsF'])
    initial = norm(s['rhsF']-base['Fimage'])/norm(s['rhsF'])
    assert error < 1e-12 and abs(residual-s['detail']['relativeFullFResidual']) < 1e-12
    f = np.atleast_1d(graph['originalPhysicalCVs']).astype(np.int64)-1
    expected = base['value'].copy()
    expected[f] = s['valueF']
    np.testing.assert_array_equal(expected, s['value'])
    packed = read(Path(s['packedFile']), ['packed'])['packed']
    parent_image = action(packed, s['value'])
    parent_error = relative(parent_image, s['image'])
    remaining = s['originalRHS']-parent_image
    parent_residual = norm(remaining)/norm(s['originalRHS'])
    consistency = norm(remaining[f]-(s['rhsF']-image))/max(1., norm(s['originalRHS']))
    assert parent_error < 1e-12 and consistency < 1e-8
    assert abs(parent_residual-s['diagnostic']['relativeParentOriginalCResidual']) < 1e-10
    assert not report['physicalSolution'] and report['newActualNSEvaluations'] == report['newFlowAcceptances'] == 0
    assert not s['detail']['oldApproximateSchurIdentityReplaced']
    assert s['detail']['nativeIterations'] <= 32
    prototype_file = VAL/'fullFlsqrscreen01/full_F_lsqr_screen_summary.json'
    prototype = json.loads(prototype_file.read_text())
    out = dict(status='passed-independent-real-complete-F-LSQR-initial-final-image-and-fixed-C-parent-diagnostic',
        processID=os.getpid(), finalFullFImageDifference=error, parentOriginalCImageDifference=parent_error,
        initialRelativeFullFResidual=initial, finalRelativeFullFResidual=residual,
        independentScipyPrototypeFResidual=prototype['finalRelativeFullFResidual'],
        identicalNativeIterationTrajectoryClaimed=False,
        fixedCParentOriginalCResidual=parent_residual, FRowsSubmatrixConsistency=consistency,
        fullFResidualImproved=residual < initial, fullFInternalTargetMet=residual <= 1e-3,
        actualNativeIterations=int(s['detail']['nativeIterations']),
        nativeAdjointGatePassed=s['detail']['adjointIdentityRelativeDifference'] < 1e-8,
        independentAllRealBidiagonalizationColumnsClaimed=False,
        originalFFCoefficientsAuditReusedByExactFileIdentity=True, oldApproximateSchurIdentityClaimed=False,
        newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0, physicalSolution=False)

for file, digest in bindings.items():
    assert sha(ROOT / file) == digest
out.update(inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k: v for k, v in out.items() if k != 'inputSHA256'}))
