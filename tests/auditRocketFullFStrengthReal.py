# -*- coding: utf-8 -*-
"""独立核对真实F内强连接全部缓存逆、有效连接计数及完整F/固定C诊断。"""
from pathlib import Path
from decimal import Decimal, localcontext
import os, json, sys, hashlib, gc
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
bindings = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path, variables=None):
    path = Path(path)
    path = path if path.is_absolute() else ROOT/path
    bindings[str(path.relative_to(ROOT))] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def vec(value, dtype=float):
    return np.atleast_1d(value).astype(dtype)


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


def norm(value):
    scale = float(np.max(abs(value))) if np.size(value) else 0.
    return float(np.linalg.norm(value/scale)*scale) if scale else 0.


def rel(value, reference):
    return norm(value-reference)/max(1, norm(reference))


def unpack(packed):
    n = int(packed['cellCount'])
    ptr, col, order = vec(packed['rowPointer'], np.int64)-1, vec(packed['column'], np.int64)-1, vec(packed['order'], np.int64)-1
    owner = np.repeat(np.arange(n, dtype=np.int32), np.diff(ptr))
    values = np.asarray(packed['values']).reshape((5, 5, len(col)), order='F').transpose(2, 0, 1)
    return n, col, order, owner, values


def action(packed, vector):
    n, col, order, owner, values = unpack(packed)
    result = np.zeros((n, 5))
    ordered = vector[order]
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(result, owner[ii], np.einsum('kab,kb->ka', values[ii], ordered[col[ii]]))
    output = np.empty_like(result)
    output[order] = result
    return output


def precondition(factor, rhs):
    output = np.zeros_like(rhs)
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(vec(bucket['patchIDs']))
        if not count:
            continue
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        ids = members.T
        output[ids] = np.einsum('pij,pj->pi', inverse, rhs[ids].reshape(count, 5*size)).reshape(count, size, 5)
    return output


summary_file = folder/'full_F_strength_summary.json'
assert summary_file.stat().st_size < 128*1024
summary = json.loads(summary_file.read_text(encoding='utf-8-sig'))
assert summary['innerRunComplete'] and summary['parentDiagnosticComplete'] and not summary['candidateRejected']
s = read(folder/'full_F_strength_real_evidence.mat')
fdata = read(s['factorFile'])
factor = fdata['factor']
assert s['sourceBinding'] == fdata['sourceBinding'] == summary['sourceBinding']
for key, digest in s['sourceIdentity'].items():
    assert sha(ROOT/'internalFlow/cutCell'/(key+'.m')) == digest
graph = read(s['graphFile'], ['graph'])['graph']
old = json.loads((VAL/'fullFreal01_independent_full_F.json').read_text())
assert old['status'].startswith('passed-')
graph_key = str(Path(s['graphFile']).relative_to(ROOT))
assert bindings[graph_key] == old['inputSHA256'][graph_key]
nf, column, order, owner, values = unpack(graph)
part = factor['partition']
group = vec(part['groupID'], np.int64)-1
members = vec(part['members'], np.int64)-1
pointer = vec(part['groupPointer'], np.int64)-1
counts = vec(part['sizes'], np.int64)
ng = int(part['patchCount'])
np.testing.assert_array_equal(np.sort(members), np.arange(nf))
np.testing.assert_array_equal(group[members], np.repeat(np.arange(ng), counts))
position = np.empty(nf, dtype=int)
position[members] = np.arange(nf)-np.repeat(pointer[:-1], counts)
row_cv, col_cv = order[owner], order[column]
same = group[row_cv] == group[col_cv]
selected = np.flatnonzero(same)
selected = selected[np.argsort(group[row_cv[selected]], kind='stable')]
edge_ptr = np.r_[0, np.cumsum(np.bincount(group[row_cv[selected]], minlength=ng))]
maximum_ordinary = maximum_exact = maximum_source = 0.
decimal_count = fallback_count = checked = 0
fallback_group = np.zeros(ng, bool)
for size, bucket in enumerate(rows(factor['buckets']), 1):
    groups = vec(bucket['patchIDs'], np.int64)-1
    count = len(groups)
    if not count:
        continue
    inv = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
    mem = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
    bounds, fallback = vec(bucket['certifiedBounds']), vec(bucket['fallback'], bool)
    fallback_group[groups] = fallback
    assert np.all(np.isfinite(bounds) & (bounds < 1e-8))
    page_of_group = np.full(ng, -1, dtype=int)
    page_of_group[groups] = np.arange(count)
    N = 5*size
    gamma = N*np.finfo(float).eps/(1-N*np.finfo(float).eps)
    for first in range(0, count, 128):
        last = min(count, first+128)
        batch = groups[first:last]
        edges = np.concatenate([selected[edge_ptr[g]:edge_ptr[g+1]] for g in batch])
        page = page_of_group[group[row_cv[edges]]]-first
        rr, cc = position[row_cv[edges]], position[col_cv[edges]]
        block = np.zeros((last-first, N, N))
        adjacent = np.zeros((last-first, size, size), bool)
        adjacent[page, rr, cc] = np.any(values[edges] != 0, axis=(1, 2))
        adjacent |= adjacent.transpose(0, 2, 1)
        for a in range(5):
            for b in range(5):
                block[page, 5*rr+a, 5*cc+b] = values[edges, a, b]
        for local, g in enumerate(batch):
            np.testing.assert_array_equal(mem[:, first+local], members[pointer[g]:pointer[g+1]])
            visited = {0}
            while True:
                added = {q for q in range(size) if any(adjacent[local, q, k] for k in visited)}-visited
                if not added:
                    break
                visited |= added
            assert len(visited) == size
            if fallback[first+local]:
                diagonal = np.zeros((N, N))
                for cell in range(size):
                    ii = slice(5*cell, 5*cell+5)
                    diagonal[ii, ii] = block[local, ii, ii]
                block[local] = diagonal
                fallback_count += 1
        inverse = inv[first:last]
        ordinary = np.max(abs(block@inverse-np.eye(N)), axis=(1, 2))
        rounding = gamma*np.max(np.sum(abs(block), axis=2), axis=1)*np.max(abs(inverse), axis=(1, 2))/(1-gamma)**2
        maximum_ordinary = max(maximum_ordinary, float(ordinary.max()))
        maximum_source = max(maximum_source, float(bounds[first:last].max()))
        with localcontext() as context:
            context.prec = 80
            for local in np.flatnonzero(ordinary+rounding >= 1e-8):
                AA = [[Decimal.from_float(float(x)) for x in row] for row in block[local]]
                XX = [[Decimal.from_float(float(x)) for x in row] for row in inverse[local]]
                exact = max(abs((Decimal(1) if i == j else Decimal(0))-sum((AA[i][q]*XX[q][j] for q in range(N)), Decimal(0))) for i in range(N) for j in range(N))
                assert exact < Decimal('1e-8') and exact <= Decimal.from_float(float(bounds[first+local]))
                maximum_exact = max(maximum_exact, float(exact))
                decimal_count += 1
        checked += last-first
    print('INDEPENDENT_F_STRENGTH_BUCKET '+str(size), flush=True)
assert checked == ng and fallback_count == factor['factorizationReport']['fallbackPatchCount']
topology = read(VAL/'spatialreductionpreflight01/spatial_reduction_input.mat', ['originalGroupIDOneBased'])
old_group = np.asarray(topology['originalGroupIDOneBased']).ravel()[vec(graph['originalPhysicalCVs'], np.int64)-1]
old_cross = old_group[row_cv] != old_group[col_cv]
effective = same & (~fallback_group[group[row_cv]] | (row_cv == col_cv))
coverage = fdata['coverage']
assert np.count_nonzero(old_cross) == coverage['oldCrossPatchBlocks'] == 682056
assert np.count_nonzero(effective) == coverage['newEffectiveLocalBlocksAfterFallback']
assert np.count_nonzero(old_cross & effective) == coverage['newlyIncludedOldCrossPatchBlocks']
assert np.count_nonzero(same) == coverage['newWithinStrengthPatchBlocks']
del selected, values, row_cv, col_cv, same, old_cross, effective, block
gc.collect()
initial = read(s['inputFile'], ['rhsF', 'originalRHS'])
base = read(s['baseFile'], ['value'])
np.testing.assert_array_equal(initial['rhsF'], s['rhsF'])
np.testing.assert_array_equal(initial['originalRHS'], s['originalRHS'])
expected = precondition(factor, s['rhsF'])
map_error = rel(expected, s['single'])
single_image, final_image = action(graph, s['single']), action(graph, s['valueF'])
assert map_error < 1e-8 and rel(single_image, s['singleImage']) < 1e-8 and rel(final_image, s['Fimage']) < 1e-8
single_residual, residual = norm(s['rhsF']-single_image)/norm(s['rhsF']), norm(s['rhsF']-final_image)/norm(s['rhsF'])
detail = s['detail']
assert abs(residual-detail['relativeCResidual']) < 1e-8
f = vec(graph['originalPhysicalCVs'], np.int64)-1
parent_value = base['value'].copy()
parent_value[f] = s['valueF']
np.testing.assert_array_equal(parent_value, s['value'])
packed = read(s['packedFile'], ['packed'])['packed']
parent_image = action(packed, parent_value)
parent_difference = rel(parent_image, s['image'])
assert parent_difference < 1e-8
remaining = s['originalRHS']-parent_image
parent_residual = norm(remaining)/norm(s['originalRHS'])
identity = norm(remaining[f]-(s['rhsF']-final_image))/max(1, norm(s['originalRHS']))
assert identity < 1e-8 and abs(parent_residual-s['diagnostic']['relativeParentOriginalCResidual']) < 1e-8
assert not s['diagnostic']['oldApproximateSchurIdentityClaimedForNewMap']
for cycle in rows(detail['krylov']['cycles']):
    assert 'V' not in cycle and 'Z' not in cycle and cycle['basisOrthogonalityError'] < 1e-8 and cycle['arnoldiRelationRelativeDifference'] < 1e-8
assert not summary['fullFArnoldiBasisSaved'] and bool(summary['fullFInternalTargetMet']) == bool(detail['krylov']['flag'] == 0)
assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-real-F-strength-all-local-inverses-coverage-full-F-and-fixed-C-parent-diagnostic', processID=os.getpid(),
    FCVs=nf, strengthPatches=ng, allOriginalFUnknownsExactlyOnce=True, allStrengthPatchesConnected=True,
    originalFFCoefficientAuditReusedByExactGraphFileIdentity=True, fallbackStrengthPatches=fallback_count,
    maximumOrdinaryBLASLocalInverseDefect=maximum_ordinary, maximumSourceOriginalInputInverseBound=maximum_source,
    decimal80CheckedBlocks=decimal_count, maximumDecimal80LocalInverseDefect=maximum_exact, strengthCoverage=coverage,
    fixedStrengthMapDifference=map_error, singleFullFResidual=single_residual, finalFullFResidual=residual,
    previousThirtyTwoStepFullFResidual=summary['oldThirtyTwoStepFullFResidual'], fullFResidualImproved=residual < summary['oldThirtyTwoStepFullFResidual'],
    fullFInternalTargetMet=bool(summary['fullFInternalTargetMet']), fixedCParentOriginalCResidual=parent_residual,
    fixedCUnchangedAndSameRHSVerified=True, FRowsSubmatrixConsistency=identity, originalCImageDifference=parent_difference,
    oldApproximateSchurIdentityNotClaimed=True, nativeFullFArnoldiGatesPassed=True, fullFArnoldiBasisSaved=False,
    independentAllRealFullFArnoldiColumnsClaimed=False, numericalStudyFullFApplications=int(detail['CApplications'])+2,
    numericalStudyNewLocalFApplications=int(detail['localBlockApplications'])+1, numericalStudyOriginalCDiagnosticApplications=1,
    newApproximateSchurAssemblies=0, newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0,
    newProductionSolvedAreaM2=0, physicalSolution=False, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k != 'inputSHA256'}))
