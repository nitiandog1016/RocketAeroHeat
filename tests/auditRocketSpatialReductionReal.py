# -*- coding: utf-8 -*-
"""独立核对真实原子块、全部F/粗局部逆、近似Schur系数及完整原C残差。"""
from pathlib import Path
from decimal import Decimal, localcontext
import os, json, hashlib, sys, gc
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
    return n, ptr, col, order, owner, values


def block_action(packed, value):
    n, ptr, col, order, owner, values = unpack(packed)
    ordered = value[order]
    out = np.zeros((n, 5))
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(out, owner[ii], np.einsum('kab,kb->ka', values[ii], ordered[col[ii]]))
    answer = np.empty_like(out)
    answer[order] = out
    return answer


def rect_action(coupling, value):
    owner, column = vec(coupling['owner'], np.int64)-1, vec(coupling['column'], np.int64)-1
    values = np.asarray(coupling['values']).reshape((5, 5, len(column)), order='F').transpose(2, 0, 1)
    out = np.zeros((int(coupling['rowCount']), 5))
    for first in range(0, len(column), 4096):
        ii = slice(first, first+4096)
        np.add.at(out, owner[ii], np.einsum('kab,kb->ka', values[ii], value[column[ii]]))
    return out


def local_map(factor, value, group_key):
    out = np.zeros_like(value)
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(vec(bucket[group_key]))
        if not count:
            continue
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        ids = members.T
        out[ids] = np.einsum('pij,pj->pi', inverse, value[ids].reshape(count, 5*size)).reshape(count, size, 5)
    return out


def audit_local(factor, values, row_nodes, col_nodes, group, counts, position, group_key, label):
    ng = len(counts)
    selected = np.flatnonzero((row_nodes >= 0) & (col_nodes >= 0))
    selected = selected[group[row_nodes[selected]] == group[col_nodes[selected]]]
    selected = selected[np.argsort(group[row_nodes[selected]], kind='stable')]
    edge_ptr = np.r_[0, np.cumsum(np.bincount(group[row_nodes[selected]], minlength=ng))]
    maximum_ordinary = maximum_exact = maximum_source_bound = 0.
    exact_count = fallback_count = checked = 0
    covered = []
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        groups = vec(bucket[group_key], np.int64)-1
        count = len(groups)
        if not count:
            continue
        inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
        bounds, fallback = vec(bucket['certifiedBounds']), vec(bucket['fallback'], bool)
        assert np.all(np.isfinite(bounds) & (bounds < 1e-8))
        page_of_group = np.full(ng, -1, dtype=int)
        page_of_group[groups] = np.arange(count)
        N = 5*size
        gamma = N*np.finfo(float).eps/(1-N*np.finfo(float).eps)
        for first in range(0, count, 128):
            last = min(count, first+128)
            batch = groups[first:last]
            edges = np.concatenate([selected[edge_ptr[q]:edge_ptr[q+1]] for q in batch])
            page = page_of_group[group[row_nodes[edges]]]-first
            rr, cc = position[row_nodes[edges]], position[col_nodes[edges]]
            block = np.zeros((last-first, N, N))
            for a in range(5):
                for b in range(5):
                    block[page, 5*rr+a, 5*cc+b] = values[edges, a, b]
            for local, g in enumerate(batch):
                ids = members[:, first+local]
                np.testing.assert_array_equal(group[ids], np.full(size, g))
                np.testing.assert_array_equal(position[ids], np.arange(size))
                covered.extend(ids.tolist())
                if fallback[first+local]:
                    diagonal = np.zeros((N, N))
                    for cell in range(size):
                        ii = slice(5*cell, 5*cell+5)
                        diagonal[ii, ii] = block[local, ii, ii]
                    block[local] = diagonal
                    fallback_count += 1
            inv = inverse[first:last]
            ordinary = np.max(abs(block@inv-np.eye(N)), axis=(1, 2))
            rounding = gamma*np.max(np.sum(abs(block), axis=2), axis=1)*np.max(abs(inv), axis=(1, 2))/(1-gamma)**2
            maximum_ordinary = max(maximum_ordinary, float(ordinary.max()))
            maximum_source_bound = max(maximum_source_bound, float(bounds[first:last].max()))
            with localcontext() as context:
                context.prec = 80
                for local in np.flatnonzero(ordinary+rounding >= 1e-8):
                    AA = [[Decimal.from_float(float(v)) for v in row] for row in block[local]]
                    XX = [[Decimal.from_float(float(v)) for v in row] for row in inv[local]]
                    exact = max(abs((Decimal(1) if a == b else Decimal(0))-sum((AA[a][q]*XX[q][b] for q in range(N)), Decimal(0))) for a in range(N) for b in range(N))
                    assert exact < Decimal('1e-8') and exact <= Decimal.from_float(float(bounds[first+local]))
                    maximum_exact = max(maximum_exact, float(exact))
                    exact_count += 1
            checked += last-first
    np.testing.assert_array_equal(np.sort(covered), np.arange(int(factor['cellCount'])))
    assert checked == np.count_nonzero(counts)
    print('INDEPENDENT_ALL_LOCAL_INVERSES '+label, flush=True)
    return dict(localBlocks=checked, allUnknownsExactlyOnce=True, fallbackBlocks=fallback_count,
        maximumOrdinaryBLASDefect=maximum_ordinary, maximumSourceOriginalInputBound=maximum_source_bound,
        decimal80CheckedBlocks=exact_count, maximumDecimal80InverseDefect=maximum_exact)


summary_file = folder/'spatial_reduction_real_summary.json'
assert summary_file.stat().st_size < 128*1024
summary = json.loads(summary_file.read_text(encoding='utf-8-sig'))
assert summary['assemblyComplete'] and summary['schurProbesComplete'] and summary['coarseFactorAvailable'] and summary['innerRunComplete'] and not summary['candidateRejected']
a = read(folder/'spatial_reduction_assembly_evidence.mat')
reduction, binding = a['reduction'], a['sourceBinding']
assert summary['sourceBinding'] == binding and sha(binding['sourceFile']) == binding['sourceFileSHA256']
for key, digest in a['sourceIdentity'].items():
    assert sha(ROOT/'internalFlow/cutCell'/(key+'.m')) == digest
root = read(a['packedFile'], ['packed', 'sourceBinding'])
assert root['sourceBinding'] == binding
packed = root['packed']
n, ptr, col, order, owner, values = unpack(packed)
split = reduction['split']
f, c = vec(split['fineOriginalCVs'], np.int64)-1, vec(split['coarseOriginalCVs'], np.int64)-1
np.testing.assert_array_equal(np.sort(np.r_[f, c]), np.arange(n))
input_data = read(a['inputFile'])
np.testing.assert_array_equal(c, vec(input_data['coarseOriginalCVsOneBased'], np.int64)-1)
nf, nc = len(f), len(c)
fine_index, coarse_index = vec(split['fineIndex'], np.int64)-1, vec(split['coarseIndex'], np.int64)-1
row_cv, col_cv = order[owner], order[col]
row_f, col_f, row_c, col_c = fine_index[row_cv], fine_index[col_cv], coarse_index[row_cv], coarse_index[col_cv]
for key, mask, rr, cc in [('B', (row_f >= 0)&(col_c >= 0), row_f, col_c), ('E', (row_c >= 0)&(col_f >= 0), row_c, col_f), ('D', (row_c >= 0)&(col_c >= 0), row_c, col_c)]:
    coupling = reduction[key]
    np.testing.assert_array_equal(vec(coupling['owner'], np.int64)-1, rr[mask])
    np.testing.assert_array_equal(vec(coupling['column'], np.int64)-1, cc[mask])
    np.testing.assert_array_equal(np.asarray(coupling['values']).transpose(2, 0, 1), values[mask])
group = vec(split['fineGroupID'], np.int64)-1
counts = vec(split['groupFCounts'], np.int64)
position = vec(split['finePosition'], np.int64)-1
f_checks = audit_local(reduction['Ffactor'], values, row_f, col_f, group, counts, position, 'groupIDs', 'F')

# 重装每个局部消元贡献，并用两次矩阵乘法及累加的舍入界核对全部S系数。
schur = reduction['schur']
ns, sptr, scol, sorder, sowner, svalues = unpack(schur)
assert ns == nc
np.testing.assert_array_equal(sorder, np.arange(nc))
np.testing.assert_array_equal(sptr, vec(input_data['schurStructuralRowPointerZeroBased'], np.int64))
np.testing.assert_array_equal(scol, vec(input_data['schurStructuralColumnZeroBased'], np.int64))
keys = sowner.astype(np.int64)*nc+scol
expected = np.zeros_like(svalues)
absolute_sum = np.zeros_like(svalues)
product_error = np.zeros_like(svalues)
contribution_counts = np.zeros(len(scol))
D = reduction['D']
drow, dcol = vec(D['owner'], np.int64)-1, vec(D['column'], np.int64)-1
dvalues = np.asarray(D['values']).transpose(2, 0, 1)
where = np.searchsorted(keys, drow*nc+dcol)
np.testing.assert_array_equal(keys[where], drow*nc+dcol)
expected[where] = dvalues
absolute_sum[where] = abs(dvalues)
contribution_counts[where] = 1
E, B = reduction['E'], reduction['B']
erow, ecol = vec(E['owner'], np.int64)-1, vec(E['column'], np.int64)-1
brow, bcol = vec(B['owner'], np.int64)-1, vec(B['column'], np.int64)-1
evals, bvals = np.asarray(E['values']).transpose(2, 0, 1), np.asarray(B['values']).transpose(2, 0, 1)
incoming, outgoing = group[ecol], group[brow]
ei, bi = np.argsort(incoming, kind='stable'), np.argsort(outgoing, kind='stable')
ep = np.r_[0, np.cumsum(np.bincount(incoming, minlength=nc))]
bp = np.r_[0, np.cumsum(np.bincount(outgoing, minlength=nc))]
buckets = rows(reduction['Ffactor']['buckets'])
pages = vec(reduction['Ffactor']['pageOfGroup'], np.int64)-1
for g in range(nc):
    size = int(counts[g])
    if not size:
        continue
    ie, oe = ei[ep[g]:ep[g+1]], bi[bp[g]:bp[g+1]]
    if not len(ie) or not len(oe):
        continue
    r, rp = np.unique(erow[ie], return_inverse=True)
    q, qp = np.unique(bcol[oe], return_inverse=True)
    left, right = np.zeros((5*len(r), 5*size)), np.zeros((5*size, 5*len(q)))
    for aa in range(5):
        for bb in range(5):
            left[5*rp+aa, 5*position[ecol[ie]]+bb] = evals[ie, aa, bb]
            right[5*position[brow[oe]]+aa, 5*qp+bb] = bvals[oe, aa, bb]
    inv = np.asarray(buckets[size-1]['inverse']).reshape((5*size, 5*size, len(vec(buckets[size-1]['groupIDs']))), order='F')[:, :, pages[g]]
    update = -(left@inv)@right
    gamma = 5*size*np.finfo(float).eps/(1-5*size*np.finfo(float).eps)
    magnitude = (abs(left)@abs(inv))@abs(right)/(1-gamma)**2
    error = 2*(2*gamma+gamma**2)*magnitude
    rr, cc = np.repeat(r, len(q)), np.tile(q, len(r))
    positions = np.searchsorted(keys, rr*nc+cc)
    np.testing.assert_array_equal(keys[positions], rr*nc+cc)
    def blocks(matrix):
        return matrix.reshape(len(r), 5, len(q), 5).transpose(0, 2, 1, 3).reshape(-1, 5, 5)
    expected[positions] += blocks(update)
    absolute_sum[positions] += blocks(magnitude)
    product_error[positions] += blocks(error)
    contribution_counts[positions] += 1
    if (g+1) % 8000 == 0:
        print('INDEPENDENT_SCHUR_LOCAL_CONTRIBUTIONS '+str(g+1), flush=True)
gamma = contribution_counts*np.finfo(float).eps/(1-contribution_counts*np.finfo(float).eps)
error_bound = product_error/(1-gamma[:, None, None])+2*gamma[:, None, None]*absolute_sum/(1-gamma[:, None, None])**2+64*np.finfo(float).tiny
delta = abs(expected-svalues)
assert np.all(delta <= error_bound)
maximum_schur_budget_ratio = float(np.max(delta/error_bound))
maximum_schur_coefficient_difference = float(np.max(delta/np.maximum(1, abs(expected))))
del expected, absolute_sum, product_error, error_bound, delta, contribution_counts, values, row_cv, col_cv, row_f, col_f, row_c, col_c
gc.collect()
print('INDEPENDENT_ALL_APPROXIMATE_SCHUR_COEFFICIENTS_VERIFIED', flush=True)
cfdata = read(folder/'spatial_coarse_factor_evidence.mat')
coarse = cfdata['coarseFactor']
assert cfdata['sourceBinding'] == binding
part = coarse['partition']
cgroup, csizes = vec(part['groupID'], np.int64)-1, vec(part['sizes'], np.int64)
cmembers = vec(part['members'], np.int64)-1
cpointer = vec(part['groupPointer'], np.int64)-1
cposition = np.empty(nc, dtype=int)
cposition[cmembers] = np.arange(nc)-np.repeat(cpointer[:-1], csizes)
coarse_checks = audit_local(coarse, svalues, sorder[sowner], sorder[scol], cgroup, csizes, cposition, 'patchIDs', 'coarse-S')
probes = read(folder/'spatial_schur_probe_evidence.mat', ['schurProbes'])['schurProbes']
probe_results = []
for probe in rows(probes):
    right = probe['rhs']
    implicit = rect_action(D, right)-rect_action(E, local_map(reduction['Ffactor'], rect_action(B, right), 'groupIDs'))
    image = block_action(schur, right)
    implicit_error, image_error = rel(implicit, probe['implicit']), rel(image, probe['assembled'])
    assert implicit_error < 1e-8 and image_error < 1e-8 and rel(image, implicit) < 1e-8
    probe_results.append(dict(implicitMapDifference=implicit_error, assembledSImageDifference=image_error, actualSActionRelation=rel(image, implicit)))
s = read(folder/'spatial_reduction_real_evidence.mat')
parent = read(a['rhsFile'], ['remaining', 'sourceFile'])
assert s['sourceBinding'] == binding and parent['sourceFile'] == binding['sourceFile']
rhs = s['rhs']
np.testing.assert_array_equal(rhs, parent['remaining'])
bF, bC = rhs[f], rhs[c]
first = local_map(reduction['Ffactor'], bF, 'groupIDs')
coarse_rhs = bC-rect_action(E, first)
xC = local_map(coarse, coarse_rhs, 'patchIDs')
xF = local_map(reduction['Ffactor'], bF-rect_action(B, xC), 'groupIDs')
single = np.zeros_like(rhs)
single[f], single[c] = xF, xC
map_error = rel(single, s['single'])
single_image, final_image = block_action(packed, s['single']), block_action(packed, s['value'])
single_image_error, final_image_error = rel(single_image, s['singleImage']), rel(final_image, s['image'])
assert max(map_error, single_image_error, final_image_error) < 1e-8
single_residual, residual = norm(rhs-single_image)/norm(rhs), norm(rhs-final_image)/norm(rhs)
detail = s['detail']
assert abs(single_residual-s['singleSummary']['relativeCResidual'])/max(1, single_residual) < 1e-8 and abs(residual-detail['relativeCResidual']) < 1e-8
steps = int(detail['krylov']['arnoldiSteps'])
assert detail['CApplications'] == steps+1 and detail['FLocalApplications'] == 2*steps and detail['coarseLocalApplications'] == steps and detail['rectangularCouplingApplications'] == 2*steps
for cycle in rows(detail['krylov']['cycles']):
    assert 'V' not in cycle and 'Z' not in cycle and cycle['basisOrthogonalityError'] < 1e-8 and cycle['arnoldiRelationRelativeDifference'] < 1e-8
assert not summary['realFullArnoldiBasisSaved'] and bool(summary['internalCTargetMet']) == bool(detail['krylov']['flag'] == 0)
assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-real-spatial-split-all-F-and-coarse-inverses-Schur-coefficients-and-complete-C', processID=os.getpid(),
    originalCVs=n, FOriginalCVs=nf, coarseOriginalCVs=nc, allOriginalCVsAndFiveFieldsKept=True,
    allOriginalBCEDCoefficientsExactlyVerified=True, FLocalInverseChecks=f_checks, coarseLocalInverseChecks=coarse_checks,
    approximateSchurBlocks=len(scol), allApproximateSchurCoefficientsVerifiedWithinProductAndSumRoundingBudget=True,
    maximumSchurCoefficientNormalizedDifference=maximum_schur_coefficient_difference, maximumSchurRoundingBudgetRatio=maximum_schur_budget_ratio,
    coarsePatternAndRepresentativesMatchFrozenStructuralPreflight=True, threeImplicitSchurProbes=probe_results,
    completeSpatialBackSubstitutionMapDifference=map_error, singleOriginalCImageDifference=single_image_error,
    finalOriginalCImageDifference=final_image_error, singleOriginalCResidual=single_residual, finalOriginalCResidual=residual,
    internalCTargetMet=bool(summary['internalCTargetMet']), nativeArnoldiGatesPassed=True, realFullArnoldiBasisSaved=False,
    independentAllRealArnoldiColumnAuditClaimed=False, solverCApplications=steps+1, spatialMapApplications=steps,
    solverFLocalApplications=2*steps, solverCoarseLocalApplications=steps, solverRectangularApplications=2*steps,
    numericalStudyOriginalCApplications=steps+3, numericalStudyFLocalApplications=2*steps+5,
    numericalStudyCoarseLocalApplications=steps+1, numericalStudyRectangularApplications=2*steps+11,
    numericalStudyPackedSchurProbeApplications=3, fullFFInverseOrExactSchurClaimed=False,
    newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0,
    originalNSConvergenceOrRegionalQualificationProven=False, physicalSolution=False, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k != 'inputSHA256'}))
