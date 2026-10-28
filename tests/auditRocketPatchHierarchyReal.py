# -*- coding: utf-8 -*-
"""独立核对真实各层全部RCP系数、缓存逆、固定V循环与最终完整C残差。"""
from pathlib import Path
from decimal import Decimal, localcontext
import os, json, sys, hashlib, gc
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
name = sys.argv[1]
assert name.replace('_', '').isalnum()
folder = VAL / name
target = VAL / (name+'_independent_patch_hierarchy.json')
assert not target.exists()
bindings = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(8*1024**2), b''):
            h.update(part)
    return h.hexdigest()


def read(path, variables=None):
    path = Path(path)
    path = path if path.is_absolute() else ROOT/path
    bindings[str(path.relative_to(ROOT))] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


def vector(value, dtype=float):
    return np.atleast_1d(value).astype(dtype)


def norm(value):
    scale = float(np.max(abs(value))) if np.size(value) else 0.
    return float(np.linalg.norm(value/scale)*scale) if scale else 0.


def rel(value, reference):
    return norm(value-reference)/max(1, norm(reference))


def unpack(packed):
    n = int(packed['cellCount'])
    ptr = vector(packed['rowPointer'], np.int64)-1
    col = vector(packed['column'], np.int64)-1
    order = vector(packed['order'], np.int64)-1
    owners = np.repeat(np.arange(n, dtype=np.int32), np.diff(ptr))
    values = np.asarray(packed['values']).reshape((5, 5, len(col)), order='F').transpose(2, 0, 1)
    return n, ptr, col, order, owners, values


def action(packed, value):
    n, ptr, col, order, owners, values = unpack(packed)
    ordered = value[order]
    output = np.zeros((n, 5))
    for first in range(0, len(col), 4096):
        ids = slice(first, first+4096)
        np.add.at(output, owners[ids], np.einsum('kab,kb->ka', values[ids], ordered[col[ids]]))
    answer = np.empty_like(output)
    answer[order] = output
    return answer


def local_map(factor, value):
    answer = np.zeros_like(value)
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(vector(bucket['patchIDs']))
        if not count:
            continue
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
        inv = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        ids = members.T
        answer[ids] = np.einsum('pij,pj->pi', inv, value[ids].reshape(count, 5*size)).reshape(count, size, 5)
    return answer


summary_file = folder/'patch_hierarchy_real_summary.json'
assert summary_file.stat().st_size < 128*1024
summary = json.loads(summary_file.read_text(encoding='utf-8-sig'))
assert summary['hierarchyAvailable'] and summary['innerRunComplete'] and not summary['candidateRejected']
hdata = read(folder/'patch_hierarchy_factor_evidence.mat')
hierarchy = hdata['hierarchy']
levels = rows(hierarchy['levels'])
binding = hdata['sourceBinding']
assert summary['sourceBinding'] == binding and sha(binding['sourceFile']) == binding['sourceFileSHA256']
for key, digest in hdata['sourceIdentity'].items():
    assert sha(ROOT/'internalFlow/cutCell'/(key+'.m')) == digest

# 根层原图及局部逆逐数组相等；既有全部真实根块逆审计据此复用。
original = read(hdata['factorFile'], ['factor'])['factor']
for a, b in zip(rows(original['buckets']), rows(levels[0]['factor']['buckets'])):
    for key in ('patchIDs', 'members', 'inverse', 'certifiedBounds', 'ordinaryDefects', 'fallback', 'rejectedPatchBounds'):
        np.testing.assert_array_equal(a[key], b[key])
for key in ('groupID', 'groupPointer', 'members', 'sizes'):
    np.testing.assert_array_equal(original['partition'][key], levels[0]['factor']['partition'][key])
del original
gc.collect()
packed_data = read(hdata['packedFile'], ['packed', 'sourceBinding'])
assert packed_data['sourceBinding'] == binding
for key in ('rowPointer', 'column', 'diagonalPosition', 'order', 'inverseOrder', 'values'):
    np.testing.assert_array_equal(packed_data['packed'][key], levels[0]['packed'][key])
del packed_data
gc.collect()
root_audit_file = VAL/'coupledpatchreal01_independent_coupled_patches.json'
root_audit = json.loads(root_audit_file.read_text(encoding='utf-8-sig'))
assert root_audit['status'].startswith('passed-')
assert root_audit['inputSHA256'][str(Path(hdata['factorFile']).relative_to(ROOT))] == bindings[str(Path(hdata['factorFile']).relative_to(ROOT))]
transfers = []
inverse_results = []

for level_id, level in enumerate(levels[:-1]):
    fine = level['packed']
    child = levels[level_id+1]['packed']
    n, ptr, col, order, owners, values = unpack(fine)
    nc, child_ptr, child_col, child_order, _, child_values = unpack(child)
    np.testing.assert_array_equal(child_order, np.arange(nc))
    part = level['factor']['partition']
    group = vector(part['groupID'], np.int64)-1
    sizes = vector(part['sizes'], np.int64)
    members = vector(part['members'], np.int64)-1
    group_ptr = vector(part['groupPointer'], np.int64)-1
    np.testing.assert_array_equal(np.sort(members), np.arange(n))
    np.testing.assert_array_equal(group[members], np.repeat(np.arange(nc), sizes))
    np.testing.assert_array_equal(np.diff(group_ptr), sizes)
    row_cv, col_cv = order[owners], order[col]
    row, column = group[row_cv], group[col_cv]
    counts = coo_matrix((np.ones(len(col)), (row, column)), shape=(nc, nc)).tocsr()
    np.testing.assert_array_equal(counts.indptr, child_ptr)
    np.testing.assert_array_equal(counts.indices, child_col)
    gamma = counts.data*np.finfo(float).eps/(1-counts.data*np.finfo(float).eps)
    maximum_coefficient_difference = 0.
    maximum_budget_ratio = 0.
    for a in range(5):
        for b in range(5):
            input_values = values[:, a, b]/sizes[row]
            coefficient = coo_matrix((input_values, (row, column)), shape=(nc, nc)).tocsr()
            sum_abs = coo_matrix((abs(input_values), (row, column)), shape=(nc, nc)).tocsr().data
            np.testing.assert_array_equal(coefficient.indptr, child_ptr)
            np.testing.assert_array_equal(coefficient.indices, child_col)
            delta = abs(coefficient.data-child_values[:, a, b])
            rounding_bound = 2*gamma*sum_abs/(1-gamma)+64*counts.data*np.finfo(float).tiny
            assert np.all(delta <= rounding_bound)
            maximum_coefficient_difference = max(maximum_coefficient_difference, float(np.max(delta/np.maximum(1, abs(coefficient.data)))))
            maximum_budget_ratio = max(maximum_budget_ratio, float(np.max(delta/rounding_bound)))
    transfers.append(dict(level=level_id+1, fineCells=n, coarseCells=nc, fineBlocks=len(col), coarseBlocks=len(child_col),
        all25CoefficientAggregationsVerified=True, maximumRelativeCoefficientDifference=maximum_coefficient_difference,
        maximumSummationRoundingBudgetRatio=maximum_budget_ratio, algebraicRPIdentityVerified=True))
    print('INDEPENDENT_ALL_RCP_COEFFICIENTS_LEVEL '+str(level_id+1), flush=True)

    if level_id == 0:
        inverse_results.append(dict(level=1, rootCacheExactlyUnchanged=True, priorAllRealRootInverseAuditReused=True))
        continue
    # 所有新局部逆按原系数批量重装，接近阈值时独立80位核对。
    position = np.empty(n, dtype=np.int32)
    position[members] = np.arange(n)-np.repeat(group_ptr[:-1], sizes)
    selected = np.flatnonzero(group[row_cv] == group[col_cv])
    selected = selected[np.argsort(group[row_cv[selected]], kind='stable')]
    edges_ptr = np.r_[0, np.cumsum(np.bincount(group[row_cv[selected]], minlength=nc))]
    maximum_ordinary = maximum_exact = maximum_source_bound = 0.
    exact_count = fallback_count = checked = 0
    for size, bucket in enumerate(rows(level['factor']['buckets']), 1):
        groups = vector(bucket['patchIDs'], np.int64)-1
        count = len(groups)
        if not count:
            continue
        inv = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        mem = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
        fallback = vector(bucket['fallback'], bool)
        bounds = vector(bucket['certifiedBounds'])
        assert np.all(np.isfinite(bounds) & (bounds < 1e-8))
        page_of_group = np.full(nc, -1, dtype=int)
        page_of_group[groups] = np.arange(count)
        N = 5*size
        g = N*np.finfo(float).eps/(1-N*np.finfo(float).eps)
        for first in range(0, count, 128):
            last = min(count, first+128)
            batch = groups[first:last]
            edges = np.concatenate([selected[edges_ptr[q]:edges_ptr[q+1]] for q in batch])
            page = page_of_group[group[row_cv[edges]]]-first
            rr, cc = position[row_cv[edges]], position[col_cv[edges]]
            block = np.zeros((last-first, N, N))
            adjacent = np.zeros((last-first, size, size), bool)
            adjacent[page, rr, cc] = np.any(values[edges] != 0, axis=(1, 2))
            adjacent |= adjacent.transpose(0, 2, 1)
            for a in range(5):
                for b in range(5):
                    block[page, 5*rr+a, 5*cc+b] = values[edges, a, b]
            for local, q in enumerate(batch):
                np.testing.assert_array_equal(mem[:, first+local], members[group_ptr[q]:group_ptr[q+1]])
                visited = {0}
                while True:
                    added = {j for j in range(size) if any(adjacent[local, j, k] for k in visited)}-visited
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
            rounding = g*np.max(np.sum(abs(block), axis=2), axis=1)*np.max(abs(inverse), axis=(1, 2))/(1-g)**2
            maximum_ordinary = max(maximum_ordinary, float(ordinary.max()))
            maximum_source_bound = max(maximum_source_bound, float(bounds[first:last].max()))
            with localcontext() as context:
                context.prec = 80
                for local in np.flatnonzero(ordinary+rounding >= 1e-8):
                    AA = [[Decimal.from_float(float(v)) for v in row_] for row_ in block[local]]
                    XX = [[Decimal.from_float(float(v)) for v in row_] for row_ in inverse[local]]
                    exact = max(abs((Decimal(1) if a == b else Decimal(0))-sum((AA[a][q]*XX[q][b] for q in range(N)), Decimal(0))) for a in range(N) for b in range(N))
                    assert exact < Decimal('1e-8') and exact <= Decimal.from_float(float(bounds[first+local]))
                    maximum_exact = max(maximum_exact, float(exact))
                    exact_count += 1
            checked += last-first
    assert checked == nc and fallback_count == level['factor']['factorizationReport']['fallbackPatchCount']
    inverse_results.append(dict(level=level_id+1, newLocalPatches=nc, allPatchesConnected=True, originalCVDiagonalFallbackPatches=fallback_count,
        maximumOrdinaryBLASInverseDefect=maximum_ordinary, maximumSourceInverseBound=maximum_source_bound,
        decimal80CheckedPatches=exact_count, maximumDecimal80InverseDefect=maximum_exact))
    print('INDEPENDENT_ALL_NEW_INVERSES_LEVEL '+str(level_id+1), flush=True)

last = levels[-1]
n, ptr, col, order, owners, values = unpack(last['packed'])
dense = np.zeros((5*n, 5*n))
for edge in range(len(col)):
    rr = order[owners[edge]]+np.arange(5)*n
    cc = order[col[edge]]+np.arange(5)*n
    dense[np.ix_(rr, cc)] = values[edge]
terminal_defect = float(np.max(abs(dense@last['inverse']-np.eye(5*n))))
assert last['terminal'] and terminal_defect < 1e-8 and rel(last['inverse'], np.linalg.solve(dense, np.eye(5*n))) < 1e-8
del row_cv, col_cv, row, column, selected, position, block
gc.collect()

s = read(folder/'patch_hierarchy_real_evidence.mat')
parent = read(hdata['rhsFile'], ['remaining', 'sourceFile'])
assert s['sourceBinding'] == binding and parent['sourceFile'] == binding['sourceFile']
rhs = s['rhs']
np.testing.assert_array_equal(rhs, parent['remaining'])
level_trace = []
root_local = None
root_local_image = None


def cycle(level_id, b):
    global root_local, root_local_image
    level = levels[level_id]
    n = int(level['packed']['cellCount'])
    if level['terminal']:
        return (level['inverse']@b.ravel(order='F')).reshape((n, 5), order='F')
    local = local_map(level['factor'], b)
    local_image = action(level['packed'], local)
    remaining = b-local_image
    part = level['factor']['partition']
    group = vector(part['groupID'], np.int64)-1
    sizes = vector(part['sizes'])
    coarse = np.zeros((int(part['patchCount']), 5))
    for equation in range(5):
        coarse[:, equation] = np.bincount(group, weights=remaining[:, equation], minlength=len(sizes))/sizes
    correction = cycle(level_id+1, coarse)[group]
    value = local+correction
    if level_id == 0:
        root_local, root_local_image = local, local_image
    level_trace.append(dict(level=level_id+1, cells=n, relativeResidualAfterLocalMap=norm(remaining)/norm(b),
        coarseRightHandSideToFineNormRatio=norm(coarse)/norm(b),
        prolongedCoarseCorrectionToLocalValueNormRatio=norm(correction)/max(np.finfo(float).tiny, norm(local)),
        VcycleValueToRightHandSideNormRatio=norm(value)/norm(b)))
    return value


single = cycle(0, rhs)
single_map_error = rel(single, s['single'])
single_image = action(levels[0]['packed'], s['single'])
single_image_error = rel(single_image, s['singleImage'])
assert single_map_error < 1e-8 and single_image_error < 1e-8
single_residual = norm(rhs-single_image)/norm(rhs)
assert abs(single_residual-s['singleSummary']['relativeCResidual'])/max(1, single_residual) < 1e-8
final_image = action(levels[0]['packed'], s['value'])
final_image_error = rel(final_image, s['image'])
residual = norm(rhs-final_image)/norm(rhs)
detail = s['detail']
assert final_image_error < 1e-8 and abs(residual-detail['relativeCResidual']) < 1e-8
kr = detail['krylov']
steps = int(kr['arnoldiSteps'])
assert detail['solverCApplications'] == steps+1 and detail['VcycleApplications'] == steps
np.testing.assert_array_equal(detail['preconditionerCApplicationsByLevel'], [steps]*(len(levels)-1)+[0])
np.testing.assert_array_equal(detail['preconditionerLocalApplicationsByLevel'], [steps]*(len(levels)-1)+[0])
np.testing.assert_array_equal(detail['preconditionerTerminalApplicationsByLevel'], [0]*(len(levels)-1)+[steps])
assert detail['totalRootCApplications'] == 2*steps+1
assert not summary['realFullArnoldiBasisSaved']
for entry in rows(kr['cycles']):
    assert 'V' not in entry and 'Z' not in entry and entry['basisOrthogonalityError'] < 1e-8 and entry['arnoldiRelationRelativeDifference'] < 1e-8
assert bool(summary['internalCTargetMet']) == bool(kr['flag'] == 0)

# 只在保存同右端的局部/粗校正二维张成空间内诊断，不把它称为新迭代或原NS方向。
images = np.column_stack([root_local_image.ravel(order='F'), (single_image-root_local_image).ravel(order='F')])
scales = np.linalg.norm(images, axis=0)
normalized = images/scales[None, :]
coefficients = np.linalg.lstsq(normalized, rhs.ravel(order='F'), rcond=1e-13)[0]
pair_residual = norm(rhs.ravel(order='F')-normalized@coefficients)/norm(rhs)
pair_weights = coefficients/scales
assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-real-all-RCP-coefficients-new-inverses-fixed-V-cycle-and-final-full-C', processID=os.getpid(),
    levelCellCounts=[int(e['packed']['cellCount']) for e in levels], transfers=transfers, inverseChecks=inverse_results,
    rootPackedGraphAndLocalCacheExactlyUnchanged=True, priorRootAllInverseAuditReusedByExactArrayIdentity=True,
    terminalFullFiveFieldInverseDefect=terminal_defect, fixedVcycleMapRelativeDifference=single_map_error,
    singleFullCImageRelativeDifference=single_image_error, singleFullCResidual=single_residual,
    finalFullCImageRelativeDifference=final_image_error, finalFullCResidual=residual,
    previousLocalOnlyFullCResidual=root_audit['finalFullCResidual'], fullCResidualImprovedComparedWithLocalOnly=residual < root_audit['finalFullCResidual'],
    internalCTargetMet=bool(summary['internalCTargetMet']), nativeArnoldiGatesPassed=True,
    realFullArnoldiBasisSaved=False, independentAllRealArnoldiColumnAuditClaimed=False,
    singleVcycleTrace=sorted(level_trace, key=lambda row:row['level']),
    sameRhsTwoImageDiagnostic=dict(relativeFullCResidual=pair_residual, localAndCoarseScalarWeights=pair_weights.tolist(),
        diagnosticOnly=True, implementedAsNewSolver=False, extraActualNSEvaluations=0),
    solverRootCApplications=int(detail['solverCApplications']), preconditionerRootCApplications=steps,
    totalSolverPlusPreconditionerRootCApplications=int(detail['totalRootCApplications']),
    constructionFineAndCoarseCProbeApplications=2*(len(levels)-1), singleProbePreconditionerCApplications=len(levels)-1,
    extraRootCDiagnosticApplications=2, independentAuditFixedVcycleApplications=1,
    newSelectedJApplications=0, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    originalNSConvergenceOrRegionalQualificationProven=False, physicalSolution=False,
    inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k != 'inputSHA256'}))
