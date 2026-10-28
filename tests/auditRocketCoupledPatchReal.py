# -*- coding: utf-8 -*-
"""独立核对真实五量局部块覆盖、原输入逆证书、完整C及固定局部作用。"""
from pathlib import Path
from decimal import Decimal, localcontext
import os, json, hashlib, sys, gc
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
name = sys.argv[1]
assert name.replace('_', '').isalnum()
folder = VAL / name
target = VAL / (name + '_independent_coupled_patches.json')
assert not target.exists()
bindings = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path, variables=None):
    path = Path(path)
    path = path if path.is_absolute() else ROOT / path
    bindings[str(path.relative_to(ROOT))] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


def norm(value):
    scale = float(np.max(abs(value))) if np.size(value) else 0.
    return float(np.linalg.norm(value / scale) * scale) if scale else 0.


def rel(value, reference):
    return norm(value - reference) / max(1, norm(reference))


def rows(value):
    return [value] if isinstance(value, dict) else list(np.atleast_1d(value).ravel())


summary_file = folder / 'coupled_patch_real_summary.json'
assert summary_file.stat().st_size < 128 * 1024
summary = json.loads(summary_file.read_text(encoding='utf-8-sig'))
assert summary['factorAvailable']
fdata = read(folder / 'coupled_patch_factor_evidence.mat')
factor = fdata['factor']
binding = fdata['sourceBinding']
assert binding == summary['sourceBinding']
for key, digest in fdata['sourceIdentity'].items():
    assert sha(ROOT / 'internalFlow/cutCell' / (key + '.m')) == digest
assert sha(binding['sourceFile']) == binding['sourceFileSHA256']
data = read(fdata['packedFile'], ['packed', 'sourceBinding'])
packed = data['packed']
assert data['sourceBinding'] == binding
assert bindings[str(Path(fdata['packedFile']).relative_to(ROOT))] == 'beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'
n = int(packed['cellCount'])
order = packed['order'].astype(np.int64) - 1
mapping = packed['inverseOrder'].astype(np.int64) - 1
np.testing.assert_array_equal(order[mapping], np.arange(n))
ptr = packed['rowPointer'].astype(np.int64) - 1
col = packed['column'].astype(np.int64) - 1
owners = np.repeat(np.arange(n, dtype=np.int32), np.diff(ptr))
values = packed['values'].transpose(2, 0, 1)
partition = factor['partition']
group = partition['groupID'].astype(np.int64) - 1
members = partition['members'].astype(np.int64) - 1
group_ptr = partition['groupPointer'].astype(np.int64) - 1
sizes = partition['sizes'].astype(np.int64)
ng = int(partition['patchCount'])
np.testing.assert_array_equal(np.sort(members), np.arange(n))
np.testing.assert_array_equal(np.diff(group_ptr), sizes)
np.testing.assert_array_equal(group[members], np.repeat(np.arange(ng), sizes))
assert np.all((sizes >= 1) & (sizes <= 8)) and n == factor['cellCount']
position = np.empty(n, dtype=np.int32)
position[members] = np.arange(n) - np.repeat(group_ptr[:-1], sizes)
row_cv = order[owners]
col_cv = order[col]
selected = np.flatnonzero(group[row_cv] == group[col_cv])
selected = selected[np.argsort(group[row_cv[selected]], kind='stable')]
edge_ptr = np.r_[0, np.cumsum(np.bincount(group[row_cv[selected]], minlength=ng))]
assert len(selected) == factor['factorizationReport']['includedOriginalBlocksBeforeFallback']
maximum_ordinary = 0.
maximum_exact = 0.
maximum_bound = 0.
decimal_count = 0
fallback_count = 0
connectivity_checked = 0
all_groups = []

# 原块按128个局部组逐批装配，不申请全局稠密矩阵或复制全部缓存逆。
for size, bucket in enumerate(rows(factor['buckets']), 1):
    groups = np.atleast_1d(bucket['patchIDs']).astype(np.int64) - 1
    count = len(groups)
    if not count:
        continue
    all_groups.extend(groups.tolist())
    mem = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64) - 1
    inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
    source_bounds = np.atleast_1d(bucket['certifiedBounds'])
    fallback = np.atleast_1d(bucket['fallback']).astype(bool)
    assert np.all(np.isfinite(source_bounds) & (source_bounds < 1e-8))
    page_of_group = np.full(ng, -1, dtype=np.int32)
    page_of_group[groups] = np.arange(count)
    N = 5 * size
    gamma = N * np.finfo(float).eps / (1 - N * np.finfo(float).eps)
    for first in range(0, count, 128):
        last = min(count, first + 128)
        group_batch = groups[first:last]
        edges = np.concatenate([selected[edge_ptr[g]:edge_ptr[g+1]] for g in group_batch])
        page = page_of_group[group[row_cv[edges]]] - first
        rr = position[row_cv[edges]]
        cc = position[col_cv[edges]]
        block = np.zeros((last-first, N, N))
        adjacent = np.zeros((last-first, size, size), dtype=bool)
        adjacent[page, rr, cc] = np.any(values[edges] != 0, axis=(1, 2))
        adjacent |= adjacent.transpose(0, 2, 1)
        for a in range(5):
            for b in range(5):
                block[page, 5*rr+a, 5*cc+b] = values[edges, a, b]
        for local, g in enumerate(group_batch):
            np.testing.assert_array_equal(mem[:, first+local], members[group_ptr[g]:group_ptr[g+1]])
            visited = {0}
            while True:
                added = {j for j in range(size) if any(adjacent[local, j, q] for q in visited)} - visited
                if not added:
                    break
                visited |= added
            assert len(visited) == size
            connectivity_checked += 1
            if fallback[first+local]:
                diagonal = np.zeros((N, N))
                for cell in range(size):
                    ids = slice(5*cell, 5*cell+5)
                    diagonal[ids, ids] = block[local, ids, ids]
                block[local] = diagonal
                fallback_count += 1
        inv = inverse[first:last]
        ordinary = np.max(abs(block @ inv - np.eye(N)), axis=(1, 2))
        rounding = gamma * np.max(np.sum(abs(block), axis=2), axis=1) * np.max(abs(inv), axis=(1, 2)) / (1-gamma)**2
        maximum_ordinary = max(maximum_ordinary, float(ordinary.max()))
        maximum_bound = max(maximum_bound, float(source_bounds[first:last].max()))
        # 接近阈值的块以原double输入的80位乘积独立核验；普通乘积与证书保持不同身份。
        needs_exact = np.flatnonzero(ordinary + rounding >= 1e-8)
        assert np.all(ordinary[np.setdiff1d(np.arange(last-first), needs_exact)] < 1e-8)
        with localcontext() as context:
            context.prec = 80
            for local in needs_exact:
                AA = [[Decimal.from_float(float(v)) for v in row] for row in block[local]]
                XX = [[Decimal.from_float(float(v)) for v in row] for row in inv[local]]
                exact = max(abs((Decimal(1) if i == j else Decimal(0)) - sum((AA[i][q]*XX[q][j] for q in range(N)), Decimal(0))) for i in range(N) for j in range(N))
                assert exact < Decimal('1e-8') and exact <= Decimal.from_float(float(source_bounds[first+local]))
                maximum_exact = max(maximum_exact, float(exact))
                decimal_count += 1
                if decimal_count % 100 == 0:
                    print('INDEPENDENT_DECIMAL80_COUPLED_PATCHES ' + str(decimal_count), flush=True)
    print('INDEPENDENT_ALL_FIVE_FIELD_PATCH_BUCKET ' + str(size), flush=True)
np.testing.assert_array_equal(np.sort(all_groups), np.arange(ng))
assert connectivity_checked == ng and fallback_count == factor['factorizationReport']['fallbackPatchCount']
del row_cv, col_cv, selected, edge_ptr, position, block
gc.collect()


def action(vector):
    ordered = vector[order]
    output = np.zeros((n, 5))
    for first in range(0, len(col), 4096):
        ids = slice(first, first+4096)
        np.add.at(output, owners[ids], np.einsum('kab,kb->ka', values[ids], ordered[col[ids]]))
    result = np.empty_like(output)
    result[order] = output
    return result


def precondition(vector):
    output = np.zeros((n, 5))
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(np.atleast_1d(bucket['patchIDs']))
        if not count:
            continue
        mem = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64) - 1
        inv = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        physical = mem.T
        b = vector[physical].reshape(count, 5*size)
        output[physical] = np.einsum('pij,pj->pi', inv, b).reshape(count, size, 5)
    return output


assert summary['innerRunComplete'] and not summary['candidateRejected'], 'A rejected Arnoldi run requires a separate saved-basis audit.'
s = read(folder / 'coupled_patch_real_evidence.mat')
assert s['sourceBinding'] == binding
parent = read(fdata['rhsFile'], ['remaining', 'sourceFile'])
assert parent['sourceFile'] == binding['sourceFile']
rhs = s['rhs']
np.testing.assert_array_equal(rhs, parent['remaining'])
map_error = rel(precondition(rhs), s['single'])
single_image = action(s['single'])
image_error = rel(single_image, s['singleImage'])
assert map_error < 1e-8 and image_error < 1e-8
single_residual = norm(rhs - single_image) / norm(rhs)
assert abs(single_residual-s['singleSummary']['relativeCResidual'])/max(1, single_residual) < 1e-8
final_image = action(s['value'])
final_image_error = rel(final_image, s['image'])
residual = norm(rhs-final_image) / norm(rhs)
detail = s['detail']
assert final_image_error < 1e-8 and abs(residual-detail['relativeCResidual'])/max(1, residual) < 1e-8
krylov = detail['krylov']
cycle = rows(krylov['cycles'])[0]
assert cycle['arnoldiRelationRelativeDifference'] < 1e-8 and cycle['basisOrthogonalityError'] < 1e-8
assert detail['CApplications'] == krylov['modelApplications'] == krylov['arnoldiSteps']+1
assert detail['localBlockApplications'] == krylov['preconditionerApplications'] == krylov['arnoldiSteps']
assert not summary['realFullArnoldiBasisSaved'] and 'V' not in cycle and 'Z' not in cycle
assert bool(summary['internalCTargetMet']) == bool(krylov['flag'] == 0)
assert summary['extraFullCDiagnosticApplications'] == 2
for key in ('newSelectedJApplications', 'newActualNSEvaluations', 'newFlowAcceptances', 'newProductionSolvedAreaM2'):
    assert summary[key] == 0
assert not summary['physicalSolution']
assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-real-all-coupled-patch-inverses-complete-C-and-final-residual', processID=os.getpid(),
    cellCount=n, scalarUnknowns=5*n, patchCount=ng, maximumPatchCells=8, allOriginalCVsExactlyOnce=True,
    allPatchesConnected=True, originalCVDiagonalFallbackPatches=fallback_count,
    maximumOrdinaryBLASInverseDefect=maximum_ordinary, maximumCertifiedOriginalInputInverseDefect=maximum_bound,
    decimal80CheckedPatches=decimal_count, maximumDecimal80InverseDefect=maximum_exact,
    fixedFullPreconditionerMapRelativeDifference=map_error, singleFullCImageRelativeDifference=image_error,
    singleFullCResidual=single_residual, finalFullCImageRelativeDifference=final_image_error, finalFullCResidual=residual,
    internalCTargetMet=bool(summary['internalCTargetMet']), nativeArnoldiGatesPassed=True,
    realFullArnoldiBasisSaved=False, independentAllRealArnoldiColumnAuditClaimed=False,
    internalSolverCApplications=int(detail['CApplications']), internalLocalMapApplications=int(detail['localBlockApplications']),
    extraFullCDiagnosticApplications=2, independentAuditFullCApplications=2, independentAuditLocalMapApplications=1,
    newSelectedJApplications=0, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    originalNSConvergenceOrNewRegionalQualificationProven=False, physicalSolution=False,
    inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
