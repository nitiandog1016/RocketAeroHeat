# -*- coding: utf-8 -*-
"""独立核对全部真实扩展逆、冻结核心唯一回写、完整F和固定C父矩阵诊断。"""
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
target = VAL/(name+'_independent_restricted_overlap.json')
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


def action(packed, vector):
    n, ptr, col, order, owner, values = unpack(packed)
    output = np.zeros((n, 5))
    ordered = vector[order]
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(output, owner[ii], np.einsum('kab,kb->ka', values[ii], ordered[col[ii]]))
    answer = np.empty_like(output)
    answer[order] = output
    return answer


def apply_core(factor, rhs):
    output = np.zeros_like(rhs)
    for bucket in rows(factor['buckets']):
        count = len(vec(bucket['groupIDs']))
        if not count:
            continue
        k, m = int(bucket['overlapSize']), int(bucket['coreSize'])
        members = np.asarray(bucket['members']).reshape((k, count), order='F').astype(np.int64)-1
        core = np.asarray(bucket['coreMembers']).reshape((m, count), order='F').astype(np.int64)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*k, 5*k, count), order='F').transpose(2, 0, 1)
        solved = np.einsum('pij,pj->pi', inverse, rhs[members.T].reshape(count, 5*k))
        output[core.T] = solved[:, :5*m].reshape(count, m, 5)
    return output


summary_file = folder/'restricted_overlap_summary.json'
assert summary_file.stat().st_size < 128*1024
summary = json.loads(summary_file.read_text(encoding='utf-8-sig'))
assert summary['factorAvailable'] and summary['innerRunComplete'] and summary['parentDiagnosticComplete'] and not summary['candidateRejected']
s = read(folder/'restricted_overlap_real_evidence.mat')
fdata = read(s['factorFile'])
factor = fdata['factor']
assert s['sourceBinding'] == fdata['sourceBinding'] == summary['sourceBinding']
assert sha(s['sourceBinding']['sourceFile']) == s['sourceBinding']['sourceFileSHA256']
for key, digest in s['sourceIdentity'].items():
    assert sha(ROOT/'internalFlow/cutCell'/(key+'.m')) == digest
graph = read(s['graphFile'], ['graph'])['graph']
prior = json.loads((VAL/'fullFreal01_independent_full_F.json').read_text())
assert prior['status'].startswith('passed-')
graph_key = str(Path(s['graphFile']).relative_to(ROOT))
assert bindings[graph_key] == prior['inputSHA256'][graph_key]
frozen = read(s['topologyFile'])
topology = factor['topology']
for key in ('coreGroupIDOneBased', 'coreMembersOneBased', 'corePointerOneBased', 'overlapMembersOneBased', 'overlapPointerOneBased', 'coreSizes', 'overlapSizes'):
    np.testing.assert_array_equal(vec(topology[key], np.int64), vec(frozen[key], np.int64))
n, ptr, column, order, owner, values = unpack(graph)
np.testing.assert_array_equal(order, np.arange(n))
group_ptr = vec(topology['overlapPointerOneBased'], np.int64)-1
core_ptr = vec(topology['corePointerOneBased'], np.int64)-1
all_members = vec(topology['overlapMembersOneBased'], np.int64)-1
all_core = vec(topology['coreMembersOneBased'], np.int64)-1
ng = len(core_ptr)-1
np.testing.assert_array_equal(np.sort(all_core), np.arange(n))
maximum_ordinary = maximum_exact = maximum_bound = 0.
decimal_count = fallback_count = checked = original_occurrences = 0
core_coverage = []
for bucket in rows(factor['buckets']):
    count = len(vec(bucket['groupIDs']))
    if not count:
        continue
    k, m = int(bucket['overlapSize']), int(bucket['coreSize'])
    assert 1 <= m <= 8 and m <= k <= m+2
    groups = vec(bucket['groupIDs'], np.int64)-1
    members = np.asarray(bucket['members']).reshape((k, count), order='F').astype(np.int64)-1
    core = np.asarray(bucket['coreMembers']).reshape((m, count), order='F').astype(np.int64)-1
    inverse = np.asarray(bucket['inverse']).reshape((5*k, 5*k, count), order='F').transpose(2, 0, 1)
    bounds, fallback = vec(bucket['certifiedBounds']), vec(bucket['fallback'], bool)
    assert np.all(np.isfinite(bounds) & (bounds < 1e-8))
    N = 5*k
    gamma = N*np.finfo(float).eps/(1-N*np.finfo(float).eps)
    for first in range(0, count, 128):
        last = min(count, first+128)
        block = np.zeros((last-first, N, N))
        for local, page in enumerate(range(first, last)):
            q = groups[page]
            nodes = members[:, page]
            np.testing.assert_array_equal(nodes, all_members[group_ptr[q]:group_ptr[q+1]])
            np.testing.assert_array_equal(core[:, page], all_core[core_ptr[q]:core_ptr[q+1]])
            np.testing.assert_array_equal(nodes[:m], core[:, page])
            assert len(np.unique(nodes)) == k
            core_coverage.extend(core[:, page].tolist())
            edges = np.concatenate([np.arange(ptr[node], ptr[node+1]) for node in nodes])
            rr = np.repeat(np.arange(k), ptr[nodes+1]-ptr[nodes])
            sort = np.argsort(nodes)
            sorted_nodes = nodes[sort]
            found = np.searchsorted(sorted_nodes, column[edges])
            keep = (found < k) & (sorted_nodes[np.minimum(found, k-1)] == column[edges])
            edges, rr = edges[keep], rr[keep]
            cc = sort[found[keep]]
            original_occurrences += len(edges)
            row_ids = 5*rr[:, None]+np.arange(5)[None, :]
            col_ids = 5*cc[:, None]+np.arange(5)[None, :]
            block[local, row_ids[:, :, None], col_ids[:, None, :]] = values[edges]
            adjacent = np.zeros((k, k), bool)
            adjacent[rr, cc] = np.any(values[edges] != 0, axis=(1, 2))
            adjacent |= adjacent.T
            assert all(np.any(adjacent[halo, :m]) for halo in range(m, k))
            if fallback[page]:
                diagonal = np.zeros((N, N))
                for cell in range(k):
                    ii = slice(5*cell, 5*cell+5)
                    diagonal[ii, ii] = block[local, ii, ii]
                block[local] = diagonal
                fallback_count += 1
        inv = inverse[first:last]
        ordinary = np.max(abs(block@inv-np.eye(N)), axis=(1, 2))
        rounding = gamma*np.max(np.sum(abs(block), axis=2), axis=1)*np.max(abs(inv), axis=(1, 2))/(1-gamma)**2
        maximum_ordinary = max(maximum_ordinary, float(ordinary.max()))
        maximum_bound = max(maximum_bound, float(bounds[first:last].max()))
        with localcontext() as context:
            context.prec = 80
            for local in np.flatnonzero(ordinary+rounding >= 1e-8):
                AA = [[Decimal.from_float(float(x)) for x in row] for row in block[local]]
                XX = [[Decimal.from_float(float(x)) for x in row] for row in inv[local]]
                exact = max(abs((Decimal(1) if a == b else Decimal(0))-sum((AA[a][q]*XX[q][b] for q in range(N)), Decimal(0))) for a in range(N) for b in range(N))
                assert exact < Decimal('1e-8') and exact <= Decimal.from_float(float(bounds[first+local]))
                maximum_exact = max(maximum_exact, float(exact))
                decimal_count += 1
        checked += last-first
    print('INDEPENDENT_EXTENDED_BLOCKS core='+str(m)+' extent='+str(k), flush=True)
np.testing.assert_array_equal(np.sort(core_coverage), np.arange(n))
assert checked == ng and original_occurrences == 863948
assert fallback_count == factor['factorizationReport']['fallbackOverlapBlocks']
assert fdata['measuredFactor']['bytes'] <= factor['factorizationReport']['factorStorageUpperBoundBytes'] <= 512*1024**2
assert fdata['cacheBytes'] <= 768*1024**2
del block, core_coverage
gc.collect()
input_data = read(s['inputFile'], ['rhsF', 'originalRHS'])
base = read(s['baseFile'], ['value'])
np.testing.assert_array_equal(input_data['rhsF'], s['rhsF'])
np.testing.assert_array_equal(input_data['originalRHS'], s['originalRHS'])
mapped = apply_core(factor, s['rhsF'])
mapping_error = rel(mapped, s['single'])
single_image, final_image = action(graph, s['single']), action(graph, s['valueF'])
assert mapping_error < 1e-8 and rel(single_image, s['singleImage']) < 1e-8 and rel(final_image, s['Fimage']) < 1e-8
single_residual, residual = norm(s['rhsF']-single_image)/norm(s['rhsF']), norm(s['rhsF']-final_image)/norm(s['rhsF'])
detail = s['detail']
assert abs(residual-detail['relativeFullFResidual']) < 1e-8
f = vec(graph['originalPhysicalCVs'], np.int64)-1
expected = base['value'].copy()
expected[f] = s['valueF']
np.testing.assert_array_equal(expected, s['value'])
packed = read(s['packedFile'], ['packed'])['packed']
parent_image = action(packed, s['value'])
parent_error = rel(parent_image, s['image'])
remaining = s['originalRHS']-parent_image
parent_residual = norm(remaining)/norm(s['originalRHS'])
identity = norm(remaining[f]-(s['rhsF']-final_image))/max(1, norm(s['originalRHS']))
assert parent_error < 1e-8 and identity < 1e-8 and abs(parent_residual-s['diagnostic']['relativeParentOriginalCResidual']) < 1e-8
assert not s['diagnostic']['oldApproximateSchurIdentityClaimedForNewMap']
for cycle in rows(detail['krylov']['cycles']):
    assert 'V' not in cycle and 'Z' not in cycle and cycle['basisOrthogonalityError'] < 1e-8 and cycle['arnoldiRelationRelativeDifference'] < 1e-8
assert not summary['fullFArnoldiBasisSaved'] and bool(summary['fullFInternalTargetMet']) == bool(detail['krylov']['flag'] == 0)
assert detail['FApplications'] == detail['krylov']['arnoldiSteps']+1 and detail['localFApplications'] == detail['krylov']['arnoldiSteps']
assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-all-real-extended-inverses-frozen-core-map-full-F-and-parent-diagnostic', processID=os.getpid(),
    FCVs=n, corePatches=ng, overlapMemberOccurrences=len(all_members), allCoreUnknownsReturnedExactlyOnce=True,
    topologyExactlyMatchesFrozenPreflight=True, originalFFCoefficientAuditReusedByExactGraphFileIdentity=True,
    allOriginalExtendedBlockOccurrences=original_occurrences, fallbackOverlapBlocks=fallback_count,
    maximumOrdinaryBLASInverseDefect=maximum_ordinary, maximumSourceOriginalInputInverseBound=maximum_bound,
    decimal80CheckedBlocks=decimal_count, maximumDecimal80InverseDefect=maximum_exact,
    completeRestrictedCoreMapDifference=mapping_error, singleFullFResidual=single_residual, finalFullFResidual=residual,
    previousBestFullFResidual=summary['oldThirtyTwoStepFullFResidual'], fullFResidualImproved=residual < summary['oldThirtyTwoStepFullFResidual'],
    fullFInternalTargetMet=bool(summary['fullFInternalTargetMet']), fixedCParentOriginalCResidual=parent_residual,
    fixedCUnchangedAndSameRHSVerified=True, FRowsSubmatrixConsistency=identity, originalCImageDifference=parent_error,
    oldApproximateSchurIdentityNotClaimed=True, nativeFullFArnoldiGatesPassed=True, fullFArnoldiBasisSaved=False,
    independentAllRealFullFArnoldiColumnsClaimed=False,
    actualFactorArrayBytes=int(fdata['measuredFactor']['bytes']), actualGraphAndFactorCacheBytes=int(fdata['cacheBytes']), processPeakMemoryClaimed=False,
    numericalStudyFullFApplications=int(detail['FApplications'])+2, numericalStudyRestrictedOverlapApplications=int(detail['localFApplications'])+1,
    numericalStudyOriginalCDiagnosticApplications=1, newFullFGraphExtractions=0, newApproximateSchurAssemblies=0,
    newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0,
    physicalSolution=False, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k != 'inputSHA256'}))
