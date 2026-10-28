# -*- coding: utf-8 -*-
"""独立核对真实全部F系数、冻结局部逆复用、完整F残差和固定C原矩阵诊断。"""
from pathlib import Path
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
target = VAL/(name+'_independent_full_F.json')
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


def action(packed, value):
    n, col, order, owner, values = unpack(packed)
    output = np.zeros((n, 5))
    ordered = value[order]
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(output, owner[ii], np.einsum('kab,kb->ka', values[ii], ordered[col[ii]]))
    answer = np.empty_like(output)
    answer[order] = output
    return answer


def rect(coupling, value):
    owner, col = vec(coupling['owner'], np.int64)-1, vec(coupling['column'], np.int64)-1
    values = np.asarray(coupling['values']).reshape((5, 5, len(col)), order='F').transpose(2, 0, 1)
    output = np.zeros((int(coupling['rowCount']), 5))
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(output, owner[ii], np.einsum('kab,kb->ka', values[ii], value[col[ii]]))
    return output


def local(factor, value):
    output = np.zeros_like(value)
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(vec(bucket['groupIDs']))
        if not count:
            continue
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        ids = members.T
        output[ids] = np.einsum('pij,pj->pi', inverse, value[ids].reshape(count, 5*size)).reshape(count, size, 5)
    return output


summary_file = folder/'full_F_real_summary.json'
assert summary_file.stat().st_size < 128*1024
summary = json.loads(summary_file.read_text(encoding='utf-8-sig'))
assert summary['innerRunComplete'] and summary['parentDiagnosticComplete'] and not summary['candidateRejected']
s = read(folder/'full_F_real_evidence.mat')
g = read(s['graphFile'])
graph = g['graph']
binding = s['sourceBinding']
assert g['sourceBinding'] == summary['sourceBinding'] == binding
for key, digest in s['sourceIdentity'].items():
    assert sha(ROOT/'internalFlow/cutCell'/(key+'.m')) == digest
old = json.loads((VAL/'spatialreductionreal01_independent_spatial_reduction.json').read_text())
oldS = json.loads((VAL/'spatialschurreal01_independent_spatial_schur.json').read_text())
assert old['status'].startswith('passed-') and oldS['status'].startswith('passed-')
assembly = read(s['assemblyFile'], ['reduction'])
r = assembly['reduction']
factor, split, B = r['Ffactor'], r['split'], r['B']
del r, assembly
gc.collect()
assembly_key = str(Path(s['assemblyFile']).relative_to(ROOT))
assert bindings[assembly_key] == old['inputSHA256'][assembly_key]
packed = read(s['packedFile'], ['packed'])['packed']
n, col, order, owner, values = unpack(packed)
nf, gfcol, gorder, gfowner, gfvalues = unpack(graph)
f, c = vec(split['fineOriginalCVs'], np.int64)-1, vec(split['coarseOriginalCVs'], np.int64)-1
np.testing.assert_array_equal(vec(graph['originalPhysicalCVs'], np.int64)-1, f)
np.testing.assert_array_equal(gorder, np.arange(nf))
assert nf == len(f) == 187749 and n == 217272
np.testing.assert_array_equal(np.sort(np.r_[f, c]), np.arange(n))
fine_index = vec(split['fineIndex'], np.int64)-1
row_cv, col_cv = order[owner], order[col]
row_f, col_f = fine_index[row_cv], fine_index[col_cv]
selected = np.flatnonzero((row_f >= 0)&(col_f >= 0))
selected = selected[np.lexsort((col_f[selected], row_f[selected]))]
np.testing.assert_array_equal(gfowner, row_f[selected])
np.testing.assert_array_equal(gfcol, col_f[selected])
for first in range(0, len(selected), 32768):
    ii = slice(first, first+32768)
    np.testing.assert_array_equal(gfvalues[ii], values[selected[ii]])
group = vec(split['originalGroupID'], np.int64)-1
cross_count = int(np.count_nonzero(group[row_cv[selected]] != group[col_cv[selected]]))
assert len(selected) == 1285801 and cross_count == 682056
del selected, row_cv, col_cv, row_f, col_f, values, owner
gc.collect()
source = read(s['sourceCandidateFile'], ['value', 'rhs', 'diagnostic', 'sourceBinding'])
source_key = str(Path(s['sourceCandidateFile']).relative_to(ROOT))
assert bindings[source_key] == oldS['inputSHA256'][source_key]
np.testing.assert_array_equal(source['rhs'], s['originalRHS'])
np.testing.assert_array_equal(source['value'][c], s['fixedC'])
initial = read(folder/'full_F_inputs.mat')
rhsF = source['rhs'][f]-rect(B, s['fixedC'])
np.testing.assert_array_equal(source['value'][f], initial['previousFValue'])
assert rel(rhsF, s['rhsF']) < 1e-8
prior_map = local(factor, rhsF)
assert rel(prior_map, initial['previousFValue']) < 1e-8
previous_image = action(graph, initial['previousFValue'])
final_image = action(graph, s['valueF'])
assert rel(previous_image, initial['previousFImage']) < 1e-8 and rel(final_image, s['Fimage']) < 1e-8
previous_residual = norm(rhsF-previous_image)/norm(rhsF)
residual = norm(rhsF-final_image)/norm(rhsF)
d = s['detail']
assert abs(previous_residual-initial['previousFullFResidual']) < 1e-8 and abs(residual-d['relativeFullFResidual']) < 1e-8
full = np.zeros_like(s['originalRHS'])
full[f], full[c] = s['valueF'], s['fixedC']
np.testing.assert_array_equal(full, s['value'])
parent_image = action(packed, full)
parent_error = rel(parent_image, s['image'])
assert parent_error < 1e-8
remaining = s['originalRHS']-parent_image
parent_residual = norm(remaining)/norm(s['originalRHS'])
f_rows, c_rows = norm(remaining[f])/norm(s['originalRHS']), norm(remaining[c])/norm(s['originalRHS'])
identity = norm(remaining[f]-(rhsF-final_image))/max(1, norm(s['originalRHS']))
diagnostic = s['diagnostic']
assert abs(parent_residual-diagnostic['relativeParentOriginalCResidual']) < 1e-8 and identity < 1e-8
assert not diagnostic['oldApproximateSchurIdentityClaimedForNewFValue']
assert not summary['fullFArnoldiBasisSaved']
for cycle in rows(d['krylov']['cycles']):
    assert 'V' not in cycle and 'Z' not in cycle and cycle['basisOrthogonalityError'] < 1e-8 and cycle['arnoldiRelationRelativeDifference'] < 1e-8
assert bool(summary['fullFInternalTargetMet']) == bool(d['krylov']['flag'] == 0)
assert d['FApplications'] == d['krylov']['arnoldiSteps']+1 and d['localFApplications'] == d['krylov']['arnoldiSteps']
assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-all-real-FF-coefficients-frozen-local-map-final-F-and-fixed-C-parent-diagnostic', processID=os.getpid(),
    originalCVs=n, FCVs=nf, originalFFBlocks=1285801, crossPatchFFBlocks=cross_count,
    allFFCoefficientsExactlyMatchOriginalC=True, sourceOriginalPhysicalCVMapExactlyVerified=True,
    existingLocalFInverseAuditReusedByAssemblyByteIdentity=True, newLocalFInverseConstructions=0,
    sameFixedCAndSameFRhsAsPreviousCandidate=True, previousFullFResidual=previous_residual, finalFullFResidual=residual,
    fullFInternalTargetMet=bool(summary['fullFInternalTargetMet']), parentOriginalCResidual=parent_residual,
    FRowsRelativeToOriginalRHS=f_rows, CRowsRelativeToOriginalRHS=c_rows, FRowsSubmatrixConsistency=identity,
    oldApproximateSchurIdentityNotClaimedForNewFValue=True, parentOriginalCTargetMet=bool(diagnostic['parentOriginalCTargetMet']),
    originalCImageDifference=parent_error, nativeFullFArnoldiGatesPassed=True,
    fullFArnoldiBasisSaved=False, independentAllRealFullFArnoldiColumnsClaimed=False,
    fullFSolverApplications=int(d['FApplications']), fullFDiagnosticApplications=2, fixedLocalFApplications=int(d['localFApplications']),
    originalCDiagnosticApplications=1, newApproximateSchurAssemblies=0, newSelectedJApplications=0, newActualNSEvaluations=0,
    newFlowAcceptances=0, newProductionSolvedAreaM2=0, physicalSolution=False, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k != 'inputSHA256'}))
