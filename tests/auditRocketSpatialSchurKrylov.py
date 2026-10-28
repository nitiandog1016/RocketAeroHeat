# -*- coding: utf-8 -*-
"""独立重放全部粗S方向、块回代及原F行误差来源，S和原C范数严格分开。"""
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
target = VAL/(name+'_independent_spatial_schur.json')
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
    out = np.zeros((n, 5))
    ordered = value[order]
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(out, owner[ii], np.einsum('kab,kb->ka', values[ii], ordered[col[ii]]))
    answer = np.empty_like(out)
    answer[order] = out
    return answer


def rect(coupling, value):
    owner, col = vec(coupling['owner'], np.int64)-1, vec(coupling['column'], np.int64)-1
    values = np.asarray(coupling['values']).reshape((5, 5, len(col)), order='F').transpose(2, 0, 1)
    out = np.zeros((int(coupling['rowCount']), 5))
    for first in range(0, len(col), 4096):
        ii = slice(first, first+4096)
        np.add.at(out, owner[ii], np.einsum('kab,kb->ka', values[ii], value[col[ii]]))
    return out


def local(factor, value, key):
    answer = np.zeros_like(value)
    for size, bucket in enumerate(rows(factor['buckets']), 1):
        count = len(vec(bucket[key]))
        if not count:
            continue
        members = np.asarray(bucket['members']).reshape((size, count), order='F').astype(np.int64)-1
        inverse = np.asarray(bucket['inverse']).reshape((5*size, 5*size, count), order='F').transpose(2, 0, 1)
        ids = members.T
        answer[ids] = np.einsum('pij,pj->pi', inverse, value[ids].reshape(count, 5*size)).reshape(count, size, 5)
    return answer


summary = json.loads((folder/'spatial_schur_krylov_summary.json').read_text(encoding='utf-8-sig'))
assert summary['numericalRunComplete'] and not summary['candidateRejected'] and summary['fullSArnoldiBasisSaved']
s = read(folder/'spatial_schur_krylov_evidence.mat')
for key, digest in s['sourceIdentity'].items():
    assert sha(ROOT/'internalFlow/cutCell'/(key+'.m')) == digest
audit = json.loads((VAL/'spatialreductionreal01_independent_spatial_reduction.json').read_text())
assert audit['status'].startswith('passed-')
a = read(s['assemblyFile'], ['reduction', 'sourceBinding'])
coarse = read(s['coarseFactorFile'], ['coarseFactor'])['coarseFactor']
for path in (s['assemblyFile'], s['coarseFactorFile']):
    key = str(Path(path).relative_to(ROOT))
    assert bindings[key] == audit['inputSHA256'][key]
reduction = a['reduction']
split = reduction['split']
f, c = vec(split['fineOriginalCVs'], np.int64)-1, vec(split['coarseOriginalCVs'], np.int64)-1
rhs = s['rhs']
initial = read(folder/'spatial_schur_condensed_rhs.mat')
np.testing.assert_array_equal(rhs, initial['rhs'])
first = local(reduction['Ffactor'], rhs[f], 'groupIDs')
coarse_rhs = rhs[c]-rect(reduction['E'], first)
assert rel(first, initial['first']) < 1e-8 and rel(coarse_rhs, s['coarseRhs']) < 1e-8
local_coarse = local(coarse, coarse_rhs, 'patchIDs')
single_image = action(reduction['schur'], local_coarse)
single_residual = norm(coarse_rhs-single_image)/norm(coarse_rhs)
assert rel(local_coarse, initial['singleCoarse']) < 1e-8 and rel(single_image, initial['singleSImage']) < 1e-8
assert abs(single_residual-initial['singleSResidual']) < 1e-8
d = s['detail']
cycle = rows(d['krylov']['cycles'])[0]
V, Z, H = cycle['V'], cycle['Z'], cycle['smallH']
nc = len(c)
images = np.zeros_like(Z)
maximum_mapping = 0.
for j in range(Z.shape[1]):
    images[:, j] = action(reduction['schur'], Z[:, j].reshape((nc, 5), order='F')).ravel(order='F')
    expected = local(coarse, V[:, j].reshape((nc, 5), order='F'), 'patchIDs').ravel(order='F')/cycle['imageNormalizationScales'][j]
    maximum_mapping = max(maximum_mapping, rel(expected, Z[:, j]))
    if (j+1) % 8 == 0:
        print('INDEPENDENT_ALL_APPROXIMATE_S_ARNOLDI_COLUMNS '+str(j+1), flush=True)
relation = rel(images, V@H)
orth = float(np.linalg.norm(V.T@V-np.eye(V.shape[1])))
assert maximum_mapping < 1e-8 and relation < 1e-8 and orth < 1e-8
weights = cycle['weights']
candidate = cycle['xStart']+Z@weights
assert rel(candidate, s['coarseValue'].ravel(order='F')) < 1e-8
Simage = action(reduction['schur'], s['coarseValue'])
Sremaining = coarse_rhs-Simage
Sresidual = norm(Sremaining)/norm(coarse_rhs)
predicted = cycle['rStart']-V@(H@weights)
combination = norm(Sremaining.ravel(order='F')-predicted)/max(1, norm(coarse_rhs))
assert rel(Simage, s['Simage']) < 1e-8 and abs(Sresidual-d['relativeCResidual']) < 1e-8 and combination < 1e-8
value = np.zeros_like(rhs)
value[c] = s['coarseValue']
right_F = rhs[f]-rect(reduction['B'], value[c])
value[f] = local(reduction['Ffactor'], right_F, 'groupIDs')
assert rel(value, s['value']) < 1e-8
packed = read(s['packedFile'], ['packed'])['packed']
image = action(packed, s['value'])
assert rel(image, s['image']) < 1e-8
remaining = rhs-image
Cresidual = norm(remaining)/norm(rhs)
fine_residual, coarse_residual = norm(remaining[f])/norm(rhs), norm(remaining[c])/norm(rhs)
identity = norm(remaining[c]-Sremaining)/max(1, norm(rhs))
fraction = norm(remaining[f])**2/norm(remaining)**2
diagnostic = s['diagnostic']
assert abs(Cresidual-diagnostic['relativeFullCResidual']) < 1e-8 and abs(fraction-diagnostic['fineRowsRemainingSquaredNormFraction']) < 1e-8 and identity < 1e-8

# 直接使用原F-F系数分解当前失败向量：局部逆缺陷与跨组F-F贡献分别测量。
n, col, order, owner, values = unpack(packed)
row_cv, col_cv = order[owner], order[col]
fine_index = vec(split['fineIndex'], np.int64)-1
row_F, col_F = fine_index[row_cv], fine_index[col_cv]
groups = vec(split['originalGroupID'], np.int64)-1
selected = np.flatnonzero((row_F >= 0)&(col_F >= 0)&(groups[row_cv] == groups[col_cv]))
local_A_image = np.zeros((len(f), 5))
for first in range(0, len(selected), 4096):
    ids = selected[first:first+4096]
    np.add.at(local_A_image, row_F[ids], np.einsum('kab,kb->ka', values[ids], value[f][col_F[ids]]))
full_A_image = image[f]-rect(reduction['B'], value[c])
cross_FF = full_A_image-local_A_image
local_inverse_defect = right_F-local_A_image
decomposition = norm(remaining[f]-(local_inverse_defect-cross_FF))/max(1, norm(rhs))
assert decomposition < 1e-8
assert all(sha(ROOT/path) == digest for path, digest in bindings.items())
out = dict(status='passed-independent-all-spatial-S-physical-directions-condensation-and-parent-F-row-decomposition', processID=os.getpid(),
    assemblyAndFactorPriorFullCoefficientAndInverseAuditReusedByByteIdentity=True, coarseScalarUnknowns=5*nc,
    allPhysicalApproximateSArnoldiColumns=Z.shape[1], allColumnSImageRelation=relation,
    maximumFixedCoarseMapDifference=maximum_mapping, basisOrthogonalityError=orth, finalSCombinationRelation=combination,
    originalFullCDiagnosticImageDifference=rel(image, s['image']), singleApproximateSResidual=single_residual,
    finalApproximateSResidual=Sresidual, finalOriginalCResidual=Cresidual,
    FRowsRelativeToOriginalRHS=fine_residual, CRowsRelativeToOriginalRHS=coarse_residual,
    coarseRowsCondensationConsistency=identity, FRowsRemainingSquaredNormFraction=fraction,
    localFInverseEquationDefectRelativeToOriginalRHS=norm(local_inverse_defect)/norm(rhs),
    crossPatchFFImageRelativeToOriginalRHS=norm(cross_FF)/norm(rhs), FRowExactDecompositionDifference=decomposition,
    attributionScope='this-saved-same-RHS-vector-only-not-a-global-rank-or-universal-cause-proof',
    approximateSTargetMet=bool(diagnostic['approximateSTargetMet']), originalCInternalTargetMet=bool(diagnostic['originalCInternalTargetMet']),
    actualSAndOriginalCResidualDenominatorsKeptDistinct=True, newFInverseConstructions=0, newSchurAssemblies=0,
    numericalStudyPackedSApplications=int(d['CApplications'])+2, numericalStudyOriginalCDiagnosticApplications=1,
    numericalStudyCoarseLocalApplications=int(d['localBlockApplications'])+1, numericalStudyFLocalApplications=2,
    numericalStudyRectangularApplications=2, independentAuditPackedSApplications=Z.shape[1]+2,
    newSelectedJApplications=0, newActualNSEvaluations=0, newFlowAcceptances=0, newProductionSolvedAreaM2=0,
    physicalSolution=False, inputSHA256=bindings, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k != 'inputSHA256'}))
