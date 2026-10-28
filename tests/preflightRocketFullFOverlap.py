# -*- coding: utf-8 -*-
"""只读两条F候选与完整F强连接，冻结每个核心组最多两个邻点的重叠预检。"""
from pathlib import Path
import os, json, hashlib, gc
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat, savemat
from scipy.sparse import coo_matrix
ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT/'outputs/v83_validation'
folder = VAL/'fullFoverlappreflight01'
assert not folder.exists()
folder.mkdir()
bindings = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path, variables=None):
    path = Path(path)
    bindings[str(path.relative_to(ROOT))] = sha(path)
    return loadmat(path, simplify_cells=True, variable_names=variables)


e1 = VAL/'fullFreal01/full_F_real_evidence.mat'
e2 = VAL/'fullFstrengthreal02/full_F_strength_real_evidence.mat'
a = read(e1, ['rhsF', 'valueF', 'Fimage'])
b = read(e2, ['rhsF', 'valueF', 'Fimage'])
np.testing.assert_array_equal(a['rhsF'], b['rhsF'])
rhs = a['rhsF']
r1, r2 = rhs-a['Fimage'], rhs-b['Fimage']
correlation = float(np.sum(r1*r2)/(np.linalg.norm(r1)*np.linalg.norm(r2)))
images = np.column_stack([a['Fimage'].ravel(order='F'), b['Fimage'].ravel(order='F')])
scales = np.linalg.norm(images, axis=0)
normalized = images/scales[None, :]
weights = np.linalg.lstsq(normalized, rhs.ravel(order='F'), rcond=1e-13)[0]
predicted = float(np.linalg.norm(rhs.ravel(order='F')-normalized@weights)/np.linalg.norm(rhs))
pair_weights = (weights/scales).tolist()
del a, b, r2, images, normalized
gc.collect()
factor_file = VAL/'fullFstrengthreal02/full_F_strength_factor_evidence.mat'
fdata = read(factor_file, ['factor'])
part = fdata['factor']['partition']
del fdata
gc.collect()
graph_file = VAL/'fullFreal01/full_F_graph_evidence.mat'
g = read(graph_file, ['graph'])['graph']
n = int(g['cellCount'])
group = part['groupID'].astype(np.int64)-1
core = part['members'].astype(np.int64)-1
core_ptr = part['groupPointer'].astype(np.int64)-1
core_sizes = part['sizes'].astype(np.int64)
ng = int(part['patchCount'])
ptr, col, order = g['rowPointer'].astype(np.int64)-1, g['column'].astype(np.int64)-1, g['order'].astype(np.int64)-1
np.testing.assert_array_equal(order, np.arange(n))
owners = np.repeat(np.arange(n, dtype=np.int32), np.diff(ptr))
values = g['values'].transpose(2, 0, 1)
magnitude = np.empty(len(col))
for first in range(0, len(col), 4096):
    magnitude[first:first+4096] = np.linalg.norm(values[first:first+4096], axis=(1, 2))
diagonal = magnitude[g['diagonalPosition'].astype(np.int64)-1]
weights = magnitude/diagonal[owners]
weights[g['diagonalPosition'].astype(np.int64)-1] = 0
G = coo_matrix((weights, (owners, col)), shape=(n, n)).tocsr()
G = G.maximum(G.T)
G.setdiag(0)
G.eliminate_zeros()
original_ids = g['originalPhysicalCVs'].astype(np.int64)
extensions = []
ext_counts = np.zeros(ng, dtype=np.int64)
member_lists = []
selected_edges = 0
core_to_extended_edges = 0
for q in range(ng):
    nodes = core[core_ptr[q]:core_ptr[q+1]]
    pools, strengths = [], []
    for node in nodes:
        ii = slice(G.indptr[node], G.indptr[node+1])
        candidate, strength = G.indices[ii], G.data[ii]
        outside = group[candidate] != q
        pools.extend(candidate[outside].tolist())
        strengths.extend(strength[outside].tolist())
    if pools:
        candidate, inv = np.unique(pools, return_inverse=True)
        score = np.zeros(len(candidate))
        np.maximum.at(score, inv, strengths)
        rank = np.lexsort((original_ids[candidate], -score))
        extra = candidate[rank[:2]]
    else:
        extra = np.empty(0, dtype=np.int64)
    members = np.r_[nodes, extra]
    assert len(np.unique(members)) == len(members) and len(extra) <= 2
    member_lists.append(members)
    extensions.extend(extra.tolist())
    ext_counts[q] = len(extra)
    for node in members:
        edges = col[ptr[node]:ptr[node+1]]
        included = np.isin(edges, members, assume_unique=True)
        selected_edges += int(np.count_nonzero(included))
        if group[node] == q:
            core_to_extended_edges += int(np.count_nonzero(included))
member_ptr = np.r_[0, np.cumsum(core_sizes+ext_counts)]
members = np.concatenate(member_lists)
membership_count = np.bincount(members, minlength=n)
sizes = core_sizes+ext_counts
numeric_inverse = int(8*25*np.sum(sizes**2))
factor_bound = numeric_inverse+int(96*n+160*ng+8*len(members)+65536)
np.testing.assert_array_equal(np.sort(core), np.arange(n))
nf = 19752
nl = 11
assert g['sourceOriginalCellCount'] == nf*nl
layer = (original_ids-1)//nf
layer_energy = np.bincount(layer, weights=np.sum(r1**2, axis=1), minlength=nl)
equation_energy = np.sum(r1**2, axis=0)
denominator = float(np.sum(r1**2))
top = np.argsort(np.sum(r1**2, axis=1))[::-1][:64]
input_file = folder/'full_F_overlap_input.mat'
savemat(input_file, dict(coreGroupIDOneBased=group+1, coreMembersOneBased=core+1, corePointerOneBased=core_ptr+1,
    overlapMembersOneBased=members+1, overlapPointerOneBased=member_ptr+1, coreSizes=core_sizes,
    overlapSizes=sizes, graphEvidenceSHA256=sha(graph_file), sourceStrengthFactorSHA256=sha(factor_file)), do_compression=True)
out = dict(status='read-only-fixed-core-two-neighbor-restricted-overlap-topology-and-partial-storage-preflight', processID=os.getpid(),
    FCVs=n, corePatches=ng, originalFGraphBlocks=len(col), maximumAddedFNeighbors=2,
    coreAndFCRolesExactlyUnchanged=True, coreUnknownsReturnedExactlyOnce=True, overlapPhysicalCellsMerged=False,
    overlapMembers=int(len(members)), addedMembers=int(ext_counts.sum()), maximumCVOverlapMembership=int(membership_count.max()),
    minimumOverlapCells=int(sizes.min()), maximumOverlapCells=int(sizes.max()),
    originalFFBlockOccurrencesInsideOverlapSets=selected_edges, coreRowsToOverlapOriginalFFBlocks=core_to_extended_edges,
    inverseNumericArrayBytes=numeric_inverse, conservativeFactorArrayBoundBytes=factor_bound,
    proposedFactorBudgetBytes=512*1024**2, factorBoundFitsProposedBudget=factor_bound <= 512*1024**2,
    constructionAndWholeRuntimeBudgetProven=False, numericalOverlapInversesConstructed=False, convergenceOrRegionalQualificationProven=False,
    boundarySelection='two-highest-symmetric-row-diagonal-normalized-five-field-block-neighbors-original-global-CV-ID-tie',
    outputPolicy='solve-extended-local-system-return-only-original-core-rows-fixed-map-prototype-not-yet-implemented',
    savedTwoResidualCorrelation=correlation, savedTwoCandidateImageSpanPredictedFResidual=predicted,
    savedTwoCandidateWeights=pair_weights, imageSpanResultIsDiagnosticOnly=True,
    betterSavedFullFRemainingLayerSquaredNormFractions=(layer_energy/denominator).tolist(),
    betterSavedFullFRemainingEquationSquaredNormFractions=(equation_energy/denominator).tolist(),
    top64FUnknownRemainingSquaredNormFraction=float(np.sum(r1[top]**2)/denominator),
    topologyMATFile=str(input_file.relative_to(ROOT)), topologyMATFileSHA256=sha(input_file),
    newFullFApplications=0, newOriginalCApplications=0, newLocalInverseConstructions=0,
    newActualNSEvaluations=0, newProductionSolvedAreaM2=0, physicalSolution=False,
    inputSHA256=bindings, scriptSHA256=sha(__file__))
(folder/'full_F_overlap_preflight.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k != 'inputSHA256'}))
