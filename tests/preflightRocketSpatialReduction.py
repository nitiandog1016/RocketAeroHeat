# -*- coding: utf-8 -*-
"""只读原C：预估按原CV分C/F的局部消元填充，不构造新数值因子或NS。"""
from pathlib import Path
import os, json, hashlib, gc
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
folder = VAL/'spatialreductionpreflight01'
assert not folder.exists()
folder.mkdir()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


factor_file = VAL/'coupledpatchreal01/coupled_patch_factor_evidence.mat'
packed_file = VAL/'blockilureal01/steady_graph_packed_evidence.mat'
assert sha(factor_file) == '09b04be0c77f76f8c04bd717745a0dd2fbcef7cdb6809d53bc403a07ed6349e0'
assert sha(packed_file) == 'beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'
s = loadmat(factor_file, simplify_cells=True, variable_names=['factor'])
part = s['factor']['partition']
del s
gc.collect()
s = loadmat(packed_file, simplify_cells=True, variable_names=['packed', 'sourceBinding'])
packed = s['packed']
binding = s['sourceBinding']
n = int(packed['cellCount'])
ng = int(part['patchCount'])
order = packed['order'].astype(np.int64)-1
ptr = packed['rowPointer'].astype(np.int64)-1
col = packed['column'].astype(np.int64)-1
owners = np.repeat(np.arange(n, dtype=np.int32), np.diff(ptr))
values = packed['values'].transpose(2, 0, 1)
group = part['groupID'].astype(np.int64)-1
sizes = part['sizes'].astype(np.int64)
members = part['members'].astype(np.int64)-1
group_ptr = part['groupPointer'].astype(np.int64)-1
block_norm = np.empty(len(col))
for first in range(0, len(col), 4096):
    block_norm[first:first+4096] = np.linalg.norm(values[first:first+4096], axis=(1, 2))
diagonal = block_norm[packed['diagonalPosition'].astype(np.int64)-1]
strength = block_norm/diagonal[owners]
strength[packed['diagonalPosition'].astype(np.int64)-1] = 0
algebraic_score = np.zeros(n)
np.maximum.at(algebraic_score, owners, strength)
score = np.empty(n)
score[order] = algebraic_score
# 每组保留相对非对角耦合最强的CV为C，其他CV仍以完整五量作为F参与原方程。
coarse_cv = np.empty(ng, dtype=np.int64)
for q in range(ng):
    ids = members[group_ptr[q]:group_ptr[q+1]]
    coarse_cv[q] = ids[np.argmax(score[ids])]
assert len(np.unique(coarse_cv)) == ng
coarse_index = np.full(n, -1, dtype=np.int64)
coarse_index[coarse_cv] = np.arange(ng)
row_cv, col_cv = order[owners], order[col]
is_c_row, is_c_col = coarse_index[row_cv] >= 0, coarse_index[col_cv] >= 0
active = block_norm != 0
cf = active & is_c_row & ~is_c_col
fc = active & ~is_c_row & is_c_col
cc = active & is_c_row & is_c_col
CF = coo_matrix((np.ones(int(cf.sum())), (coarse_index[row_cv[cf]], group[col_cv[cf]])), shape=(ng, ng)).tocsr()
FC = coo_matrix((np.ones(int(fc.sum())), (group[row_cv[fc]], coarse_index[col_cv[fc]])), shape=(ng, ng)).tocsr()
CC = coo_matrix((np.ones(int(cc.sum())), (coarse_index[row_cv[cc]], coarse_index[col_cv[cc]])), shape=(ng, ng)).tocsr()
CF.data[:] = 1
FC.data[:] = 1
CC.data[:] = 1
# 仅原CV组内的F-F块作为局部近似逆，S=D-C*M_F*B的所有可能填充都保留。
pattern = (CF@FC+CC).tocsr()
pattern.data[:] = 1
pattern.sort_indices()
f_sizes = sizes-1
inverse_numeric = int(8*25*np.sum(f_sizes**2))
schur_numeric = int(8*25*pattern.nnz)
schur_csr = int(schur_numeric+8*pattern.nnz+8*(ng+1))
pair_paths = int(np.sum(np.diff(CF.tocsc().indptr)*np.diff(FC.indptr)))
evidence_file = folder/'spatial_reduction_topology.npz'
np.savez(evidence_file, coarseOriginalCVsOneBased=coarse_cv+1, originalGroupIDOneBased=group+1,
    eliminatedCVCountsByPatch=f_sizes, schurStructuralRowPointerZeroBased=pattern.indptr,
    schurStructuralColumnZeroBased=pattern.indices)
out = dict(status='read-only-original-C-spatial-CF-local-elimination-structural-resource-preflight', processID=os.getpid(),
    sourceBinding=binding, originalCVs=n, originalScalarUnknowns=5*n,
    coarseOriginalCVs=ng, coarseScalarUnknowns=5*ng, FOriginalCVs=n-ng, FScalarUnknowns=5*(n-ng),
    representativePolicy='maximum-original-row-offdiagonal-block-Frobenius-to-diagonal-Frobenius-in-each-fixed-patch-lowest-ID-tie',
    allOriginalCVsAndFiveFieldsKeptInCompleteOperator=True, physicalGridChanged=False,
    fineBlocks=len(col), originalCFBlocks=int(cf.sum()), originalFCBlocks=int(fc.sum()), originalCCBlocks=int(cc.sum()),
    CFGroupIncidences=CF.nnz, FCGroupIncidences=FC.nnz,
    possibleLocalEliminationCoarsePairPaths=pair_paths, retainedSchurStructuralBlocks=pattern.nnz,
    localFFInverseNumericBytes=inverse_numeric, coarseSchurNumericBytes=schur_numeric,
    coarseSchurCSRArrayEstimateBytes=schur_csr, inversePlusCoarseCSRBytes=inverse_numeric+schur_csr,
    notIncludedInEstimate=['original-C-storage', 'CF-FC-coefficient-caches', 'construction-temporaries', 'deeper-levels', 'Krylov-basis', 'process-runtime'],
    processPeakMemoryClaimed=False, numericalLocalFFInversesConstructed=False,
    numericalSchurCoefficientsConstructed=False, coarseInverseOrConvergenceProven=False,
    physicalSolution=False, newCApplications=0, newSelectedJApplications=0, newActualNSEvaluations=0,
    newFlowAcceptances=0, newProductionSolvedAreaM2=0,
    topologyEvidenceFile=str(evidence_file.relative_to(ROOT)), topologyEvidenceSHA256=sha(evidence_file),
    inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (factor_file, packed_file)}, scriptSHA256=sha(__file__))
(folder/'spatial_reduction_preflight.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps(out))
