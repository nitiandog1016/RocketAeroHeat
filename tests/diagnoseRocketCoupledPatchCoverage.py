# -*- coding: utf-8 -*-
"""只读新局部组与未收敛候选，量化法向/跨面覆盖及剩余分布。"""
from pathlib import Path
import os, json, hashlib, sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
VAL = ROOT / 'outputs/v83_validation'
name = sys.argv[1]
assert name.replace('_', '').isalnum()
folder = VAL / name
target = folder / 'patch_coverage_and_remaining.json'
assert not target.exists()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(8*1024**2), b''):
            h.update(part)
    return h.hexdigest()


factor_file = folder / 'coupled_patch_factor_evidence.mat'
s = loadmat(factor_file, simplify_cells=True, variable_names=['factor', 'sourceBinding'])
partition = s['factor']['partition']
n = int(partition['cellCount'])
nf = 19752
nl = 11
assert n == nf*nl
group = partition['groupID'].astype(np.int64)-1
members = partition['members'].astype(np.int64)-1
ptr = partition['groupPointer'].astype(np.int64)-1
sizes = partition['sizes'].astype(np.int64)
ng = int(partition['patchCount'])
distinct_faces = np.zeros(ng, dtype=int)
distinct_layers = np.zeros(ng, dtype=int)
one_face_partial = 0
cross_face_groups = 0
for g in range(ng):
    ids = members[ptr[g]:ptr[g+1]]
    faces = np.unique(ids % nf)
    layers = np.unique(ids // nf)
    distinct_faces[g] = len(faces)
    distinct_layers[g] = len(layers)
    one_face_partial += len(faces) == 1 and len(layers) < nl
    cross_face_groups += len(faces) > 1
line_group_counts = np.array([len(np.unique(group[f::nf])) for f in range(nf)])
full_lines_in_single_patch = int(np.count_nonzero(line_group_counts == 1))
assert full_lines_in_single_patch == 0  # 本配置每组最多8CV，完整法向列有11CV。
source = s['sourceBinding']
del s
evidence_file = folder / 'coupled_patch_real_evidence.mat'
k = loadmat(evidence_file, simplify_cells=True, variable_names=['rhs', 'image', 'detail'])
rhs = k['rhs']
remaining = rhs-k['image']
denominator = float(np.sum(remaining**2))
before_denominator = float(np.sum(rhs**2))
layer_energy = np.sum(remaining.reshape(nl, nf, 5)**2, axis=(1, 2))
equation_energy = np.sum(remaining**2, axis=0)
cv_energy = np.sum(remaining**2, axis=1)
top = np.argsort(cv_energy)[::-1][:64]
patch_energy = np.bincount(group, weights=cv_energy, minlength=ng)
assert abs(np.sqrt(denominator/before_denominator)-k['detail']['relativeCResidual']) < 1e-12
out = dict(status='read-only-real-patch-physical-members-and-unresolved-full-C-remaining-diagnostic', processID=os.getpid(),
    sourceBinding=source, cellCount=n, originalFaces=nf, layers=nl, patchCount=ng,
    patchSizeHistogram=np.bincount(sizes, minlength=9)[1:].tolist(),
    patchesSpanningMultipleOriginalFaces=int(cross_face_groups),
    patchesWithinOnePartialNormalColumn=int(one_face_partial),
    completeElevenLayerColumnsInOnePatch=full_lines_in_single_patch,
    groupsPerOriginalNormalColumnRange=[int(line_group_counts.min()), int(line_group_counts.max())],
    meanGroupsPerOriginalNormalColumn=float(line_group_counts.mean()),
    maximumDistinctFacesInPatch=int(distinct_faces.max()), maximumDistinctLayersInPatch=int(distinct_layers.max()),
    finalFullCRelativeResidual=float(np.sqrt(denominator/before_denominator)),
    finalRemainingLayerSquaredNormFractions=(layer_energy/denominator).tolist(),
    finalRemainingEquationSquaredNormFractions=(equation_energy/denominator).tolist(),
    finalRemainingFirstFourLayerSquaredNormFraction=float(layer_energy[:4].sum()/denominator),
    top64CVRemainingSquaredNormFraction=float(cv_energy[top].sum()/denominator),
    top64OriginalCVsOneBased=(top+1).tolist(),
    top64OriginalFaceIDsOneBased=(top % nf+1).tolist(), top64LayersOneBased=(top//nf+1).tolist(),
    maximumSinglePatchRemainingSquaredNormFraction=float(patch_energy.max()/denominator),
    observedResidualReductionProvesGlobalOrUniqueFailureCause=False,
    pressureOrHeatFluxErrorClaimed=False, newCApplications=0, newLocalMapApplications=0,
    newSelectedJApplications=0, newActualNSEvaluations=0, newProductionSolvedAreaM2=0,
    inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (factor_file, evidence_file)}, scriptSHA256=sha(__file__))
target.write_text(json.dumps(out, indent=2)+'\n', encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if not k.startswith('top64Original') and k != 'top64LayersOneBased'}))
