# -*- coding: utf-8 -*-
"""只读完整父C的选定局部像列，诊断标量归一后的相关性和局部QR可行性。"""
from pathlib import Path
import os,json
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'): os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha
from screenRocketParentImageScale import setup

target=VAL/'parent_local_image_diagnostic01.json';assert not target.exists();bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


s=read(VAL/'coupledpatchreal01/coupled_patch_factor_evidence.mat',['factor','packedFile']);factor=s['factor']
packed=read(s['packedFile'],['packed'])['packed'];A,plan,order=setup(packed,factor)
norms=read(VAL/'parentCimagescalescreen02/parent_image_scale_evidence.mat',['columnNormAlgebraic'])['columnNormAlgebraic']
groups=[]
entries=factor['buckets'];entries=[entries] if isinstance(entries,dict) else list(np.atleast_1d(entries).ravel())
for b in entries:
    ids=np.atleast_1d(b['patchIDs']).astype(np.int64)
    if ids.size: groups.append(ids)
group_ids=np.concatenate(groups);nodes=[];inverses=[]
for members,inv in plan:
    nodes.extend(list(members));inverses.extend(list(inv))
assert len(nodes)==len(group_ids)==29523
scores=np.array([norms[g].max() for g in nodes]);rank=np.argsort(-scores,kind='stable')
selected=np.unique(np.r_[rank[:16],np.linspace(0,len(nodes)-1,16).astype(np.int64)])
owner=np.repeat(np.arange(len(order)),np.diff(A.indptr));col=A.indices;sorting=np.argsort(col,kind='stable')
counts=np.bincount(col,minlength=len(order));ptr=np.r_[0,np.cumsum(counts)];results=[]
for index in selected:
    g=nodes[index];inv=inverses[index]
    blocks=np.concatenate([sorting[ptr[j]:ptr[j+1]] for j in g]);row,rr=np.unique(owner[blocks],return_inverse=True)
    cc=np.repeat(np.arange(len(g)),counts[g]);panel=np.zeros((5*len(row),5*len(g)))
    ri=5*rr[:,None]+np.arange(5)[None,:];ci=5*cc[:,None]+np.arange(5)[None,:]
    panel[ri[:,:,None],ci[:,None,:]]=A.data[blocks]
    image=panel@inv;normalized=image/np.linalg.norm(image,axis=0)[None,:]
    gram=normalized.T@normalized;np.fill_diagonal(gram,0)
    singular=np.linalg.svd(normalized,compute_uv=False)
    Q,R=np.linalg.qr(image,mode='reduced');right_map=np.linalg.solve(R.T,inv.T).T
    actual_image=panel@right_map
    orth=float(np.linalg.norm(actual_image.T@actual_image-np.eye(actual_image.shape[1])))
    agreement=float(np.linalg.norm(actual_image-Q)/max(1.,np.linalg.norm(Q)))
    results.append(dict(originalPatchIDOneBased=int(group_ids[index]),originalCVsOneBased=(order[g]+1).tolist(),
        completeImageRows=len(row),columnNormMaximum=float(norms[g].max()),
        scalarNormalizedMaximumAbsoluteColumnCorrelation=float(np.max(abs(gram))),
        scalarNormalizedLocalSingularRange=[float(singular.min()),float(singular.max())],
        scalarNormalizedLocalCondition=float(singular.max()/singular.min()),
        proposedQRRightMapCompleteImageOrthogonality=orth,proposedQRRightMapAgreement=agreement,
        localQRNumericalGateMet=orth<1e-8 and agreement<1e-8))
assert all(sha(ROOT/path)==digest for path,digest in bindings.items())
out=dict(status='read-only-selected-complete-parent-local-image-correlation-and-QR-feasibility-diagnostic',processID=os.getpid(),
    selectedGroups=len(results),selection='16-largest-complete-image-column-norm-groups-plus-16-uniform-factor-indices',
    allGroupsOrGlobalSpectrumClaimed=False,fullGlobalQRPreconditionerConstructed=False,
    newProductionQualificationProven=False,localGroupResults=results,
    maximumSelectedNormalizedLocalCondition=max(x['scalarNormalizedLocalCondition'] for x in results),
    maximumSelectedColumnCorrelation=max(x['scalarNormalizedMaximumAbsoluteColumnCorrelation'] for x in results),
    allSelectedQRNumericalGatesMet=all(x['localQRNumericalGateMet'] for x in results),
    newCompleteCApplications=0,separateReadOnlyLocalImageProducts=2*len(results),
    newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('inputSHA256','localGroupResults')}))
