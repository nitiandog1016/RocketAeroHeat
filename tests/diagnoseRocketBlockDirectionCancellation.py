# -*- coding: utf-8 -*-
"""逐个读取已拒绝的方向，区分小H病态、组合相消与实际C剩余；不调用模型。"""
from pathlib import Path
import os,json,hashlib
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'block_direction_cancellation01.json';assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def positive_ratio(a,b):return float(a/max(np.finfo(float).tiny,b))
results=[];bindings={}
for case,name in [('blockilufailurediag01','arnoldi_failure_evidence.mat'),('scaledilureal01','scaled_arnoldi_failure_evidence.mat'),('blockfillreal01','scaled_arnoldi_failure_evidence.mat')]:
 file=VAL/case/name;f=loadmat(file,simplify_cells=True)['failure'];Z=f['Z'];V=f['V'];H=f['smallH'];w=np.atleast_1d(f['weights']);s=np.linalg.svd(H,compute_uv=False)
 combined=Z@w;image=V@(H@w);z_bound2=v_bound2=0.
 for first in range(0,len(combined),4096):
  ids=slice(first,first+4096);z_bound2+=float(np.sum((abs(Z[ids])@abs(w))**2));v_bound2+=float(np.sum((abs(V[ids])@(abs(H)@abs(w)))**2))
 active=V.shape[1]-int(np.linalg.norm(V[:,-1])==0);normb=np.linalg.norm(f['b'])
 results.append(dict(case=case,arnoldiSteps=int(f['arnoldiSteps']),HSmallestSingularValue=float(s.min()),HLargestSingularValue=float(s.max()),
  HConditionEstimate=positive_ratio(s.max(),s.min()),HRetainedRankAtUnchangedRelativeThreshold=int(np.count_nonzero(s>1e-13*s.max())),
  coefficientNormRelativeToRhs=positive_ratio(np.linalg.norm(w),normb),
  physicalDirectionCombinationCancellation=positive_ratio(np.sqrt(z_bound2),np.linalg.norm(combined)),
  projectedImageCombinationCancellation=positive_ratio(np.sqrt(v_bound2),np.linalg.norm(image)),
  candidateRecomputationRelativeToFormation=positive_ratio(np.linalg.norm(f['candidate']-f['xStart']-combined),np.sqrt(z_bound2)),
  actualRelativeCResidual=positive_ratio(np.linalg.norm(f['actualRemaining']),normb),
  savedOrthogonalityError=float(f['orthogonalityError']),savedCombinedRelationDifference=float(f['modelRelationRelativeDifference']),
  basisActiveColumns=active,usableDirection=False))
 bindings[str(file.relative_to(ROOT))]=sha(file);del f,Z,V,H,w,combined,image
out=dict(status='diagnosed-saved-rejected-direction-conditioning-and-cancellation-no-new-model',processID=os.getpid(),cases=results,inputSHA256=bindings,
 newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
