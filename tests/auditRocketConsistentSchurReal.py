# -*- coding: utf-8 -*-
"""独立核对真实一致Schur五条作用、六次块回代及全部局部对角弱耦合。"""
from pathlib import Path
import os,json,hashlib,sys,gc
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1];assert name.replace('_','').isalnum()
target=VAL/(name+'_independent_consistent_schur.json');assert not target.exists();bindings={}
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(path,variables=None):
 path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[str(path.relative_to(ROOT))]=sha(path)
 return loadmat(path,simplify_cells=True,variable_names=variables)
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def stable_norm(a):
 scale=float(np.max(abs(a))) if np.size(a) else 0.
 return float(np.linalg.norm(a/scale)*scale) if scale else 0.
def rel(a,b):return float(stable_norm(a-b)/max(1,stable_norm(b)))
file=VAL/name/'consistent_schur_real_evidence.mat';s=read(file);report=s['report'];binding=s['sourceBinding'];n=int(report['cellCount']);order=s['algebraicOrder'].astype(int)-1
refinement_verified=False
if 'refinementReport' in s:
 refinement=s['refinementReport'];previous=Path(refinement['sourceFile']);previous=previous if previous.is_absolute() else ROOT/previous
 assert sha(previous)==refinement['sourceSHA256'];assert sha(ROOT/'internalFlow/cutCell/rocketTwoFieldSingularRatio.m')==refinement['helperSHA256']
 assert refinement['newMomentumSolves']==refinement['newActualNSEvaluations']==refinement['newProductionSolvedAreaM2']==0
 old=read(previous)
 def unchanged(a,b,path=''):
  if path in ('sigmaRatio','report.localDiagnostics'):return
  if isinstance(a,dict):
   for key in a:
    if key.startswith('__'):continue
    assert key in b;unchanged(a[key],b[key],path+'.'+key if path else key)
  elif isinstance(a,str):assert a==b
  else:np.testing.assert_array_equal(a,b)
 unchanged(old,s);refinement_verified=True;del old
for key,value in s['sourceIdentity'].items():assert sha(ROOT/'internalFlow/cutCell'/(key+'.m'))==value
p=read(binding['packedFile'],['packed'])['packed'];assert bindings[str(Path(binding['packedFile']).relative_to(ROOT))]==binding['packedFileSHA256']
np.testing.assert_array_equal(p['order'].astype(int)-1,order)
ptr=p['rowPointer'].astype(np.int64)-1;col=p['column'].astype(np.int64)-1;cv=p['values'].transpose(2,0,1);owners=np.repeat(np.arange(n,dtype=np.int32),np.diff(ptr));diag=p['diagonalPosition'].astype(np.int64)-1
mf=read(s['momentumFile'],['M','sourceBinding']);assert mf['sourceBinding']==binding;mf=mf['M'];mv=mf['values'].transpose(2,0,1);mi=mf['inverseDiagonal'].transpose(2,0,1)
np.testing.assert_array_equal(mf['rowPointer'].astype(np.int64)-1,ptr);np.testing.assert_array_equal(mf['column'].astype(np.int64)-1,col)
md=mf['diagonalPosition'].astype(np.int64)-1
scale=read(s['scaleFile'],['rowScale','columnScale','assemblySHA256']);R=scale['rowScale'];C=scale['columnScale']
schur=read(s['assemblyFile'],['schurMatrix'])['schurMatrix'];assert bindings[str(Path(s['assemblyFile']).relative_to(ROOT))]==scale['assemblySHA256']
sptr=schur['rowPointer'].astype(np.int64)-1;scol=schur['column'].astype(np.int64)-1;sv=schur['values'].transpose(2,0,1);so=np.repeat(np.arange(n,dtype=np.int32),np.diff(sptr));sd=schur['diagonalPosition'].astype(np.int64)-1
def product(values,columns,row_ids,vector,output,select_rows=None,select_columns=None):
 result=np.zeros((n,output))
 for first in range(0,len(columns),4096):
  ids=slice(first,first+4096);blocks=values[ids]
  if select_rows is not None:blocks=blocks[:,select_rows,:]
  if select_columns is not None:blocks=blocks[:,:,select_columns]
  np.add.at(result,row_ids[ids],np.einsum('kab,kb->ka',blocks,vector[columns[ids]]))
 return result
def base(vector,r,c):return product(cv,col,owners,vector,len(r),r,c)
def momentum(rhs,observe=False):
 forward=np.zeros((n,3));answer=np.zeros((n,3))
 with np.errstate(over='ignore',invalid='ignore'):
  for i in range(n):
   ids=slice(ptr[i],md[i]);forward[i]=rhs[i]-np.einsum('kab,kb->a',mv[ids],forward[col[ids]])
   if observe and not np.all(np.isfinite(forward[i])):return None,dict(stage='forward-lower-sum',algebraicRow=i+1,completedForwardRows=i,completedBackwardRows=0)
  for i in range(n-1,-1,-1):
   ids=slice(md[i]+1,ptr[i+1]);current=forward[i]-np.einsum('kab,kb->a',mv[ids],answer[col[ids]])
   if observe and not np.all(np.isfinite(current)):return None,dict(stage='backward-upper-sum',algebraicRow=i+1,completedForwardRows=n,completedBackwardRows=n-1-i)
   answer[i]=mi[i]@current
   if observe and not np.all(np.isfinite(answer[i])):return None,dict(stage='backward-diagonal-multiply',algebraicRow=i+1,completedForwardRows=n,completedBackwardRows=n-1-i)
 return (answer,None) if observe else answer

local= cv[diag][:,1:4,1:4];local_inverse=s['localInverse'].transpose(2,0,1)
local_defect=float(np.max(abs(np.einsum('kab,kbc->kac',local,local_inverse)-np.eye(3))))
assert local_defect<1e-8
independent_inverse=np.linalg.inv(local);inverse_difference=float(np.linalg.norm(independent_inverse-local_inverse)/np.linalg.norm(independent_inverse));assert inverse_difference<1e-8
del independent_inverse,local;gc.collect()
rhs=read(s['rhsFile'],['remaining','sourceFile']);assert rhs['sourceFile']==binding['sourceFile'];rm=rhs['remaining'][order,1:4];ym=momentum(rm);condensed=rhs['remaining'][order][:,[0,4]]-base(ym,[0,4],[1,2,3])
np.testing.assert_allclose(rm,s['momentumRhs'],rtol=0,atol=0);assert rel(ym,s['condensedMomentum'])<1e-8 and rel(condensed,s['condensedRhs'])<1e-8
condensed_momentum_gain=stable_norm(ym)/stable_norm(rm);condensed_momentum_residual=stable_norm(base(ym,[1,2,3],[1,2,3])-rm)/stable_norm(rm)
assert abs(condensed_momentum_gain/report['condensedMomentumToRhsGainSeconds']-1)<1e-8
failure_result=None
if report.get('momentumMapRejected',False):
 failed=s['failedCase'];bx=base(failed['value'],[1,2,3],[0,4]);assert rel(bx,failed['Bvalue'])<1e-8
 _,observed=momentum(bx,True);saved=failed['observation'];f=saved['firstFailure'];assert not saved['valueIsUsable'] and observed is not None
 for key in ('stage','algebraicRow'):assert observed[key]==f[key]
 for key in ('completedForwardRows','completedBackwardRows'):assert observed[key]==saved[key]
 block=f['factorBlocks'];dependency=f['dependencyValues'];initial=f['initialRowValue'].ravel()
 with np.errstate(over='ignore',invalid='ignore'):
  if f['stage']=='backward-diagonal-multiply':replayed=f['inverseDiagonal']@initial
  else:replayed=initial-np.einsum('abk,bk->a',block,dependency)
 expected=f['nonfiniteRowValue'].ravel();np.testing.assert_array_equal(np.isfinite(replayed),np.isfinite(expected));np.testing.assert_array_equal(np.isnan(replayed),np.isnan(expected))
 i=int(f['algebraicRow'])-1;original=int(order[i]+1);assert original==failed['originalCV']
 nf=int(binding['sourceFaces']);assert (original-1)%nf+1==failed['originalFace'] and (original-1)//nf+1==failed['layer']
 failure_result=dict(name=failed['name'],stage=observed['stage'],firstAlgebraicRow=i+1,originalCV=original,originalFace=int(failed['originalFace']),layer=int(failed['layer']),
   independentFirstNonfiniteStageAndRowExactlyMatch=True,storedFailureLedgerNonfiniteMasksExactlyMatch=True,BvalueRelativeDifference=rel(bx,failed['Bvalue']),
   maximumFiniteInputMagnitude=float(np.max(abs(bx))),maximumFiniteRecordedPrefixMagnitude=float(saved['maximumFiniteIntermediateMagnitude']))
 print('INDEPENDENT_FIXED_MOMENTUM_FIRST_NONFINITE_VERIFIED',flush=True)
result=[]
for case in rows(s['cases']):
 x=case['value'];bx=base(x,[1,2,3],[0,4]);ym=momentum(bx);dx=base(x,[0,4],[0,4]);ey=base(ym,[0,4],[1,2,3]);image=dx-ey
 local_y=np.einsum('kab,kb->ka',local_inverse,bx);local_coupling=base(local_y,[0,4],[1,2,3]);saved_image=product(sv,scol,so,x,2);aimage=base(ym,[1,2,3],[1,2,3]);detail=case['detail']
 errors=[rel(bx,detail['Bvalue']),rel(ym,detail['momentumValue']),rel(dx,detail['Dvalue']),rel(ey,detail['EmomentumValue']),rel(image,case['image']),rel(local_y,case['localMomentum']),rel(local_coupling,case['localCoupling']),rel(saved_image,case['savedLocalImage']),rel(aimage,case['momentumImage'])]
 assert max(errors)<1e-8
 discrepancy=float(np.linalg.norm(R*(saved_image-image))/np.linalg.norm(R*image));momentum_residual=float(np.linalg.norm(aimage-bx)/np.linalg.norm(bx));summary=case['summary']
 assert abs(discrepancy-summary['balancedLocalVersusConsistentRelativeDifference'])/max(1,discrepancy)<1e-8
 assert abs(momentum_residual-summary['rawMomentumRelativeResidual'])/max(1,momentum_residual)<1e-8
 balanced_direction=np.linalg.norm(x/C);local_gain=float(np.linalg.norm(R*saved_image)/balanced_direction);consistent_gain=float(np.linalg.norm(R*image)/balanced_direction)
 assert abs(local_gain-summary['balancedLocalImageToDirection'])/max(1,local_gain)<1e-8 and abs(consistent_gain-summary['balancedConsistentImageToDirection'])/max(1,consistent_gain)<1e-8
 result.append(dict(name=summary['name'],maximumCompleteIntermediateRelativeDifference=max(errors),balancedLocalVersusConsistentRelativeDifference=discrepancy,rawMomentumRelativeResidual=momentum_residual,balancedLocalImageToDirection=local_gain,balancedConsistentImageToDirection=consistent_gain))
 print('INDEPENDENT_CONSISTENT_SCHUR '+summary['name'],flush=True)

diagonal=sv[sd]*R[:,:,None]*C[:,None,:];np.testing.assert_allclose(diagonal,s['diagonal'].transpose(2,0,1),rtol=1e-14,atol=0)
singular=np.linalg.svd(diagonal,compute_uv=False);ratio=singular[:,1]/singular[:,0];difference=float(np.max(abs(ratio-s['sigmaRatio'])));assert difference<1e-14
d=report['localDiagnostics'];assert int(np.count_nonzero(ratio<1e-8))==d['ratiosBelow1e8'] and int(np.count_nonzero(ratio<1e-12))==d['ratiosBelow1e12']
worst=d['worstAlgebraicRows'].astype(int)-1;original=order[worst]+1;np.testing.assert_array_equal(original,d['worstOriginalCV']);nf=int(binding['sourceFaces']);np.testing.assert_array_equal((original-1)%nf+1,d['worstOriginalFace']);np.testing.assert_array_equal((original-1)//nf+1,d['worstLayers'])
if failure_result is not None:assert report['fixedMomentumSolves']==2 and report['consistentSchurActions']==0 and report['observedMomentumSolves']==1
else:assert report['fixedMomentumSolves']==11 and report['consistentSchurActions']==5 and report['observedMomentumSolves']==5
assert report['newActualNSEvaluations']==report['newProductionSolvedAreaM2']==0
assert all(sha(ROOT/path)==value for path,value in bindings.items())
out=dict(status='passed-independent-real-consistent-Schur-rejection-first-momentum-overflow-and-all-scaled-local-diagonals' if failure_result is not None else 'passed-independent-real-five-consistent-Schur-actions-and-all-scaled-local-diagonals',processID=os.getpid(),cellCount=n,cases=result,
 fixedMomentumFailure=failure_result,consistentSchurActionsCompleted=int(report['consistentSchurActions']),
 diagnosticOnlyRefinementKeepsAllOtherPayloadValues=refinement_verified,
 sameQCondensedMomentumInverseGainSeconds=condensed_momentum_gain,sameQRawMomentumResidualOfFactorSolve=condensed_momentum_residual,
 maximumOriginalLocalInverseDefect=local_defect,independentLocalInverseRelativeDifference=inverse_difference,allLocalScaledDiagonalSingularRatioMaximumDifference=difference,
 localScaledDiagonalSingularRatioRange=[float(ratio.min()),float(ratio.max())],localRatiosBelow1eMinus8=int(np.count_nonzero(ratio<1e-8)),localRatiosBelow1eMinus12=int(np.count_nonzero(ratio<1e-12)),
 completeGlobalNumericalRankOrUsableInverseProven=False,newSelectedJApplications=0,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
