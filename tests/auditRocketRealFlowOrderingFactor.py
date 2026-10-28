# -*- coding: utf-8 -*-
"""独立流向排序、完整模式/系数及所有保留LU位置、局部逆和三条实际C作用。"""
from pathlib import Path
import os,json,hashlib,time
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';folder=VAL/'floworderfactor01';target=VAL/'floworderfactor01_independent_real_factor.json';assert not target.exists();bindings={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(p,variables=None):
 p=Path(p);p=p if p.is_absolute() else ROOT/p;bindings[str(p.relative_to(ROOT))]=sha(p);return loadmat(p,simplify_cells=True,variable_names=variables)
s=read(folder/'flow_ordered_factor_evidence.mat',['factor','geometryFile','packedFile']);f=s['factor'];g=read(s['geometryFile'],['centroid','unitFlowDirection']);p=read(s['packedFile'],['packed'])['packed'];del s
n=int(p['cellCount']);order=f['order'].astype(int)-1;projection=g['centroid']@g['unitFlowDirection'];expectedOrder=np.lexsort((np.arange(n),projection));np.testing.assert_array_equal(order,expectedOrder);del g
inverse=np.empty(n,int);inverse[order]=np.arange(n);oldOrder=p['order'].astype(int)-1;oldptr=p['rowPointer'].astype(np.int64)-1;oldcol=p['column'].astype(np.int64)-1
oldOwners=np.repeat(np.arange(n),np.diff(oldptr));rows=inverse[oldOrder[oldOwners]];columns=inverse[oldOrder[oldcol]];mapping=np.lexsort((columns,rows));expectedCol=columns[mapping];expectedPtr=np.r_[0,np.cumsum(np.bincount(rows,minlength=n))]
np.testing.assert_array_equal(f['flowReorderingReport']['newToOldBlockPosition'].astype(int)-1,mapping);ptr=f['rowPointer'].astype(np.int64)-1;col=f['column'].astype(np.int64)-1;diag=f['diagonalPosition'].astype(np.int64)-1
np.testing.assert_array_equal(ptr,expectedPtr);np.testing.assert_array_equal(col,expectedCol);del rows,columns,expectedCol,expectedPtr,oldOwners
pv=p['values'].transpose(2,0,1);fv=f['values'].transpose(2,0,1);den2=err2=outside2=0.;maximum_backward=inverse_defect=0.;retained=discarded=outside_blocks=0;timer=time.monotonic()
for i in range(n):
 begin,end=ptr[i],ptr[i+1];columns=col[begin:end];positions={int(c):j for j,c in enumerate(columns)};C=pv[mapping[begin:end]];product=np.zeros_like(C);bounds=np.linalg.norm(C,axis=(1,2));outside={}
 for low in range(begin,diag[i]):
  j=int(col[low]);mult=fv[low]
  for up in range(diag[j],ptr[j+1]):
   k=int(col[up]);v=mult@fv[up];where=positions.get(k)
   if where is not None:
    product[where]+=v;bounds[where]+=np.linalg.norm(mult)*np.linalg.norm(fv[up]);retained+=int(k>j)
   else:
    discarded+=int(k>j)
    if k in outside:outside[k]+=v
    else:outside[k]=v.copy()
 for up in range(diag[i],end):
  where=positions[int(col[up])];product[where]+=fv[up];bounds[where]+=np.linalg.norm(fv[up])
 error=product-C;err2+=float(np.sum(error**2));den2+=float(np.sum(C*C));maximum_backward=max(maximum_backward,float(np.max(np.linalg.norm(error,axis=(1,2))/np.maximum(bounds,np.finfo(float).tiny))))
 outside2+=sum(float(np.sum(v*v)) for v in outside.values());outside_blocks+=len(outside);inverse_defect=max(inverse_defect,float(np.max(abs(fv[diag[i]]@f['inverseDiagonal'][:,:,i]-np.eye(5)))))
 if (i+1)%50000==0:print('INDEPENDENT_FLOW_LU rows=%d/%d'%(i+1,n),flush=True)
relative=float(np.sqrt(err2/max(1,den2)));assert relative<1e-8 and maximum_backward<1e-8 and inverse_defect<1e-8
assert retained==f['factorizationReport']['retainedBlockProducts'] and discarded==f['factorizationReport']['discardedOutOfPatternFillProducts']
assert f['flowReorderingReport']['persistentFactorStorageUpperBoundBytes']<=p['options']['maximumPersistentBytes']
del f,fv;originalOwners=np.repeat(np.arange(n),np.diff(oldptr));a=read(folder/'flow_ordering_actions_evidence.mat',['cases'])['cases'];probes=[]
for c in np.atleast_1d(a).ravel():
 v=c['value'][oldOrder];out=np.zeros((n,5))
 for first in range(0,len(oldcol),4096):
  ids=slice(first,first+4096);np.add.at(out,originalOwners[ids],np.einsum('kab,kb->ka',pv[ids],v[oldcol[ids]]))
 physical=np.zeros((n,5));physical[oldOrder]=out;difference=float(np.linalg.norm(physical-c['Cvalue'])/max(1,np.linalg.norm(physical)));assert difference<2e-10
 residual=float(np.linalg.norm(c['rhs']-physical)/np.linalg.norm(c['rhs']));assert abs(residual-c['summary']['relativeCResidual'])/max(1,residual)<2e-10
 probes.append(dict(name=c['summary']['name'],independentCActionRelativeDifference=difference,relativeCResidual=residual,directionNormToRhsSeconds=float(np.linalg.norm(c['value'])/np.linalg.norm(c['rhs']))))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-real-flow-order-full-pattern-all-retained-LU-pivots-and-three-C-actions',processID=os.getpid(),cellCount=n,
 completePatternAndOriginalValuesKept=True,retainedLUProductRelativeDifference=relative,maximumCancellationAwareRetainedBlockBackwardError=maximum_backward,allPivotInverseDefect=inverse_defect,
 fullLUApproximationRelativeDifference=float(np.sqrt((err2+outside2)/max(1,den2))),retainedBlockProducts=retained,discardedFillProducts=discarded,outsidePatternLUBlocks=outside_blocks,probes=probes,auditSeconds=time.monotonic()-timer,
 correctFactorArithmeticDoesNotCertifyStableInverse=True,newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
