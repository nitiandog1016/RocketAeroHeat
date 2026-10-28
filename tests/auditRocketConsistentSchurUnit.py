# -*- coding: utf-8 -*-
"""独立重建小域块三角因子，核对一致Schur及其每个中间量。"""
from pathlib import Path
import os,json,hashlib,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1];assert name.replace('_','').isalnum()
target=VAL/(name+'_independent_consistent_schur.json');assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def dense(blocks,ptr,col,n):
 out=np.zeros((n*blocks.shape[0],n*blocks.shape[1]))
 for i in range(n):
  for k in range(ptr[i],ptr[i+1]):out[np.ix_(i+np.arange(blocks.shape[0])*n,col[k]+np.arange(blocks.shape[1])*n)]=blocks[:,:,k]
 return out
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
file=VAL/name/'consistent_schur_unit_evidence.mat';s=loadmat(file,simplify_cells=True);result=[]
for case in rows(s['cases']):
 op=case['operator'];n=int(op['cellCount']);ptr=op['rowPointer'].astype(int)-1;col=op['column'].astype(int)-1;mf=op['momentumFactor']
 factor=dense(mf['values'],ptr,col,n);L=np.eye(3*n);U=np.zeros_like(L)
 for i in range(n):
  for j in range(n):
   ids=np.ix_(i+np.arange(3)*n,j+np.arange(3)*n)
   if j<i:L[ids]=factor[ids]
   else:U[ids]=factor[ids]
 B=dense(op['B'],ptr,col,n);E=dense(op['E'],ptr,col,n);D=dense(op['D'],ptr,col,n)
 reference=D-E@np.linalg.solve(U,np.linalg.solve(L,B));x=case['value'].ravel(order='F');detail=case['detail'];bx=B@x;ym=np.linalg.solve(U,np.linalg.solve(L,bx));image=D@x-E@ym
 errors=[rel(reference,case['reference']),rel(image,case['image'].ravel(order='F')),rel(bx,detail['Bvalue'].ravel(order='F')),rel(ym,detail['momentumValue'].ravel(order='F')),rel(D@x,detail['Dvalue'].ravel(order='F')),rel(E@ym,detail['EmomentumValue'].ravel(order='F'))]
 assert max(errors)<1e-10 and case['linearity']<1e-10;result.append(dict(maximumCompleteIntermediateRelativeDifference=max(errors)))
if 'lower' in s:
 assert s['report']['finiteObserverExactlyMatchesDefault'] and s['report']['threeFirstOverflowStagesVerified'] and s['report']['oldNonfiniteRejectionPreserved']
 for key,stage,row in [('lower','forward-lower-sum',3),('upper','backward-upper-sum',1),('diagonalFailure','backward-diagonal-multiply',1)]:
  observed=s[key];assert not observed['valueIsUsable'];f=observed['firstFailure'];assert f['stage']==stage and f['algebraicRow']==row
  with np.errstate(over='ignore',invalid='ignore'):
   if stage=='backward-diagonal-multiply':value=f['inverseDiagonal']@f['initialRowValue']
   else:value=f['initialRowValue']-np.einsum('abk,bk->a',f['factorBlocks'].reshape((2,2,-1)),f['dependencyValues'].reshape((2,-1)))
  np.testing.assert_array_equal(np.isfinite(value),np.isfinite(f['nonfiniteRowValue']))
if 'diagonalTests' in s:
 singular=np.linalg.svd(s['diagonalTests'].transpose(2,0,1),compute_uv=False);reference=np.divide(singular[:,1],singular[:,0],out=np.zeros(len(singular)),where=singular[:,0]>0)
 assert np.max(abs(reference-s['ratios']))<1e-14 and np.max(abs(reference-s['referenceRatios']))<1e-14
out=dict(status='passed-independent-small-fixed-momentum-inverse-consistent-Schur-all-intermediates',processID=os.getpid(),cases=result,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
