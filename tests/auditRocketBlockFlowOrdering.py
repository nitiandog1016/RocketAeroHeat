# -*- coding: utf-8 -*-
"""独立稠密C、给定方向投影、全部原块和重排后三角解核对。"""
from pathlib import Path
import hashlib,json,os,sys
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_flow_order_unit.json');assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
file=VAL/case/'block_flow_order_unit_evidence.mat';s=loadmat(file,simplify_cells=True);n=len(s['centroid']);d=s['direction']/np.linalg.norm(s['direction']);projection=s['centroid']@d;expected=np.lexsort((np.arange(n),projection));results=[]
for c in np.atleast_1d(s['cases']).ravel():
 p=c['ordered'];o=c['original'];f=c['factor'];order=p['order'].astype(int)-1;np.testing.assert_array_equal(order,expected)
 mapping=p['flowReorderingReport']['newToOldBlockPosition'].astype(int)-1;np.testing.assert_array_equal(np.sort(mapping),np.arange(len(mapping)));np.testing.assert_array_equal(p['values'],o['values'][:,:,mapping])
 ptr=p['rowPointer'].astype(int)-1;col=p['column'].astype(int)-1;A=np.zeros((5*n,5*n));L=np.eye(5*n);U=np.zeros_like(A);mask=np.zeros((n,n),bool)
 for i in range(n):
  for entry in range(ptr[i],ptr[i+1]):
   j=col[entry];rows=slice(5*i,5*i+5);cols=slice(5*j,5*j+5);A[rows,cols]=p['values'][:,:,entry];mask[i,j]=True
   if j<i:L[rows,cols]=f['values'][:,:,entry]
   else:U[rows,cols]=f['values'][:,:,entry]
 indices=(order[:,None]+np.arange(5)*n).ravel();reference=s['C'][np.ix_(indices,indices)];error=float(np.linalg.norm(A-reference)/max(1,np.linalg.norm(reference)));assert error<2e-10
 referenceLU=reference.copy()
 for i in range(n):
  rows=slice(5*i,5*i+5)
  for j in np.flatnonzero(mask[i,:i]):
   js=slice(5*j,5*j+5);m=np.linalg.solve(referenceLU[js,js].T,referenceLU[rows,js].T).T;referenceLU[rows,js]=m
   for k in np.flatnonzero(mask[j,j+1:])+j+1:
    if mask[i,k]:referenceLU[rows,5*k:5*k+5]-=m@referenceLU[js,5*k:5*k+5]
 factors=L+U-np.eye(5*n);factor_error=float(np.linalg.norm(factors-referenceLU)/max(1,np.linalg.norm(referenceLU)));assert factor_error<1e-8
 rhs=c['rhs'].ravel(order='F')[indices];observed=c['solution'].ravel(order='F')[indices];direct=np.linalg.solve(U,np.linalg.solve(L,rhs));solve_error=float(np.linalg.norm(observed-direct)/max(1,np.linalg.norm(direct)));assert solve_error<1e-8
 results.append(dict(originalOrdering=o['options']['ordering'],completeCRelativeDifference=error,independentFactorRelativeDifference=factor_error,denseTriangularSolveRelativeDifference=solve_error))
assert np.all(s['report']['rejected']);np.testing.assert_array_equal(s['tied']['order'].astype(int)-1,np.arange(n))
out=dict(status='passed-independent-flow-ordering-projection-full-C-all-values-factors-and-triangular-solves',processID=os.getpid(),cases=results,
 newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
