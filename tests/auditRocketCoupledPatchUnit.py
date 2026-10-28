# -*- coding: utf-8 -*-
"""独立核对强连接组覆盖/连通、五量缓存逆和完整C的物理Arnoldi方向。"""
from pathlib import Path
import os,json,hashlib,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1];assert name.replace('_','').isalnum();target=VAL/(name+'_independent_coupled_patches.json');assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
file=VAL/name/'coupled_patches_unit_evidence.mat';s=loadmat(file,simplify_cells=True);results=[]
for case in rows(s['cases']):
 packed=case['packed'];n=int(packed['cellCount']);order=packed['order'].astype(int)-1;ptr=packed['rowPointer'].astype(int)-1;col=packed['column'].astype(int)-1;C=np.zeros((5*n,5*n));adj=np.zeros((n,n),bool)
 for i in range(n):
  for e in range(ptr[i],ptr[i+1]):
   j=col[e];ir=order[i]+np.arange(5)*n;ic=order[j]+np.arange(5)*n;C[np.ix_(ir,ic)]=packed['values'][:,:,e]
   if i!=j and np.any(packed['values'][:,:,e]):adj[order[i],order[j]]=adj[order[j],order[i]]=True
 partition=case['partition'];members=partition['members'].astype(int)-1;np.testing.assert_array_equal(np.sort(members),np.arange(n));P=np.zeros_like(C);maximum=0.
 for size,bucket in enumerate(rows(case['factor']['buckets']),1):
  groups=np.atleast_1d(bucket['patchIDs']);count=len(groups)
  if count==0:continue
  mem=np.asarray(bucket['members']).reshape((size,count),order='F').astype(int)-1;inverse=np.asarray(bucket['inverse']).reshape((5*size,5*size,count),order='F');fallback=np.atleast_1d(bucket['fallback']).astype(bool)
  for page in range(count):
   node=mem[:,page];visited={int(node[0])}
   while True:
    added={int(x) for x in node if any(adj[x,y] for y in visited)}-visited
    if not added:break
    visited|=added
   assert len(visited)==size and size<=partition['maximumPatchCells']
   ids=(node[None,:]+np.arange(5)[:,None]*n).ravel(order='F');block=C[np.ix_(ids,ids)]
   if fallback[page]:
    diagonal=np.zeros_like(block)
    for k in range(size):ii=slice(5*k,5*k+5);diagonal[ii,ii]=block[ii,ii]
    block=diagonal
   expected=np.linalg.solve(block,np.eye(5*size));assert rel(expected,inverse[:,:,page])<1e-8;P[np.ix_(ids,ids)]=expected;maximum=max(maximum,float(np.max(abs(block@inverse[:,:,page]-np.eye(5*size)))))
 assert maximum<1e-8 and rel(P,case['P'])<1e-8
 b=case['rhs'].ravel(order='F');assert rel(P@b,case['single'].ravel(order='F'))<1e-8
 x=case['value'].ravel(order='F');residual=float(np.linalg.norm(C@x-b)/np.linalg.norm(b));assert residual<1.1e-10 and rel(x,np.linalg.solve(C,b))<1e-8
 d=case['detail'];cycle=rows(d['krylov']['cycles'])[0];Z=cycle['Z'];V=cycle['V'];H=cycle['smallH'];relation=rel(C@Z,V@H);mapping=rel(P@V[:,:Z.shape[1]]/cycle['imageNormalizationScales'][None,:],Z);assert relation<1e-8 and mapping<1e-8
 results.append(dict(fullCResidual=residual,maximumLocalInverseDefect=maximum,allPhysicalArnoldiRelationDifference=relation,allLocalMapDifference=mapping,steps=int(d['krylov']['arnoldiSteps'])))
if 'toy' in s:
 np.testing.assert_array_equal(s['toyRhs'],s['toyValue']);np.testing.assert_array_equal(s['toyImage'],np.tile(s['toyRhs'].sum(axis=0),(2,1)));assert s['toyFactor']['factorizationReport']['fallbackPatchCount']==1
out=dict(status='passed-independent-small-connected-all-five-field-patches-and-complete-C-Krylov-reference',processID=os.getpid(),cases=results,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
