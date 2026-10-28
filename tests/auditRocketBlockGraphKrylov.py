# -*- coding: utf-8 -*-
"""独立稠密C与L/U验证块ILU内部Krylov的物理Z、最小残差和工作量。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_block_krylov.json');assert not target.exists()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
def rows(x):return [x] if isinstance(x,dict) else list(np.atleast_1d(x).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
file=VAL/case/'block_krylov_unit_evidence.mat';s=loadmat(file,simplify_cells=True);C=s['C'];f=s['factor'];n=len(s['Q']);order=f['order'].astype(int)-1
indices=(order[:,None]+np.arange(5)*n).ravel();L=np.eye(5*n);U=np.zeros((5*n,5*n));ptr=f['rowPointer'].astype(int)-1;cols=f['column'].astype(int)-1
for i in range(n):
    for p in range(ptr[i],ptr[i+1]):
        j=cols[p]
        if j<i:L[5*i:5*i+5,5*j:5*j+5]=f['values'][:,:,p]
        else:U[5*i:5*i+5,5*j:5*j+5]=f['values'][:,:,p]
b=s['rhs'].ravel(order='F');x=s['value'].ravel(order='F');r=s['detail']['krylov'];maximum_relation=maximum_P=0.
for c in rows(r['cycles']):
    j=int(c['arnoldiSteps']);V=np.asarray(c['V']).reshape((5*n,j+1),order='F');Z=np.asarray(c['Z']).reshape((5*n,j),order='F');H=np.asarray(c['smallH']).reshape((j+1,j),order='F')
    maximum_relation=max(maximum_relation,rel(C@Z,V@H));assert maximum_relation<2e-10
    for col in range(j):
        ref=np.zeros(5*n);ref[indices]=np.linalg.solve(U,np.linalg.solve(L,V[indices,col]));maximum_P=max(maximum_P,rel(ref,Z[:,col]))
    assert maximum_P<2e-10
    weights=np.linalg.lstsq(C@Z,c['rStart'],rcond=1e-13)[0];assert np.linalg.norm(C@Z@(weights-np.atleast_1d(c['weights'])))/max(1,np.linalg.norm(b))<2e-10
    assert np.linalg.norm(c['actualRemaining']-(b-C@c['candidate']))/max(1,np.linalg.norm(b))<2e-10
residual=float(np.linalg.norm(b-C@x)/np.linalg.norm(b));assert residual<1.1e-10 and rel(x,np.linalg.solve(C,b))<2e-8
assert r['modelApplications']==r['arnoldiSteps']+1 and r['preconditionerApplications']==r['arnoldiSteps']
assert r['basisStorageUpperBoundBytes']<=r['maximumBasisBytes'] and s['report']['budgetRejected'] and np.all(s['zero']==0)
assert not s['detail']['outerFixedLinearMapClaimed'] and s['detail']['requiresFlexibleOuterMethod']
source=Path(s['sourceFile']);out=dict(status='passed-independent-complete-C-block-ILU-Krylov-Arnoldi-and-dense-reference',case=case,processID=os.getpid(),
    relativeCResidual=residual,maximumCArnoldiRelationDifference=maximum_relation,maximumILUPhysicalDirectionDifference=maximum_P,
    arnoldiSteps=int(r['arnoldiSteps']),CApplications=int(r['modelApplications']),ILUApplications=int(r['preconditionerApplications']),
    newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (file,source)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
