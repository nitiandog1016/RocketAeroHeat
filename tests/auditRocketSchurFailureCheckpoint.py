# -*- coding: utf-8 -*-
"""独立核对小域组装检查点、拒绝主元及局部消元账；不放行部分因子。"""
from pathlib import Path
import os,json,hashlib,sys
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum();folder=VAL/case;target=VAL/(case+'_independent_checkpoint.json');assert not target.exists();bindings={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(p,variables=None):
 p=Path(p);p=p if p.is_absolute() else ROOT/p;bindings[str(p.relative_to(ROOT))]=sha(p);return loadmat(p,simplify_cells=True,variable_names=variables)
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def matrix(f):
 n=int(f['cellCount']);d=int(f['dimension']);a=np.zeros((d*n,d*n));ptr=f['rowPointer'].astype(int)-1;col=f['column'].astype(int)-1
 for i in range(n):
  for p in range(ptr[i],ptr[i+1]):a[np.ix_(i+np.arange(d)*n,col[p]+np.arange(d)*n)]=f['values'][:,:,p]
 return a
u=read(folder/'schur_checkpoint_unit_evidence.mat');parent=read(u['sourceFile']);results=[]
for j,c in enumerate(rows(parent['cases']),1):
 a=read(folder/('assembly%d'%j)/'schur_assembly_evidence.mat');np.testing.assert_array_equal(a['algebraicOrder'],c['pre']['order'])
 A=matrix(a['momentumMatrix']);S=matrix(a['schurMatrix']);assert np.linalg.norm(A-c['A'])/max(1,np.linalg.norm(A))<1e-12
 error=float(np.linalg.norm(S-c['referenceS'])/max(1,np.linalg.norm(S)));assert error<2e-10
 assert a['assemblyReport']['checkpointOnly'] and not a['assemblyReport']['physicalSolution'];results.append(dict(case=j,completeSchurCheckpointDifference=error))
f=read(folder/'synthetic_failure.mat')['failure'];raw=u['raw'];row=int(f['algebraicRow'])-1;ptr=raw['rowPointer'].astype(int)-1;diag=raw['diagonalPosition'].astype(int)-1
work=raw['values'][:,:,ptr[row]:ptr[row+1]].copy();np.testing.assert_array_equal(work,f['rowInput']);maximum=0.
for step in rows(f['retainedEliminationLedger']):
 lower=int(step['lowerEntry'])-1-ptr[row];before=work[:,:,lower].copy();mult=before@step['cachedPivotInverse'];np.testing.assert_allclose(mult,step['multiplier'],rtol=1e-12,atol=1e-14);work[:,:,lower]=mult
 targets=np.atleast_1d(step['targetEntries']).astype(int)-1-ptr[row];upper=np.asarray(step['upperBlocks']).reshape((2,2,len(targets)),order='F');products=np.einsum('ab,bck->ack',mult,upper)
 expected=np.asarray(step['retainedProducts']).reshape(products.shape,order='F');maximum=max(maximum,float(np.linalg.norm(products-expected)));assert maximum<1e-12
 work[:,:,targets]-=products
assert np.linalg.norm(work-f['actualPartialRow'])<1e-12 and f['ledgerExactlyReproducesPartialRow']
np.testing.assert_allclose(work[:,:,diag[row]-ptr[row]],f['pivotD'],rtol=1e-12,atol=1e-14)
defect=float(np.max(abs(f['pivotD']@f['inverseCandidate']-np.eye(2))));assert defect>=1e-8 and f['maximumInverseDefect']>=f['originalInverseDefectThreshold']==1e-8
assert f['partialFactor']['kind']=='partial-field-ILU-before-rejected-inverse-commit' and not f['preconditionerUsable'] and not f['physicalSolution']
assert u['report']['oldDefaultAndExplicitOutputsExactlyReplayed'] and u['report']['partialFactorSolveRejected']
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-small-assembly-checkpoints-and-original-rejected-pivot-ledger',processID=os.getpid(),assemblyChecks=results,
 maximumRetainedProductDifference=maximum,independentSyntheticInverseDefect=defect,originalThreshold=1e-8,partialFactorNotUsable=True,
 newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
