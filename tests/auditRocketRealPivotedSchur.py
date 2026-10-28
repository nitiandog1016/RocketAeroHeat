# -*- coding: utf-8 -*-
"""独立核对真实主元Schur的等价排列、完整五量映射与三条原C作用。"""
from pathlib import Path
import os,json,hashlib,sys,gc
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
from scipy.sparse import triu,tril
from scipy.sparse.linalg import spsolve_triangular
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum();folder=VAL/case;target=VAL/(case+'_independent_pivoted_schur.json');assert not target.exists();bindings={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(p,variables=None):
 p=Path(p);p=p if p.is_absolute() else ROOT/p;bindings[str(p.relative_to(ROOT))]=sha(p);return loadmat(p,simplify_cells=True,variable_names=variables)
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
report=json.loads((folder/'pivoted_schur_real_summary.json').read_text(encoding='utf-8-sig'));assert report['factorReady'] and report['wholePreconditionerReady'] and not report['candidateRejected']
t=read(folder/'pivoted_schur_factor_evidence.mat');f=t['factor'];binding=t['sourceBinding'];assembly=Path(t['assemblyFile']);n=int(f['cellCount']);N=2*n;scalar=f['scalarOrder'].astype(int)-1;perm=f['rowPermutation'].astype(int)-1;L=f['L'];U=f['U']
np.testing.assert_array_equal(np.sort(scalar),np.arange(N));np.testing.assert_array_equal(np.sort(perm),np.arange(N));assert triu(L,1).nnz==0 and tril(U,-1).nnz==0 and np.all(L.diagonal()==1)
assert ('rowScale' in f)==('columnScale' in f)
rowScale=np.asarray(f.get('rowScale',np.ones(N)));columnScale=np.asarray(f.get('columnScale',np.ones(N)));scaled='rowScale' in f
assert np.all(np.isfinite(rowScale)&(rowScale>0)) and np.all(np.isfinite(columnScale)&(columnScale>0))
assert np.all(np.isfinite(U.diagonal())) and np.all(U.diagonal()!=0);inverseDefect=float(np.max(abs(U.diagonal()*(1/U.diagonal())-1)));assert inverseDefect<1e-8
configuration=f['factorizationReport']['configuration'];assert configuration['milu']=='off' and configuration['udiag']==0 and not f['factorizationReport']['oldRejectedBlockInverseClaimedValid']
s=read(assembly,['schurMatrix'])['schurMatrix'];ptr=s['rowPointer'].astype(np.int64)-1;col=s['column'].astype(np.int64)-1;owners=np.repeat(np.arange(n,dtype=np.int32),np.diff(ptr));sv=s['values'].transpose(2,0,1)
if scaled:
 rmax=np.zeros((n,2));cmax=np.zeros((n,2))
 for a in range(2):
  for b in range(2):np.maximum.at(rmax[:,a],owners,abs(sv[:,a,b]))
 r=1/rmax
 for a in range(2):
  for b in range(2):np.maximum.at(cmax[:,b],col,abs(sv[:,a,b])*r[owners,a])
 np.testing.assert_allclose(rowScale,r.ravel(order='F')[scalar],rtol=1e-14,atol=0);np.testing.assert_allclose(columnScale,(1/cmax).ravel(order='F')[scalar],rtol=1e-14,atol=0)
 del rmax,cmax,r
def blockAction(values,columns,owners,vector,output):
 out=np.zeros((n,output))
 for first in range(0,len(columns),4096):
  ids=slice(first,first+4096);np.add.at(out,owners[ids],np.einsum('kab,kb->ka',values[ids],vector[columns[ids]]))
 return out
probes=[]
for j,saved in enumerate(rows(f['factorizationReport']['probes']),1):
 v=np.sin(np.arange(1,N+1)*.37) if j==1 else np.cos(np.arange(1,N+1)*.23);physical=np.zeros(N);physical[scalar]=columnScale*v
 image=(rowScale*blockAction(sv,col,owners,physical.reshape((n,2),order='F'),2).ravel(order='F')[scalar])[perm];product=L@(U@v)
 error=float(np.linalg.norm(image-product)/max(np.finfo(float).tiny,np.linalg.norm(image)));assert abs(error-saved['relativePAminusLUAction'])/max(1,error)<1e-8;probes.append(error)
factorReport=f['factorizationReport'];del s,sv,owners,t,f,L,U;gc.collect()
p=read(folder/'pivoted_field_schur_preconditioner.mat');pre=p['pre'];assert p['sourceBinding']==binding;order=pre['order'].astype(int)-1
sf=pre['schurFactor'];np.testing.assert_array_equal(sf['rowPermutation'].astype(int)-1,perm);np.testing.assert_array_equal(sf['scalarOrder'].astype(int)-1,scalar);L=sf['L'];U=sf['U']
np.testing.assert_array_equal(sf.get('rowScale',np.ones(N)),rowScale);np.testing.assert_array_equal(sf.get('columnScale',np.ones(N)),columnScale)
c=read(binding['packedFile'],['packed'])['packed'];assert sha(binding['packedFile'])==binding['packedFileSHA256'];np.testing.assert_array_equal(c['order'].astype(int)-1,order)
cp=c['rowPointer'].astype(np.int64)-1;cc=c['column'].astype(np.int64)-1;cv=c['values'].transpose(2,0,1);co=np.repeat(np.arange(n,dtype=np.int32),np.diff(cp))
for first in range(0,len(cc),4096):
 ids=slice(first,first+4096);np.testing.assert_array_equal(pre['B'][:,:,ids].transpose(2,0,1),cv[ids,1:4,:][:,:,[0,4]]);np.testing.assert_array_equal(pre['E'][:,:,ids].transpose(2,0,1),cv[ids][:,[0,4],1:4])
mf=pre['momentumFactor'];mp=mf['rowPointer'].astype(np.int64)-1;mc=mf['column'].astype(np.int64)-1;md=mf['diagonalPosition'].astype(np.int64)-1;mv=mf['values'].transpose(2,0,1);mi=mf['inverseDiagonal'].transpose(2,0,1)
np.testing.assert_array_equal(mp,cp);np.testing.assert_array_equal(mc,cc)
def momentum(b):
 forward=np.zeros((n,3));solution=np.zeros((n,3))
 for i in range(n):
  ids=slice(mp[i],md[i]);forward[i]=b[i]-np.einsum('kab,kb->a',mv[ids],forward[mc[ids]])
 for i in range(n-1,-1,-1):
  ids=slice(md[i]+1,mp[i+1]);solution[i]=mi[i]@(forward[i]-np.einsum('kab,kb->a',mv[ids],solution[mc[ids]]))
 return solution
def fullSolve(rhs):
 r=rhs[order];ym=momentum(r[:,1:4]);rs=r[:,[0,4]]-blockAction(pre['E'].transpose(2,0,1),mc,co,ym,2);b=(rowScale*rs.ravel(order='F')[scalar])[perm]
 y=spsolve_triangular(U,spsolve_triangular(L,b,lower=True,unit_diagonal=True),lower=False);x=np.zeros(N);x[scalar]=columnScale*y;ys=x.reshape((n,2),order='F')
 ym-=momentum(blockAction(pre['B'].transpose(2,0,1),mc,co,ys,3));sorted=np.zeros((n,5));sorted[:,1:4]=ym;sorted[:,[0,4]]=ys;physical=np.zeros_like(sorted);physical[order]=sorted;return physical
actions=read(folder/'pivoted_field_schur_actions.mat',['cases'])['cases'];result=[]
for j,item in enumerate(rows(actions)):
 v=item['value'];image=blockAction(cv,cc,co,v[order],5);physical=np.zeros_like(image);physical[order]=image
 error=float(np.linalg.norm(physical-item['Cvalue'])/max(1,np.linalg.norm(physical)));assert error<1e-8
 residual=float(np.linalg.norm(item['rhs']-physical)/np.linalg.norm(item['rhs']));assert abs(residual-item['summary']['relativeCResidual'])/max(1,residual)<1e-8
 result.append(dict(name=item['summary']['name'],independentCActionRelativeDifference=error,relativeCResidual=residual,directionNormToRhsSeconds=float(np.linalg.norm(v)/np.linalg.norm(item['rhs']))))
 if j==0:
  expected=fullSolve(item['rhs']);mapping=float(np.linalg.norm(expected-v)/max(1,np.linalg.norm(expected)));assert mapping<1e-8;result[-1]['independentWholeFiveFieldSolveRelativeDifference']=mapping;print('INDEPENDENT_PIVOTED_FULL_MAP complete',flush=True)
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-real-pivoted-Schur-scalar-ordering-five-field-map-and-three-complete-C-actions',processID=os.getpid(),cellCount=n,scalarUnknowns=N,
 rowAndScalarPermutationsVerified=True,reversibleEquilibrationIndependentlyVerified=scaled,maximumScalarInverseDefect=inverseDefect,noDiagonalReplacementOrCompensation=True,scalarPAminusLUProbeErrors=probes,
 actualFactorStorageBoundBytes=int(factorReport['factorStorageUpperBoundBytes']),caseResults=result,newSelectedJApplications=0,newActualNSEvaluations=0,
 originalNSAcceptanceClaimed=False,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
