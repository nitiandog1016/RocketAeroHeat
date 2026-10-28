# -*- coding: utf-8 -*-
"""独立核对标量主元方向、零丢弃完整S及五量回代，不把近似因子当新物理解。"""
from pathlib import Path
import os,json,hashlib,sys
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum();target=VAL/(case+'_independent_pivoted_schur.json');assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))

def factor_matrix(f,S):
 n=int(f['cellCount']);scalar=f['scalarOrder'].astype(int)-1 if 'scalarOrder' in f else np.arange(2*n)
 ordered=S[np.ix_(scalar,scalar)]
 assert ('rowScale' in f)==('columnScale' in f)
 if 'rowScale' in f:
  r=np.asarray(f['rowScale']);c=np.asarray(f['columnScale']);assert np.all(np.isfinite(r)&(r>0)) and np.all(np.isfinite(c)&(c>0))
  referenceR=1/np.max(abs(S),axis=1);referenceC=1/np.max(abs(referenceR[:,None]*S),axis=0)
  np.testing.assert_allclose(r,referenceR[scalar],rtol=1e-14,atol=0);np.testing.assert_allclose(c,referenceC[scalar],rtol=1e-14,atol=0)
  ordered=r[:,None]*ordered*c[None,:]
 return scalar,ordered

def independent_solve(f,b):
 scalar=f['scalarOrder'].astype(int)-1;perm=f['rowPermutation'].astype(int)-1;v=b[scalar]
 if 'rowScale' in f:v=f['rowScale']*v
 y=np.linalg.solve(f['U'].toarray(),np.linalg.solve(f['L'].toarray(),v[perm]))
 if 'columnScale' in f:y=f['columnScale']*y
 x=np.zeros_like(b);x[scalar]=y;return x
file=VAL/case/'pivoted_schur_unit_evidence.mat';s=loadmat(file,simplify_cells=True);parentFile=Path(s['parentFile']);parent=loadmat(parentFile,simplify_cells=True);results=[]
for j,c in enumerate(rows(s['cases'])):
 f=c['factor'];S=c['S'];n=int(f['cellCount']);perm=f['rowPermutation'].astype(int)-1;L=f['L'].toarray();U=f['U'].toarray();np.testing.assert_array_equal(np.sort(perm),np.arange(2*n))
 scalar,orderedS=factor_matrix(f,S);np.testing.assert_array_equal(np.sort(scalar),np.arange(2*n))
 relation=rel(L@U,orderedS[perm]);assert relation<1e-8 and np.allclose(L,np.tril(L)) and np.allclose(U,np.triu(U)) and np.all(np.diag(L)==1)
 b=c['rhs'].ravel(order='F');observed=c['value'].ravel(order='F');assert rel(observed,np.linalg.solve(S,b))<1e-8
 assert rel(S@observed,b)<1e-8 and rel(independent_solve(f,b),observed)<1e-8
 p=c['pre'];mf=p['momentumFactor'];ML=np.eye(3*n);MU=np.zeros((3*n,3*n));ptr=mf['rowPointer'].astype(int)-1;col=mf['column'].astype(int)-1
 for i in range(n):
  for entry in range(ptr[i],ptr[i+1]):
   k=col[entry];ir=i+np.arange(3)*n;kc=k+np.arange(3)*n
   if k<i:ML[np.ix_(ir,kc)]=mf['values'][:,:,entry]
   else:MU[np.ix_(ir,kc)]=mf['values'][:,:,entry]
 original=rows(parent['cases'])[j];B=original['B'];E=original['E'];order=p['order'].astype(int)-1;r=c['fullRhs'][order]
 ym=np.linalg.solve(MU,np.linalg.solve(ML,r[:,1:4].ravel(order='F')));ys=np.linalg.solve(S,r[:,[0,4]].ravel(order='F')-E@ym);ym-=np.linalg.solve(MU,np.linalg.solve(ML,B@ys))
 expected=np.zeros((n,5));expected[:,1:4]=ym.reshape((n,3),order='F');expected[:,[0,4]]=ys.reshape((n,2),order='F');physical=np.zeros_like(expected);physical[order]=expected
 full_error=rel(c['fullValue'],physical);assert full_error<1e-8 and c['linearity']<1e-9
 approximate=c['approximate'];aL=approximate['L'];aU=approximate['U'];ap=approximate['rowPermutation'].astype(int)-1;aso,aS=factor_matrix(approximate,S);errors=[]
 for k,probe in enumerate(rows(approximate['factorizationReport']['probes']),1):
  v=np.sin(np.arange(1,2*n+1)*.37) if k==1 else np.cos(np.arange(1,2*n+1)*.23);actual=(aS@v)[ap];error=float(np.linalg.norm(actual-aL@(aU@v))/max(np.finfo(float).tiny,np.linalg.norm(actual)));assert abs(error-probe['relativePAminusLUAction'])<1e-10;errors.append(error)
 assert not f['factorizationReport']['zeroDiagonalReplacementEnabled'] and not f['factorizationReport']['diagonalCompensationApplied'] and not f['factorizationReport']['oldRejectedBlockInverseClaimedValid']
 results.append(dict(case=j,completePSminusLURelativeDifference=relation,fullFiveConservationSolveDifference=full_error,approximatePAminusLUProbeErrors=errors))
assert s['report']['forcedRowPivotVerified'] and s['report']['budgetRejected']
forced=s['forced'];assert np.any(forced['rowPermutation'].astype(int)-1!=np.arange(2));assert np.linalg.norm(s['toy']['values']@s['toyValue'].ravel()-np.array([2,3]))<1e-12
stress=[]
if 'stressCases' in s:
 assert s['report']['defaultFullOutputExceptTimingReplayed'] and s['report']['oldFieldMajorSolveReplayed'] and s['report']['zeroRowsAndColumnsRejected'] and s['report']['nonfiniteScalesAndUnderflowRejected']
 assert s['report']['incompleteScaleStateRejected'] and s['report']['batchReport']['blockCount']>4096
 for c in rows(s['stressCases']):
  f=c['factor'];_,balanced=factor_matrix(f,s['stressS']);x=independent_solve(f,s['stressRhs']);np.testing.assert_allclose(x,c['value'].ravel(order='F'),rtol=1e-10,atol=0)
  err=rel(s['rawColumns']*x,s['knownBalanced']);back=float(np.max(abs(s['stressS']@x-s['stressRhs'])/(abs(s['stressS'])@abs(x)+abs(s['stressRhs']))))
  assert err<1e-9 and back<1e-12;stress.append(dict(balancedKnownSolutionRelativeDifference=err,componentwiseOriginalBackwardError=back))
out=dict(status='passed-independent-pivoted-Schur-permutation-full-S-five-conservation-and-no-diagonal-replacement',processID=os.getpid(),cases=results,
 equilibrationIndependentlyVerified=all('rowScale' in c['factor'] for c in rows(s['cases'])),stressCases=stress,
 newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (file,parentFile)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
