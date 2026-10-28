# -*- coding: utf-8 -*-
"""独立完整C字段分解、二跳Schur装配和完整块对角三角解核对。"""
from pathlib import Path
import os,json,hashlib,sys
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum();target=VAL/(case+'_independent_field_schur.json');assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def matrices(f):
 n=int(f['cellCount']);d=int(f['dimension']);ptr=f['rowPointer'].astype(int)-1;col=f['column'].astype(int)-1;L=np.eye(d*n);U=np.zeros_like(L);mask=np.zeros((n,n),bool)
 for i in range(n):
  for p in range(ptr[i],ptr[i+1]):
   j=col[p];ir=i+np.arange(d)*n;jc=j+np.arange(d)*n;mask[i,j]=True
   if j<i:L[np.ix_(ir,jc)]=f['values'][:,:,p]
   else:U[np.ix_(ir,jc)]=f['values'][:,:,p]
 return L,U,mask
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
file=VAL/case/'field_schur_unit_evidence.mat';s=loadmat(file,simplify_cells=True);C=s['C'];results=[]
for c in np.atleast_1d(s['cases']).ravel():
 pre=c['pre'];n=int(pre['cellCount']);order=pre['order'].astype(int)-1;mi=np.r_[order+n,order+2*n,order+3*n];ei=np.r_[order,order+4*n]
 A=C[np.ix_(mi,mi)];B=C[np.ix_(mi,ei)];E=C[np.ix_(ei,mi)];D=C[np.ix_(ei,ei)]
 ML,MU,P=matrices(pre['momentumFactor']);SL,SU,P2=matrices(pre['schurFactor']);np.testing.assert_array_equal(P2,(P.astype(int)@P.astype(int))>0)
 local=np.zeros_like(A)
 for i in range(n):
  ids=i+np.arange(3)*n;local[np.ix_(ids,ids)]=np.linalg.inv(A[np.ix_(ids,ids)])
 S=D-E@local@B;maskM=np.tile(P,(3,3));maskS=np.tile(P2,(2,2))
 merr=float(np.linalg.norm((ML@MU-A)*maskM)/max(1,np.linalg.norm(A)));serr=float(np.linalg.norm((SL@SU-S)*maskS)/max(1,np.linalg.norm(S)))
 assert max(merr,serr)<1e-8 and rel(S,c['referenceS'])<2e-10
 storedB=np.zeros_like(B);storedE=np.zeros_like(E);ptr=pre['rowPointer'].astype(int)-1;col=pre['column'].astype(int)-1
 for i in range(n):
  for p in range(ptr[i],ptr[i+1]):
   j=col[p];storedB[np.ix_(i+np.arange(3)*n,j+np.arange(2)*n)]=pre['B'][:,:,p];storedE[np.ix_(i+np.arange(2)*n,j+np.arange(3)*n)]=pre['E'][:,:,p]
 assert rel(storedB,B)<2e-10 and rel(storedE,E)<2e-10
 rhs=c['rhs'][order];ym=np.linalg.solve(MU,np.linalg.solve(ML,rhs[:,1:4].ravel(order='F')));ys=np.linalg.solve(SU,np.linalg.solve(SL,rhs[:,[0,4]].ravel(order='F')-E@ym));ym-=np.linalg.solve(MU,np.linalg.solve(ML,B@ys))
 expected=np.zeros((n,5));expected[:,1:4]=ym.reshape((n,3),order='F');expected[:,[0,4]]=ys.reshape((n,2),order='F');physical=np.zeros_like(expected);physical[order]=expected;error=rel(c['value'],physical);assert error<1e-8 and c['linearity']<1e-9
 assert pre['constructionReport']['allOriginalCellsAndFiveEquationsKept'] and not pre['constructionReport']['physicalEquationChanged']
 results.append(dict(ordering=c['packed']['options']['ordering'],momentumRetainedLUProductDifference=merr,schurRetainedLUProductDifference=serr,
  completeSchurAssemblyDifference=rel(S,c['referenceS']),fullConservativeSolveDifference=error,scalarTriangularAssumptionUsed=False,baseBlocks=int(P.sum()),twoHopBlocks=int(P2.sum())))
assert s['report']['budgetRejected'] and s['report']['badRhsRejected']
out=dict(status='passed-independent-small-complete-five-field-Schur-assembly-and-full-block-solves',processID=os.getpid(),cases=results,
 newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
