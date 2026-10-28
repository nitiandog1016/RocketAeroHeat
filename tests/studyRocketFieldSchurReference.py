# -*- coding: utf-8 -*-
"""冻结小域完整C的三动量/密度能量Schur分裂参考，不修改物理方程或构造新J。"""
from pathlib import Path
import os,json,hashlib,sys,inspect
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat,savemat
from scipy.linalg import solve_triangular
from scipy.sparse.linalg import gmres
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1] if len(sys.argv)>1 else 'fieldschurref01';assert name.replace('_','').isalnum();folder=VAL/name;assert not folder.exists();folder.mkdir()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
file=VAL/'blockiluunit01/block_ilu_unit_evidence.mat';s=loadmat(file,simplify_cells=True);C=s['C'];n=int(s['cases'][0]['packed']['cellCount']);p=s['cases'][0]['packed'];ptr=p['rowPointer'].astype(int)-1;col=p['column'].astype(int)-1
P=np.zeros((n,n),bool)
for i in range(n):P[i,col[ptr[i]:ptr[i+1]]]=True
P2=(P.astype(int)@P.astype(int))>0;assert np.all(P2[P]);mi=np.r_[np.arange(n)+n,np.arange(n)+2*n,np.arange(n)+3*n];ei=np.r_[np.arange(n),np.arange(n)+4*n]
A=C[np.ix_(mi,mi)];B=C[np.ix_(mi,ei)];E=C[np.ix_(ei,mi)];D=C[np.ix_(ei,ei)];Mlocal=np.zeros_like(A)
for i in range(n):
 ids=i+np.arange(3)*n;Mlocal[np.ix_(ids,ids)]=np.linalg.inv(A[np.ix_(ids,ids)])
Slocal=D-E@Mlocal@B;Ainv=np.linalg.inv(A);Sexact=D-E@Ainv@B
def ilu_inverse(matrix,mask,dim):
 index=(np.arange(n)[:,None]+np.arange(dim)*n).ravel();W=matrix[np.ix_(index,index)].copy();L=np.eye(n*dim)
 for i in range(n):
  ir=slice(dim*i,dim*i+dim)
  for j in np.flatnonzero(mask[i,:i]):
   jr=slice(dim*j,dim*j+dim);mult=np.linalg.solve(W[jr,jr].T,W[ir,jr].T).T;L[ir,jr]=mult;W[ir,jr]=0
   for k in np.flatnonzero(mask[j,j+1:])+j+1:
    if mask[i,k]:W[ir,dim*k:dim*k+dim]-=mult@W[jr,dim*k:dim*k+dim]
 left=solve_triangular(L,np.eye(n*dim),lower=True,unit_diagonal=True)
 # W是块上三角；每个对角块为完整矩阵，不能当标量上三角跳过块内下半部。
 inv=np.linalg.solve(W,left);assert np.linalg.norm(W@inv-left)/max(1,np.linalg.norm(left))<1e-8
 physical=np.zeros_like(inv);physical[np.ix_(index,index)]=inv
 return physical,L,W
Milu,ML,MU=ilu_inverse(A,P,3);Silu,SL,SU=ilu_inverse(Slocal,P2,2)
def inverse(M,S):
 out=np.zeros_like(C);out[np.ix_(ei,ei)]=S;out[np.ix_(ei,mi)]=-S@E@M
 out[np.ix_(mi,ei)]=-M@B@S;out[np.ix_(mi,mi)]=M+M@B@S@E@M;return out
old=s['cases'][1];oldIndex=(old['factor']['order'].astype(int)[:,None]-1+np.arange(5)*n).ravel();oldInverse=np.zeros_like(C)
oldInverse[np.ix_(oldIndex,oldIndex)]=np.linalg.solve(old['U'],solve_triangular(old['L'],np.eye(5*n),lower=True,unit_diagonal=True))
savedRhs=old['rhs'].ravel(order='F');savedValue=old['value'].ravel(order='F');oldReplay=float(np.linalg.norm(oldInverse@savedRhs-savedValue)/max(1,np.linalg.norm(savedValue)));assert oldReplay<2e-10
refs=[('exactCReference',np.linalg.inv(C)),('originalFiveBlockILU0',oldInverse),('exactFullSchurReference',inverse(Ainv,np.linalg.inv(Sexact))),('localMomentumExactLocalSchur',inverse(Mlocal,np.linalg.inv(Slocal))),('exactMomentumExactLocalSchur',inverse(Ainv,np.linalg.inv(Slocal))),('ILUMomentumExactLocalSchur',inverse(Milu,np.linalg.inv(Slocal))),('ILUMomentumILULocalSchur',inverse(Milu,Silu))]
cases=[];rhs0=s['cases'][0]['rhs'].ravel(order='F');rhs1=np.sin(np.arange(1,5*n+1)*.37)
for name,pre in refs:
 CP=C@pre;sv=np.linalg.svd(CP,compute_uv=False);trials=[]
 for rhs in (rhs0,rhs1):
  history=[];tolerance={'rtol':1e-8} if 'rtol' in inspect.signature(gmres).parameters else {'tol':1e-8}
  y,flag=gmres(CP,rhs,restart=40,maxiter=1,atol=0,callback=history.append,callback_type='pr_norm',**tolerance);x=pre@y;res=np.linalg.norm(rhs-C@x)/np.linalg.norm(rhs)
  trials.append(dict(iterations=len(history),flag=int(flag),actualRelativeCResidual=float(res)))
 cases.append(dict(name=name,conditionOfCPreconditioned=float(sv[0]/sv[-1]),inverseNorm=float(np.linalg.norm(pre,2)),trials=trials))
direct_difference=float(np.linalg.norm(refs[0][1]-refs[2][1])/max(1,np.linalg.norm(refs[0][1])));assert direct_difference<2e-8
savemat(folder/'field_schur_dense_evidence.mat',dict(C=C,A=A,B=B,E=E,D=D,Mlocal=Mlocal,Slocal=Slocal,Sexact=Sexact,Milu=Milu,Silu=Silu,ML=ML,MU=MU,SL=SL,SU=SU,basePattern=P,schurPattern=P2,momentumIndices=mi+1,scalarIndices=ei+1),do_compression=True)
out=dict(status='tested-small-dense-five-conservation-momentum-density-energy-Schur-reference-only',processID=os.getpid(),cellCount=n,
 basePatternBlocks=int(P.sum()),twoHopSchurPatternBlocks=int(P2.sum()),oldFiveBlockSolveReplayDifference=oldReplay,scalarUpperTriangularAssumptionUsed=False,exactSchurInverseDifference=direct_difference,localSchurCondition=float(np.linalg.cond(Slocal)),exactSchurCondition=float(np.linalg.cond(Sexact)),cases=cases,
 notImplementedForRealDomain=True,matlabSchurBackendImplemented=False,newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
 inputSHA256={str(file.relative_to(ROOT)):sha(file)},outputEvidenceSHA256=sha(folder/'field_schur_dense_evidence.mat'),scriptSHA256=sha(__file__))
(folder/'field_schur_dense_summary.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
