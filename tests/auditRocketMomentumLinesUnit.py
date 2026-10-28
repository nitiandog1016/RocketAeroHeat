# -*- coding: utf-8 -*-
"""独立核对原面法向列、完整A和保存的物理Arnoldi方向，不把内部解当物理解。"""
from pathlib import Path
import os,json,hashlib,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from decimal import Decimal,localcontext
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1];assert name.replace('_','').isalnum();target=VAL/(name+'_independent_momentum_lines.json');assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
file=VAL/name/'momentum_lines_unit_evidence.mat';s=loadmat(file,simplify_cells=True);results=[]
for case in rows(s['cases']):
 f=case['factor'];m=case['matrix'];n=int(f['cellCount']);nf=int(f['sourceFaces']);nl=int(f['layers']);inverse=f['inverseLine'].transpose(2,0,1);mapping=f['inverseOrder'].astype(int)-1
 np.testing.assert_array_equal(np.sort(mapping),np.arange(n));ptr=m['rowPointer'].astype(int)-1;col=m['column'].astype(int)-1;A=np.zeros((3*n,3*n))
 for i in range(n):
  for k in range(ptr[i],ptr[i+1]):A[np.ix_(i+np.arange(3)*n,col[k]+np.arange(3)*n)]=m['values'][:,:,k]
 assert rel(A,case['A'])<1e-14;P=np.zeros_like(A);maximum=0.
 for face in range(nf):
  cells=mapping[face+np.arange(nl)*nf];ids=(cells[None,:]+np.arange(3)[:,None]*n).ravel(order='F');line=A[np.ix_(ids,ids)];exact=np.linalg.solve(line,np.eye(3*nl));P[np.ix_(ids,ids)]=exact
  maximum=max(maximum,float(np.max(abs(line@inverse[face]-np.eye(3*nl)))));assert rel(inverse[face],exact)<1e-9
 assert maximum<1e-8 and rel(P,case['P'])<1e-9
 b=case['rhs'].ravel(order='F');x=case['singleValue'].ravel(order='F');assert rel(x,P@b)<1e-10 and rel(case['singleImage'].ravel(order='F'),A@x)<1e-11
 solved=case['solved'].ravel(order='F');assert np.linalg.norm(A@solved-b)/np.linalg.norm(b)<1.1e-10 and rel(solved,np.linalg.solve(A,b))<1e-8
 k=case['detail']['krylov'];assert k['flag']==0 and case['linearity']<1e-10
 cycle=rows(k['cycles'])[0];Z=cycle['Z'];V=cycle['V'];H=cycle['smallH'];relation=rel(A@Z,V@H);assert relation<1e-8
 results.append(dict(maximumOriginalLineInverseDefect=maximum,fullAArnoldiRelationRelativeDifference=relation,relativeAResidual=float(np.linalg.norm(A@solved-b)/np.linalg.norm(b)),steps=int(k['arnoldiSteps'])))
assert s['report']['budgetRejected'] and s['report']['singularLineRejectedWithoutCompensation']
precision_result=None
if 'frozenFactor' in s:
 precision=loadmat(s['precisionFile'],simplify_cells=True)
 for k in range(3):assert np.max(abs(s['compensatedResiduals'][:,:,k]-precision['referenceResiduals'][:,:,k]))<=s['bounds'][k]
 f=s['frozenFactor'];A=precision['A'];X=np.asarray(f['inverseLine']).reshape(A.shape);N=len(A)
 with localcontext() as context:
  context.prec=80;AA=[[Decimal.from_float(float(x)) for x in row] for row in A];XX=[[Decimal.from_float(float(x)) for x in row] for row in X]
  exact=max(abs((Decimal(1) if i==j else Decimal(0))-sum((AA[i][q]*XX[q][j] for q in range(N)),Decimal(0))) for i in range(N) for j in range(N))
 assert float(exact)<1e-8 and float(exact)<=f['factorizationReport']['originalUnitMaximumInverseDefect']
 assert f['factorizationReport']['refinedLines']==1
 precision_result=dict(decimal80FinalInverseDefect=float(exact),certifiedUpperBound=float(f['factorizationReport']['originalUnitMaximumInverseDefect']),ordinaryBLASDefect=float(np.max(abs(A@X-np.eye(N)))),threshold=1e-8)
out=dict(status='passed-independent-two-orderings-normal-line-inverses-complete-A-and-physical-Arnoldi-reference',processID=os.getpid(),cases=results,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
if precision_result is not None:out['decimal80RefinedFrozenLine']=precision_result
if 'fallbackFactor' in s:
 f=s['fallbackFactor'];assert f['factorizationReport']['fallbackColumnCount']==1 and not f['factorizationReport']['rejectedFullLineInversesClaimedValid']
 np.testing.assert_array_equal(f['inverseLine'],np.eye(6));np.testing.assert_array_equal(s['fallbackValue'],s['fallbackRhs']);np.testing.assert_array_equal(s['fallbackImage'],np.tile(s['fallbackRhs'].sum(axis=0),(2,1)))
 out['explicitDiagonalFallbackKeepsAllOriginalACouplings']=True
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
