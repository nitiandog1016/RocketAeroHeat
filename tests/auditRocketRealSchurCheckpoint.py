# -*- coding: utf-8 -*-
"""独立复算全部真实Schur输入块、拒绝行消元账及高精度2量主元；不求新流场。"""
from pathlib import Path
from decimal import Decimal,localcontext
import os,json,hashlib,gc,time
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';folder=VAL/'schurcheckpointreal01/checkpoint';target=VAL/'schurcheckpointreal01_independent_real_checkpoint.json';assert not target.exists();bindings={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def read(p,variables=None):
 p=Path(p);p=p if p.is_absolute() else ROOT/p;bindings[str(p.relative_to(ROOT))]=sha(p);return loadmat(p,simplify_cells=True,variable_names=variables)
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
assembly=folder/'schur_assembly_evidence.mat';a=read(assembly,['sourceBinding','algebraicOrder','assemblyReport','sourceIdentity']);binding=a['sourceBinding'];packedFile=Path(binding['packedFile']);assert sha(packedFile)==binding['packedFileSHA256']
p=read(packedFile,['packed'])['packed'];n=int(p['cellCount']);cp=p['rowPointer'].astype(np.int64)-1;cc=p['column'].astype(np.int64)-1;cd=p['diagonalPosition'].astype(np.int64)-1;cv=p['values'].transpose(2,0,1)
np.testing.assert_array_equal(a['algebraicOrder'],p['order']);nb=len(cc)
m=read(assembly,['momentumMatrix'])['momentumMatrix']
for key in ('rowPointer','column','diagonalPosition'):np.testing.assert_array_equal(m[key],p[key])
for first in range(0,nb,4096):np.testing.assert_array_equal(m['values'][:,:,first:first+4096].transpose(2,0,1),cv[first:first+4096,1:4,1:4])
del m;gc.collect()
for field,rr,cols in [('couplingB',[1,2,3],[0,4]),('couplingE',[0,4],[1,2,3])]:
 v=read(assembly,[field])[field]
 for first in range(0,nb,4096):np.testing.assert_array_equal(v[:,:,first:first+4096].transpose(2,0,1),cv[first:first+4096][:,rr,:][:,:,cols])
 del v;gc.collect()
s=read(assembly,['schurMatrix'])['schurMatrix'];sp=s['rowPointer'].astype(np.int64)-1;sc=s['column'].astype(np.int64)-1;sd=s['diagonalPosition'].astype(np.int64)-1;sv=s['values'].transpose(2,0,1)
localInverse=np.zeros((n,3,3));localDefect=0.
for first in range(0,n,4096):
 ids=slice(first,min(first+4096,n));D=cv[cd[ids],1:4,1:4];inv=np.linalg.inv(D);localInverse[ids]=inv;localDefect=max(localDefect,float(np.max(abs(D@inv-np.eye(3)))))
assert localDefect<1e-8
err2=den2=formation2=0.;maximumRowRelative=maximumFormationRelative=0.;timer=time.monotonic()
for i in range(n):
 begin,end=sp[i],sp[i+1];cols=sc[begin:end];locations={int(c):j for j,c in enumerate(cols)};expected=np.zeros((end-begin,2,2));formation=0.;reached=set()
 for e in range(cp[i],cp[i+1]):
  expected[locations[int(cc[e])]]=cv[e][np.ix_([0,4],[0,4])]
 formation=float(np.linalg.norm(expected))
 for e in range(cp[i],cp[i+1]):
  j=int(cc[e]);ids=slice(cp[j],cp[j+1]);targets=np.fromiter((locations[int(k)] for k in cc[ids]),dtype=np.int64);reached.update(map(int,cc[ids]))
  left=cv[e][np.ix_([0,4],[1,2,3])]@localInverse[j];right=cv[ids,1:4,:][:,:,[0,4]];product=np.einsum('ab,kbc->kac',left,right)
  expected[targets]-=product;formation+=float(np.linalg.norm(product))
 assert reached==set(map(int,cols))
 actual=sv[begin:end];difference=float(np.linalg.norm(expected-actual));norm=float(np.linalg.norm(actual));err2+=difference*difference;den2+=norm*norm;formation2+=formation*formation
 maximumRowRelative=max(maximumRowRelative,difference/max(1,norm));maximumFormationRelative=max(maximumFormationRelative,difference/max(1,formation))
 if (i+1)%50000==0:print('INDEPENDENT_REAL_SCHUR_ASSEMBLY rows=%d/%d'%(i+1,n),flush=True)
assemblyRelative=float(np.sqrt(err2/max(1,den2)));assert assemblyRelative<1e-10 and maximumFormationRelative<1e-10
del localInverse;gc.collect()
f=read(folder/'density_energy_failure.mat')['failure'];partial=f['partialFactor'];pv=partial['values'].transpose(2,0,1);rp=partial['rowPointer'].astype(np.int64)-1;pc=partial['column'].astype(np.int64)-1;pd=partial['diagonalPosition'].astype(np.int64)-1;row=int(f['algebraicRow'])-1
np.testing.assert_array_equal(rp,sp);np.testing.assert_array_equal(pc,sc);np.testing.assert_array_equal(pd,sd);assert f['assemblyFileSHA256']==sha(assembly) and f['sourceBinding']==binding
begin,end=sp[row],sp[row+1];raw=sv[begin:end].copy();np.testing.assert_array_equal(raw.transpose(1,2,0),f['rowInput']);work=raw.copy();locations={int(c):i for i,c in enumerate(sc[begin:end])};maximumProductBackward=0.
for step in rows(f['retainedEliminationLedger']):
 entry=int(step['lowerEntry'])-1;pivot=int(step['pivotAlgebraicRow'])-1;cache=f['cachedInverseDiagonal'][:,:,pivot];np.testing.assert_array_equal(cache,step['cachedPivotInverse']);mult=work[entry-begin]@cache;work[entry-begin]=mult
 upper=np.atleast_1d(step['upperEntries']).astype(int)-1;targets=np.atleast_1d(step['targetEntries']).astype(int)-1-begin;blocks=pv[upper];np.testing.assert_array_equal(blocks,np.asarray(step['upperBlocks']).reshape((2,2,len(upper)),order='F').transpose(2,0,1))
 products=np.einsum('ab,kbc->kac',mult,blocks);saved=np.asarray(step['retainedProducts']).reshape((2,2,len(upper)),order='F').transpose(2,0,1)
 backward=float(np.linalg.norm(products-saved)/max(1,np.linalg.norm(abs(mult)@abs(blocks))));maximumProductBackward=max(maximumProductBackward,backward);assert backward<1e-10
 work[targets]-=products
rowDifference=float(np.linalg.norm(work-pv[begin:end]));rowFormation=float(f['diagonalFormationNormBound']);assert rowDifference/max(1,rowFormation)<1e-10
assert f['ledgerExactlyReproducesPartialRow'] and f['originalInverseDefectThreshold']==1e-8 and f['maximumInverseDefect']>=1e-8 and not f['preconditionerUsable']
D=np.asarray(f['pivotD']);inv=np.asarray(f['inverseCandidate']);numpyDefect=float(np.max(abs(D@inv-np.eye(2))))
with localcontext() as ctx:
 ctx.prec=80;DD=[[Decimal.from_float(float(x)) for x in rr] for rr in D];II=[[Decimal.from_float(float(x)) for x in rr] for rr in inv]
 det=DD[0][0]*DD[1][1]-DD[0][1]*DD[1][0];exactInv=[[DD[1][1]/det,-DD[0][1]/det],[-DD[1][0]/det,DD[0][0]/det]]
 def residual(candidate):return max(abs(sum(DD[i][k]*candidate[k][j] for k in range(2))-Decimal(int(i==j))) for i in range(2) for j in range(2))
 exactStored=float(residual(II));rounded=[[Decimal.from_float(float(exactInv[i][j])) for j in range(2)] for i in range(2)];exactRounded=float(residual(rounded));detText=str(det)
 rawDiag=sv[sd[row]];growth=float(np.linalg.norm(D)/max(np.finfo(float).tiny,np.linalg.norm(rawDiag)))
assert int(p['order'][row])==f['originalControlVolume']==61925 and f['originalFace']==2669 and f['layer']==4
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-all-real-Schur-input-coefficients-and-rejected-pivot-ledger-audit',processID=os.getpid(),cellCount=n,baseBlocks=nb,schurBlocks=len(sc),
 allMomentumAndCrossCouplingInputValuesExactlyKept=True,allSchurPatternRowsExactlyTwoHop=True,allSchurCoefficientRelativeDifference=assemblyRelative,
 maximumSchurRowRelativeDifference=maximumRowRelative,maximumSchurFormationRelativeDifference=maximumFormationRelative,maximumIndependentLocalMomentumInverseDefect=localDefect,
 maximumRejectedRowProductBackwardDifference=maximumProductBackward,rejectedRowRecomputeRelativeToFormation=rowDifference/max(1,rowFormation),
 rejectedAlgebraicRow=row+1,originalControlVolume=int(f['originalControlVolume']),originalFace=int(f['originalFace']),layer=int(f['layer']),
 actualMatlabInverseDefect=float(f['maximumInverseDefect']),independentFloatInverseProductDefect=numpyDefect,originalThreshold=1e-8,pivotRcond=float(f['pivotRcond']),
 decimal80StoredInverseResidual=exactStored,decimal80RoundedExactInverseResidual=exactRounded,exactDoubleInputDeterminant=detText,
 pivotNormGrowthFromOriginalSchurDiagonal=growth,originalSchurDiagonal=rawDiag.tolist(),rejectedPivot=D.tolist(),diagonalCancellationRatio=float(f['diagonalCancellationRatio']),
 priorPartialFactorAllRowsIndependentlyReplayed=False,fullUsablePreconditionerProduced=False,newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
 auditSeconds=time.monotonic()-timer,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
