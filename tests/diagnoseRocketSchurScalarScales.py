# -*- coding: utf-8 -*-
"""核对原S零行/列和尺度跨度，判断ILUTP零主元后应检查的数值条件。"""
from pathlib import Path
import os,json,hashlib
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';target=VAL/'schur_scalar_scales01.json';assert not target.exists()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
file=VAL/'schurcheckpointreal01/checkpoint/schur_assembly_evidence.mat';s=loadmat(file,simplify_cells=True,variable_names=['schurMatrix','algebraicOrder','sourceBinding']);m=s['schurMatrix'];n=int(m['cellCount']);N=2*n;col=m['column'].astype(int)-1;ptr=m['rowPointer'].astype(np.int64)-1;diag=m['diagonalPosition'].astype(int)-1;owners=np.repeat(np.arange(n,dtype=np.int32),np.diff(ptr));v=m['values'].transpose(2,0,1)
rmax=np.zeros((2,n));cmax=np.zeros((2,n));rnorm2=rmax.copy();cnorm2=cmax.copy();rnnz=np.zeros((2,n),np.int64);cnnz=rnnz.copy()
for a in range(2):
 for b in range(2):
  absolute=abs(v[:,a,b]);np.maximum.at(rmax[a],owners,absolute);np.maximum.at(cmax[b],col,absolute)
  rnorm2[a]+=np.bincount(owners,weights=absolute**2,minlength=n);cnorm2[b]+=np.bincount(col,weights=absolute**2,minlength=n)
  rnnz[a]+=np.bincount(owners,weights=(absolute>0),minlength=n).astype(np.int64);cnnz[b]+=np.bincount(col,weights=(absolute>0),minlength=n).astype(np.int64)
rowMax=rmax.ravel();columnMax=cmax.ravel();zeroRows=np.flatnonzero(rowMax==0);zeroCols=np.flatnonzero(columnMax==0);dd=v[diag]
cellMatching=((dd[:,0,0]!=0)&(dd[:,1,1]!=0))|((dd[:,0,1]!=0)&(dd[:,1,0]!=0));proof=bool(np.all(cellMatching))
def span(x):
 positive=x[x>0];return dict(minimumPositive=float(positive.min()),medianPositive=float(np.median(positive)),maximum=float(positive.max()),maximumToMinimum=float(positive.max()/positive.min()))
potentialSmall=offdiag=0;columnNorm=np.sqrt(cnorm2)
for first in range(0,len(col),8192):
 ids=slice(first,first+8192)
 for a in range(2):
  for b in range(2):
   absolute=abs(v[ids,a,b]);mask=(absolute>0)&((owners[ids]!=col[ids])|(a!=b));offdiag+=int(np.count_nonzero(mask));potentialSmall+=int(np.count_nonzero(mask&(absolute<.01*columnNorm[b,col[ids]])))
examples=[];order=s['algebraicOrder'].astype(int);nf=int(s['sourceBinding']['sourceFaces'])
for index in zeroRows[:50]:
 cv=int(order[index%n]);examples.append(dict(scalarRow=int(index+1),field='density' if index<n else 'total-energy',originalCV=cv,originalFace=(cv-1)%nf+1,layer=(cv-1)//nf+1))
out=dict(status='diagnosed-complete-raw-Schur-scalar-zero-rows-columns-and-scale-spread',processID=os.getpid(),cellCount=n,scalarUnknowns=N,
 zeroScalarRows=len(zeroRows),zeroScalarColumns=len(zeroCols),zeroRowExamples=examples,zeroScalarDiagonalEntries=int(np.count_nonzero(np.r_[dd[:,0,0],dd[:,1,1]]==0)),
 fullStructuralRankProvedByDisjointCellPerfectMatchings=proof,provenStructuralRank=N if proof else None,numericalFullRankClaimed=False,
 rowInfinityNormSpan=span(rowMax),columnInfinityNormSpan=span(columnMax),rowNonzeroCountRange=[int(rnnz.min()),int(rnnz.max())],columnNonzeroCountRange=[int(cnnz.min()),int(cnnz.max())],
 nonzeroScalarOffDiagonalEntries=offdiag,potentialInputOffDiagonalEntriesBelowPoint01InitialColumnNormThreshold=potentialSmall,
 potentialFraction=potentialSmall/max(1,offdiag),actualILUTPDroppedEntriesOrZeroPivotCauseProved=False,
 interpretation='Initial column-norm threshold count is only a diagnostic; updated ILUTP coefficients and mandatory U diagonal handling differ. Structural rank does not prove numerical rank. No threshold or physical equation changed.',
 newSelectedJApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
