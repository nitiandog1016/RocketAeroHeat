# -*- coding: utf-8 -*-
"""独立验证R*A与P/R的完整小域映射，同时核对两种残差的真实身份。"""
from pathlib import Path
import os,json,hashlib,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';name=sys.argv[1];assert name.replace('_','').isalnum();target=VAL/(name+'_independent_row_momentum.json');assert not target.exists()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as s:
  for b in iter(lambda:s.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def rows(v):return [v] if isinstance(v,dict) else list(np.atleast_1d(v).ravel())
def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))
file=VAL/name/'row_momentum_unit_evidence.mat';s=loadmat(file,simplify_cells=True);results=[]
for case in rows(s['cases']):
 A=case['A'];P=case['P'];R=1/np.max(abs(A),axis=1);B=R[:,None]*A;PR=P/R[None,:];b=case['rhs'].ravel(order='F');x=case['value'].ravel(order='F');d=case['detail']
 np.testing.assert_allclose(case['rowScale'].ravel(order='F'),R,rtol=1e-13,atol=0);assert rel(B,case['B'])<1e-13 and rel(PR,case['PRinverse'])<1e-13
 actual=b-A@x;unscaled=float(np.linalg.norm(actual)/np.linalg.norm(b));scaled=float(np.linalg.norm(R*actual)/np.linalg.norm(R*b));assert unscaled<1e-8 and rel(x,np.linalg.solve(A,b))<1e-8
 assert abs(unscaled-d['relativeAResidual'])<1e-10 and abs(scaled-d['relativeRowScaledResidual'])<1e-10
 assert bool(d['originalUnitTargetMet'])==(d['relativeAResidual']<=d['options']['relativeTolerance']) and bool(d['rowScaledTargetMet'])==(d['krylov']['flag']==0)
 cycle=rows(d['krylov']['cycles'])[0];V=cycle['V'];Z=cycle['Z'];H=cycle['smallH'];scales=cycle['imageNormalizationScales'];mapping=rel((PR@V[:,:Z.shape[1]])/scales[None,:],Z);relation=rel(B@Z,V@H)
 assert mapping<1e-10 and relation<1e-8
 combination=cycle['xStart']+Z@cycle['weights'];assert rel(combination,x)<1e-10
 results.append(dict(originalAResidual=unscaled,rowScaledResidual=scaled,completePRinverseMappingDifference=mapping,completeRowSystemArnoldiDifference=relation,steps=int(d['krylov']['arnoldiSteps'])))
assert s['report']['budgetRejected'] and s['report']['zeroRowRejected']
out=dict(status='passed-independent-two-orderings-equivalent-row-system-PRinverse-and-original-unit-residuals',processID=os.getpid(),cases=results,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,inputSHA256={str(file.relative_to(ROOT)):sha(file)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
