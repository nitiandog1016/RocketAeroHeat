# -*- coding: utf-8 -*-
"""用标准库80位小矩阵残差，独立验证原单位逆精化；不重新装配真实域。"""
from pathlib import Path
from decimal import Decimal,localcontext
import os,json,hashlib
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat,savemat
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';folder=VAL/'momentumline_precision01';assert not folder.exists();folder.mkdir()
source=VAL/'momentumlinereal01/line_inverse_failure.mat';s=loadmat(source,simplify_cells=True);A=s['originalLine'];X=s['candidateInverse'].copy();N=len(A);L=s['L'];U=s['U'];R=s['R'];C=s['C'];order=s['permutation'].astype(int)-1;trials=[];candidates=[];residuals=[]
with localcontext() as ctx:
 ctx.prec=80;AA=[[Decimal.from_float(float(x)) for x in row] for row in A]
 for step in range(3):
  XX=[[Decimal.from_float(float(x)) for x in row] for row in X]
  residual=np.array([[float((Decimal(1) if i==j else Decimal(0))-sum((AA[i][q]*XX[q][j] for q in range(N)),Decimal(0))) for j in range(N)] for i in range(N)])
  candidates.append(X.copy());residuals.append(residual);trials.append(dict(refinementSteps=step,decimal80InverseDefect=float(np.max(abs(residual))),ordinaryBLASDefect=float(np.max(abs(A@X-np.eye(N))))))
  right=R[:,None]*residual;X+=C[:,None]*solve_triangular(U,solve_triangular(L,right[order],lower=True,unit_diagonal=True),lower=False)
binding=dict(sourceFile=str(source),sourceSHA256=hashlib.sha256(source.read_bytes()).hexdigest(),face=int(s['face']))
savemat(folder/'precision_reference.mat',dict(A=A,candidates=np.stack(candidates,axis=2),referenceResiduals=np.stack(residuals,axis=2),binding=binding))
out=dict(status='independent-decimal80-original-double-input-inverse-residual-reference',processID=os.getpid(),binding=binding,trials=trials,originalThreshold=1e-8,newActualNSEvaluations=0,newProductionSolvedAreaM2=0,scriptSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(folder/'precision_reference.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
