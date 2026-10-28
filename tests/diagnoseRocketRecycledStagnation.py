# -*- coding: utf-8 -*-
"""读取已存周期和图对角，定位停滞及各方程时间正则占比；不运行NS或新模型。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/case/'recycled_stagnation_diagnostic.json';assert not target.exists()
source=VAL/case/'wall_velocity_krylov_evidence.mat';r=loadmat(source,simplify_cells=True,variable_names=['result'])['result']
rows=lambda x:[x] if isinstance(x,dict) else list(np.atleast_1d(x))
a=rows(r['linearLog'])[-1];g=a['initialGuessEvidence'];rr=a['recycledKrylovReport'];rhs=r['initialRate']/r['residualScale'];remaining=rhs+a['steadyDirectionEvidence']['value']
initial=g['correctionRhs'];energy=np.sum(remaining**2,axis=1);total=float(energy.sum());D=r['firstOperator']['diagonal'];temporal=1/r['firstOperator']['pseudoTimeStepSeconds']
names=['density','xMomentum','yMomentum','zMomentum','totalEnergy'];diagonal=[]
for i,name in enumerate(names):
    ratio=np.abs(D[i,i]-temporal)/temporal;weight=remaining[:,i]**2;den=float(weight.sum());order=np.argsort(ratio);cumulative=np.cumsum(weight[order])/max(den,np.finfo(float).tiny)
    median=float(ratio[order[min(np.searchsorted(cumulative,.5),len(order)-1)]])
    diagonal.append(dict(equation=name,scaledRemainingSquaredNormFraction=den/total,
        remainingWeightedMedianPhysicsToTemporalDiagonalRatio=median,
        remainingWeightedMeanTemporalShareOfAbsoluteDiagonalSplit=float(np.sum(weight/(1+ratio))/max(den,np.finfo(float).tiny)),
        equationRemainingFractionAtPhysicsToTemporalDiagonalRatioBelowPointOne=float(weight[ratio<.1].sum()/max(den,np.finfo(float).tiny))))
cycles=[]
for c in rows(rr['cycles']):
    H=np.asarray(c['smallH']);H=H.reshape((int(c['arnoldiSteps'])+1,int(c['arnoldiSteps'])),order='F');singular=np.linalg.svd(H,compute_uv=False)
    cycles.append(dict(cycle=int(c['cycle']),arnoldiSteps=int(c['arnoldiSteps']),retainedBefore=int(c['retainedBefore']),retainedAfter=int(c['retainedAfter']),
        startResidualNorm=float(c['startResidualNorm']),actualResidualNorm=float(c['actualResidualNorm']),
        projectedHSingularRange=[float(singular.min()),float(singular.max())],projectedHSingularCondition=float(singular.max()/singular.min()),
        arnoldiRelationRelativeDifference=float(c['arnoldiRelationRelativeDifference']),basisOrthogonalityError=float(c['basisOrthogonalityError'])))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

out=dict(status='source-bound-recycled-stagnation-and-approximate-graph-diagonal-diagnostic',processID=os.getpid(),
    initialCurrentModelDefect=float(np.linalg.norm(initial)/np.linalg.norm(rhs)),actualFinalTotalNSDefect=float(a['fullSystemRelativeDefect']),
    remainingSquaredNormFractionRemoved=float(1-np.linalg.norm(remaining)**2/np.linalg.norm(initial)**2),cycles=cycles,diagonalSplit=diagonal,
    diagonalInterpretation='Approximate preconditioner B diagonal minus I/tau only; not full-J diagonal, global conditioning, physics time or proof of a better parameter.',
    newModelApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(source.relative_to(ROOT)):sha(source)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
