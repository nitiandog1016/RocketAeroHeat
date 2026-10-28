# -*- coding: utf-8 -*-
"""读取原门槛拒绝的基，定位正交丢失与小矩阵尺度；零新模型/NS。"""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/case/'rejected_basis_numerics.json';assert not target.exists();file=VAL/case/'arnoldi_failure_evidence.mat';summary=VAL/case/'block_ilu_failure_diagnostic_summary.json'
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()
f=loadmat(file,simplify_cells=True)['failure'];V=f['V'];H=f['smallH'];G=V.T@V;errors=[float(np.linalg.norm(G[:j,:j]-np.eye(j))) for j in range(1,G.shape[0]+1)]
singular=np.linalg.svd(H,compute_uv=False);column=np.linalg.norm(H,axis=0);positive=column[column>0];orth=float(np.linalg.norm(G-np.eye(len(G))))
np.testing.assert_allclose(orth,f['orthogonalityError'],rtol=2e-10,atol=1e-11)
weights=np.atleast_1d(f['weights']);y=H@weights;predicted=f['rStart']-V@y
recompute=float(np.linalg.norm(predicted-f['predictedRemaining']));rhs_norm=float(np.linalg.norm(f['b']))
# 大项相消时，跨库复算不保证相对小剩余的精度；只审核形成过程的舍入界。
# 原RHS阈值失败仍明确记录，不把此形成误差界用于放行原Arnoldi/物理方向。
formation=float(np.linalg.norm(f['rStart'])+np.linalg.norm(abs(V)@abs(y))+np.linalg.norm(abs(V)@(abs(H)@abs(weights))))
eps=np.finfo(float).eps;gamma=lambda k:k*eps/(1-k*eps)
formation_bound=4*(gamma(H.shape[1])+gamma(V.shape[1])+eps)*formation
assert recompute<=formation_bound
assert f['modelRelationThreshold']==f['orthogonalityThreshold']==1e-8 and not f['usableDirectionClaimed']
out=dict(status='source-bound-rejected-real-ILU-basis-orthogonality-and-H-scaling-diagnostic',processID=os.getpid(),
    orthogonalityError=orth,originalOrthogonalityThreshold=1e-8,prefixBasisOrthogonalityErrors=errors,
    firstBasisColumnCountOverOriginalGate=next((j+1 for j,x in enumerate(errors) if x>=1e-8),None),
    basisColumnNormRange=[float(np.sqrt(np.diag(G)).min()),float(np.sqrt(np.diag(G)).max())],
    HColumnNormRange=[float(column.min()),float(column.max())],HPositiveColumnNormRatio=float(positive.max()/positive.min()),
    HSingularValueRange=[float(singular.min()),float(singular.max())],HNumericalRankAtExistingSmallSVDThreshold=int(np.count_nonzero(singular>1e-13*singular.max())),
    savedActualResidualRelativeToRhs=float(np.linalg.norm(f['actualRemaining'])/np.linalg.norm(f['b'])),
    savedProjectedResidualRelativeToRhs=float(np.linalg.norm(f['predictedRemaining'])/np.linalg.norm(f['b'])),
    modelRelationRelativeDifference=float(f['modelRelationRelativeDifference']),
    projectedRecomputationDifference=recompute,projectedRecomputationRelativeToRhs=recompute/max(1,rhs_norm),
    originalRecomputationRhsThreshold=1e-10,originalRecomputationRhsCheckPassed=bool(recompute/max(1,rhs_norm)<1e-10),
    formationNormBound=formation,projectedRecomputationRoundoffBound=formation_bound,
    projectedRecomputationWithinFormationBound=True,physicalOrArnoldiAcceptanceFromFormationBound=False,
    proposedNumericalChecks='Whole-direction image normalization and adaptive extra reorthogonalization, with the original final gates unchanged; default old output must replay. This does not prove the ILU0 inverse will become useful.',
    newModelApplications=0,newActualNSEvaluations=0,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(p.relative_to(ROOT)):sha(p) for p in (file,summary)},scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps(out))
