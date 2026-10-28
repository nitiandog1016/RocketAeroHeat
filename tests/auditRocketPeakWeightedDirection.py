# -*- coding: utf-8 -*-
"""独立复核冻结局部峰值权重、小域稠密解及当前真实状态全部坐标/完整父C像。"""
from pathlib import Path
import os,json,sys,math
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha,local
from screenRocketParentImageScale import setup

mode,case=sys.argv[1:3];assert mode in ('unit','real') and case.replace('_','').isalnum()
folder=VAL/case;target=VAL/(case+'_independent_peak_direction.json');assert not target.exists();bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))


def primitive_metric(Q,scale,gas,v):
    rho=Q[:,0];u=Q[:,1:4]/rho[:,None];cv=gas['R']/(gas['gamma']-1)
    T=(Q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv);dq=v*scale;dr=dq[:,0]
    du=(dq[:,1:4]-u*dr[:,None])/rho[:,None]
    dt=(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*dr)/(rho*cv)
    return np.column_stack([dr/rho,du/np.sqrt(gas['gamma']*gas['R']*T)[:,None],dt/T])


checks=[]
if mode=='unit':
    s=read(folder/'peak_weighted_unit_evidence.mat')
    for a in s['cases']:
        b=a['base'];options=a['weighting']['options'];n=len(b['Q']);N=5*n
        z=primitive_metric(b['Q'],b['scale'],b['gas'],a['reference'])
        weights=np.minimum(options['maximumRowWeight'],np.hypot(1,z/options['transitionFraction']))
        np.testing.assert_allclose(weights,a['metric']['peakRowWeights'],rtol=1e-12)
        W=weights.ravel(order='F')[:,None]*b['W'];matrix_error=rel(W,a['W']);assert matrix_error<1e-12
        A=np.row_stack([b['C'],a['detail']['penaltyInOriginalResidualUnits']*W]);rhs=np.r_[b['rhs'].ravel(order='F'),np.zeros(N)]
        dense=np.linalg.lstsq(A,rhs,rcond=None)[0];dense_error=rel(a['value'].ravel(order='F'),dense);assert dense_error<1e-8
        value=b['initial'].ravel(order='F')+b['M']@(b['columnScale']*a['coordinates']).ravel(order='F')
        lift_error=rel(a['value'].ravel(order='F'),value);assert lift_error<1e-12
        assert a['weighting']['fixedWithinLinearSolve'] and not a['weighting']['hardStateUpdateBoundEnforced']
        checks.append(dict(allWeightedMatrixDifference=matrix_error,denseAugmentedSolutionDifference=dense_error,
                           allCoordinateLiftDifference=lift_error,adjointDifference=float(a['detail']['adjointIdentityRelativeDifference'])))
else:
    s=read(folder/'current_peak_candidate_evidence.mat');binding=s['sourceBinding']
    assert sha(binding['sourceFile'])==binding['sourceFileSHA256']
    source=read(binding['sourceFile'],['result'])['result'];cache=read(s['cacheFile'],['Q','rate','scale','gas'])
    np.testing.assert_array_equal(source['Q'],cache['Q']);np.testing.assert_array_equal(source['rate'],cache['rate'])
    np.testing.assert_array_equal(s['rhs'],cache['rate']/cache['scale']);del source
    seed=read(s['seedFile'],['value']);np.testing.assert_array_equal(seed['value'],s['initial'])
    packed=read(s['packedFile'],['packed'])['packed'];factor=read(s['factorFile'],['factor'])['factor']
    A,plan,order=setup(packed,factor);del packed,factor
    scale=read(s['scaleFile'],['columnScale'])['columnScale'];n=len(order)
    correction=local(plan,(s['coordinates'][order]*scale[order]).ravel()).reshape(n,5);lift=s['initial'].copy();lift[order]+=correction
    lift_error=rel(lift,s['value']);assert lift_error<1e-8
    image=np.empty_like(lift);image[order]=(A@s['value'][order].ravel()).reshape(n,5)
    c_residual=float(np.linalg.norm(s['rhs']-image)/np.linalg.norm(s['rhs']))
    np.testing.assert_allclose(c_residual,s['detail']['relativeCResidual'],rtol=1e-12)
    raw_initial=primitive_metric(cache['Q'],cache['scale'],cache['gas'],s['initial'])
    raw_final=primitive_metric(cache['Q'],cache['scale'],cache['gas'],s['value'])
    assert rel(raw_initial,s['unweightedInitial'])<1e-12 and rel(raw_final,s['unweightedFinal'])<1e-12
    options=s['weighting']['options'];weights=np.minimum(options['maximumRowWeight'],np.hypot(1,raw_initial/options['transitionFraction']))
    np.testing.assert_allclose(weights,s['metric']['peakRowWeights'],rtol=1e-12)
    detail=s['detail'];assert detail['primitiveUpdateNormIncludesRowWeights']
    expected_norm=np.linalg.norm(weights*raw_final);np.testing.assert_allclose(expected_norm,detail['primitiveUpdateNorm'],rtol=1e-12)
    # 加权优化目标与未加权物理相对变化分开，不能把权重乘积当真实状态变化。
    for key,expected in [('maximumDensityFraction',np.max(abs(raw_final[:,0]))),
                         ('maximumVelocityOverSound',np.max(np.linalg.norm(raw_final[:,1:4],axis=1))),
                         ('maximumTemperatureFraction',np.max(abs(raw_final[:,4])))]:
        np.testing.assert_allclose(detail[key],expected,rtol=1e-12)
    penalty=detail['penaltyInOriginalResidualUnits'];denom=np.linalg.norm(s['rhs']);objective=c_residual*c_residual+(penalty*expected_norm/denom)**2
    dc=abs(detail['relativeCResidual'])*1e-12;dw=abs(detail['primitiveUpdateNorm'])*1e-12
    bound=2*abs(detail['relativeCResidual'])*dc+dc*dc+(penalty/denom)**2*(2*abs(detail['primitiveUpdateNorm'])*dw+dw*dw)+16*np.finfo(float).eps*detail['normalizedAugmentedObjective']
    assert abs(objective-detail['normalizedAugmentedObjective'])<=bound
    model=read(folder/'current_peak_model_evidence.mat');np.testing.assert_array_equal(model['value'],s['value'])
    full=float(np.linalg.norm(s['rhs']+model['fullModel'])/denom)
    accurate_norm=lambda a:math.sqrt(math.fsum(float(x)*float(x) for x in a.ravel()))
    seed_defect=accurate_norm(s['rhs']+model['seedModel'])/accurate_norm(s['rhs'])
    np.testing.assert_allclose(full,model['report']['completeSelectedModelRelativeDefect'],rtol=1e-12)
    # 百万项正平方和采用实际长度的前向舍入界；保持原C、候选模型和0.1资格门槛。
    count=s['rhs'].size;gamma=count*np.finfo(float).eps/(1-count*np.finfo(float).eps)
    seed_bound=abs(seed_defect)*2*gamma/(1-gamma)
    assert abs(seed_defect-model['report']['currentSeedCompleteModelRelativeDefect'])<=seed_bound
    assert not model['report']['oldNumericMatricesFactorsOrDerivativesReused'] and model['report']['currentStateCompleteNSExactlyReplayed']
    assert model['report']['newActualNSEvaluations']==1 and model['report']['newFlowAcceptances']==0
    checks.append(dict(originalCVs=n,totalCoordinateLiftDifference=lift_error,currentCRelativeResidual=c_residual,
                       currentCompleteModelRelativeDefect=full,currentSeedCompleteModelRelativeDefect=seed_defect,
                       seedNormReferenceAbsoluteBound=seed_bound,
                       weightedPrimitiveMetricNorm=float(expected_norm),maximumRawRelativeComponents=np.max(abs(raw_final),axis=0).tolist(),
                       hardStateUpdateBoundEnforced=False,allRealBidiagonalizationColumnsIndependentlyVerified=False))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-frozen-local-peak-metric-and-total-direction',mode=mode,case=case,processID=os.getpid(),checks=checks,
         newProductionSolvedAreaM2=0,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
