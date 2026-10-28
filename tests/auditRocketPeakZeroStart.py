# -*- coding: utf-8 -*-
"""独立核对固定同一增广目标的零初猜：稠密解、完整C像及当前真实状态幅度。"""
from pathlib import Path
import os,json,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat
from screenRocketFullFLeastSquares import ROOT,VAL,sha,local
from screenRocketParentImageScale import setup

mode,case=sys.argv[1:3];assert mode in ('unit','real') and case.replace('_','').isalnum()
folder=VAL/case;target=VAL/(case+'_independent_zero_start.json');assert not target.exists();bindings={}


def read(path,variables=None):
    path=Path(path);path=path if path.is_absolute() else ROOT/path;bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


def rel(a,b):return float(np.linalg.norm(a-b)/max(1,np.linalg.norm(b)))


checks=[]
if mode=='unit':
    s=read(folder/'peak_zero_unit_evidence.mat')
    for a in s['cases']:
        b=a['source'];base=b['base'];n=len(base['Q']);N=n*5
        zr,wr=a['zeroReport'],a['warmReport'];assert zr['penaltyInOriginalResidualUnits']==wr['penaltyInOriginalResidualUnits']
        assert zr['metricNormalization']==wr['metricNormalization'] and zr['zeroDirectionNormalizedObjective']==1
        A=np.row_stack([base['C'],zr['penaltyInOriginalResidualUnits']*b['W']]);rhs=np.r_[base['rhs'].ravel(order='F'),np.zeros(N)]
        dense=np.linalg.lstsq(A,rhs,rcond=None)[0];zero_error=rel(a['zeroValue'].ravel(order='F'),dense);warm_error=rel(a['warmValue'].ravel(order='F'),dense)
        assert zero_error<1e-8 and warm_error<1e-8
        lift=base['M']@(base['columnScale']*a['zeroCoordinates']).ravel(order='F');assert rel(lift,a['zeroValue'].ravel(order='F'))<1e-12
        raw=(base['W']@a['zeroValue'].ravel(order='F')).reshape(n,5,order='F')
        assert rel(raw,a['unweightedFinal'])<1e-12
        np.testing.assert_allclose(zr['maximumDensityFraction'],np.max(abs(raw[:,0])),rtol=1e-12)
        assert zr['primitiveUpdateNormIncludesRowWeights'] and zr['beatsZeroDirectionObjective']
        checks.append(dict(zeroDenseDifference=zero_error,warmDenseDifference=warm_error,rawStateMetricDifference=rel(raw,a['unweightedFinal']),samePenalty=True))
else:
    s=read(folder/'peak_zero_candidate_evidence.mat');model=read(folder/'peak_zero_model_evidence.mat');detail=s['detail'];old=s['previousDetail']
    assert detail['penaltyInOriginalResidualUnits']==old['penaltyInOriginalResidualUnits'] and detail['metricNormalization']==old['metricNormalization']
    assert detail['zeroDirectionNormalizedObjective']==1
    if 'initialEvidenceFile' in s:
        previous=read(s['initialEvidenceFile'],['value','detail','sourceBinding']);np.testing.assert_array_equal(s['initial'],previous['value'])
        assert previous['sourceBinding']==s['sourceBinding']
        np.testing.assert_allclose(detail['initialNormalizedAugmentedObjective'],previous['detail']['normalizedAugmentedObjective'],rtol=1e-12)
        assert not model['report']['unrestartedLSQRClaimed']
    else:assert np.all(s['initial']==0) and detail['initialRelativeCResidual']==1
    packed=read(s['packedFile'],['packed'])['packed'];factor=read(s['factorFile'],['factor'])['factor'];A,plan,order=setup(packed,factor);del packed,factor
    scale=read(s['scaleFile'],['columnScale'])['columnScale'];n=len(order)
    lift=s['initial'].copy();lift[order]+=local(plan,(s['coordinates'][order]*scale[order]).ravel()).reshape(n,5)
    lift_error=rel(lift,s['value']);assert lift_error<1e-8
    image=np.empty_like(lift);image[order]=(A@s['value'][order].ravel()).reshape(n,5)
    relative_c=float(np.linalg.norm(s['rhs']-image)/np.linalg.norm(s['rhs']));np.testing.assert_allclose(relative_c,detail['relativeCResidual'],rtol=1e-12)
    cache=read(s['cacheFile'],['Q','scale','gas']);q=cache['Q'];gas=cache['gas'];rho=q[:,0];u=q[:,1:4]/rho[:,None];cv=gas['R']/(gas['gamma']-1)
    T=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv);dq=s['value']*cache['scale'];dr=dq[:,0]
    du=(dq[:,1:4]-u*dr[:,None])/rho[:,None];dt=(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*dr)/(rho*cv)
    raw=np.column_stack([dr/rho,du/np.sqrt(gas['gamma']*gas['R']*T)[:,None],dt/T]);raw_error=rel(raw,s['unweightedFinal']);assert raw_error<1e-12
    weighted=raw*s['metric']['peakRowWeights'];normed=np.linalg.norm(weighted);np.testing.assert_allclose(normed,detail['primitiveUpdateNorm'],rtol=1e-12)
    for key,expected in [('maximumDensityFraction',np.max(abs(raw[:,0]))),('maximumVelocityOverSound',np.max(np.linalg.norm(raw[:,1:4],axis=1))),('maximumTemperatureFraction',np.max(abs(raw[:,4])))]:
        np.testing.assert_allclose(detail[key],expected,rtol=1e-12)
    denom=np.linalg.norm(s['rhs']);objective=relative_c**2+(detail['penaltyInOriginalResidualUnits']*normed/denom)**2
    dc=abs(detail['relativeCResidual'])*1e-12;dw=abs(detail['primitiveUpdateNorm'])*1e-12
    bound=2*abs(detail['relativeCResidual'])*dc+dc*dc+(detail['penaltyInOriginalResidualUnits']/denom)**2*(2*detail['primitiveUpdateNorm']*dw+dw*dw)+16*np.finfo(float).eps*detail['normalizedAugmentedObjective']
    assert abs(objective-detail['normalizedAugmentedObjective'])<=bound
    np.testing.assert_array_equal(model['value'],s['value']);full=float(np.linalg.norm(s['rhs']+model['fullModel'])/denom)
    np.testing.assert_allclose(full,model['report']['completeSelectedModelRelativeDefect'],rtol=1e-12)
    assert model['report']['modelPointOneMet']==(full<=.1) and model['report']['newActualNSEvaluations']==0
    checks.append(dict(originalCVs=n,coordinateLiftDifference=lift_error,rawStateMetricDifference=raw_error,relativeCResidual=relative_c,
                       normalizedAugmentedObjective=objective,completeModelRelativeDefect=full,samePenaltyAsWarmRun=True,
                       candidateBeatsKnownZeroObjective=bool(detail['beatsZeroDirectionObjective']),allRealBidiagonalizationColumnsIndependentlyVerified=False))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-fixed-objective-zero-start-and-raw-state-metric',mode=mode,case=case,processID=os.getpid(),checks=checks,
         newProductionSolvedAreaM2=0,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
