# -*- coding: utf-8 -*-
"""核对同初态、同实际NS方向的候选坐标对照，保持并列历史和真实调用量。"""
from pathlib import Path
import os,json,hashlib
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
target=VAL/'primitivereservoir01_independent_same_direction_trial.json';assert not target.exists()
bindings={}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def read(path,variables=None):
    bindings[path.relative_to(ROOT).as_posix()]=sha(path)
    return loadmat(path,simplify_cells=True,variable_names=variables)


a=read(VAL/'krylov22/wall_velocity_krylov_evidence.mat',['result','report','gas'])
b=read(VAL/'krylov23/wall_velocity_krylov_evidence.mat',['result','report','gas'])
rows=lambda x:[x] if isinstance(x,dict) else list(x)
x,y=a['result'],b['result'];ea=rows(x['linearLog'])[0]['steadyDirectionEvidence'];eb=rows(y['linearLog'])[0]['steadyDirectionEvidence']
for key in ['initialQ','initialRate','residualScale']:np.testing.assert_array_equal(x[key],y[key])
for key in ['directionU','value','ratePlus','rateMinusOrBase','actualScalarStep']:np.testing.assert_array_equal(ea[key],eb[key])
assert x['options'].get('trialStateCoordinates','conservative')=='conservative'
assert y['options']['trialStateCoordinates']==y['nonlinearTrialStateCoordinates']=='primitive-tangent'
for key in ['steadyResidualTolerance','steadyFluxTolerance','spatialConservationTolerance','linearRelativeTolerance','sufficientDecrease']:
    assert x['options'][key]==y['options'][key]
gas=b['gas'];cv=gas['R']/(gas['gamma']-1);tr=rows(y['trace']);assert len(tr)==1
t=tr[0];q=t['initialQ'];dq=t['linearUpdateU']*y['residualScale'];rho=q[:,0];u=q[:,1:4]/rho[:,None]
T=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv);dr=dq[:,0]
du=(dq[:,1:4]-u*dr[:,None])/rho[:,None]
dt=(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*dr)/(rho*cv)
rr=rho+t['alpha']*dr;uu=u+t['alpha']*du;tt=T+t['alpha']*dt
expected=np.column_stack([rr,rr[:,None]*uu,rr*(cv*tt+.5*np.sum(uu*uu,axis=1))])
difference=float(np.linalg.norm(expected-y['Q'])/max(1,np.linalg.norm(y['Q'])));assert difference<5e-15
assert np.linalg.norm(expected-(q+t['alpha']*dq))>0
assert not y['steadyConverged'] and not b['report']['physicalSolution'] and b['report']['newProductionSolvedAreaM2']==0
assert a['report']['totalNewEquationAcceptedIterations']==b['report']['totalNewEquationAcceptedIterations']==2
for case in ['krylov22','krylov23']:
    public=VAL/(case+'_independent_riemann_reservoir.json');data=json.loads(public.read_text(encoding='utf8'));bindings[public.relative_to(ROOT).as_posix()]=sha(public)
    assert data['status'].startswith('passed') and data['checks']['continuation']['predecessorCase']=='krylov20'
summary=VAL/'primitivereservoir01/primitive_reservoir_trial_summary.json';bindings[summary.relative_to(ROOT).as_posix()]=sha(summary)
r=json.loads(summary.read_text(encoding='utf8'));assert r['sameBaseSameDirectionComparison'] and not r['alternateHistoriesMerged']
assert r['acceptedAlpha']==t['alpha'] and r['newMainActualNS']==y['residualEvaluations'] and r['newStandaloneProbeNS']==0
assert r['main']['maximumResidual']==y['maximumResidual'] and r['main']['rmsMerit']==y['rmsMerit']
for key,stat in [('rmsBetterThanReference','rmsMerit'),('maximumBetterThanReference','maximumResidual')]:assert r[key]==bool(y[stat]<x[stat])
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-same-base-same-original-NS-direction-and-primitive-candidate-branch',processID=os.getpid(),
         initialStateAndRateExactlyIdentical=True,originalNSDirectionAndJvpExactlyIdentical=True,
         originalAcceptanceGatesIdentical=True,primitiveCandidateRelativeDifference=difference,
         conservativeAcceptedAlpha=float(rows(x['trace'])[0]['alpha']),primitiveAcceptedAlpha=float(t['alpha']),
         conservativeRMS=float(x['rmsMerit']),primitiveRMS=float(y['rmsMerit']),
         conservativeMaximum=float(x['maximumResidual']),primitiveMaximum=float(y['maximumResidual']),
         alternateHistoriesMerged=False,newStandaloneProbeNS=0,newMainNS=int(y['residualEvaluations']),
         newProductionSolvedAreaM2=0,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))
