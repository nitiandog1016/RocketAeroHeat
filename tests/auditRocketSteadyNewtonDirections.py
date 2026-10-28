# -*- coding: utf-8 -*-
"""核对每个直接稳态方向的原NS差分、(-J)缺陷、下降界与所有候选，包含失败方向。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1];kind=sys.argv[2] if len(sys.argv)>2 else 'rocket'
variant=sys.argv[3] if len(sys.argv)>3 else 'result'
revision=sys.argv[4] if len(sys.argv)>4 else '';assert not revision or revision.isdigit()
assert case.replace('_','').isalnum() and kind in ('unit','rocket') and variant in ('result','limited','ready')
assert kind=='unit' or variant=='result'
target=VAL/(case+('_'+variant if variant!='result' else '')+'_independent_newton_directions'+revision+'.json');assert not target.exists()
file=VAL/case/('steady_newton_unit_evidence.mat' if kind=='unit' else 'wall_velocity_krylov_evidence.mat')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
bindings={str(file.relative_to(ROOT)):sha(file)};s=loadmat(file,simplify_cells=True);r=s[variant];scale=r['residualScale'];gas=s['gas']
assert r['linearEquation']=='steady-newton' and r['options']['linearEquation']=='steady-newton' and r['preconditionerTimeStepOnly']
rows=lambda x:[x] if isinstance(x,dict) else list(x)
trace=rows(r['trace']);logs=rows(r['linearLog']);trials=rows(r['lineLog']);cv=gas['R']/(gas['gamma']-1)
def base(e):
    index=int(e['baseAcceptedIteration'])
    assert 0<=index<=len(trace)
    return (r['initialQ'],r['initialRate']) if index==0 else (trace[index-1]['Q'],trace[index-1]['rate'])
def cell_state(q):
    if not np.isfinite(q).all() or np.any(q[:,0]<=0):return 'CutNS:State'
    rho=q[:,0];u=q[:,1:4]/rho[:,None];internal=q[:,4]-.5*rho*np.sum(u*u,axis=1)
    if np.any(internal<=0):return 'CutNS:InternalEnergy'
    T=internal/(rho*cv)
    return 'CutNS:Temperature' if np.any((T<gas['temperatureRange'][0])|(T>gas['temperatureRange'][1])) else ''
relative=lambda a,b:float(np.max(abs(a-b))/max(1,np.max(abs(b))))
linear_checks=[]
for a in logs:
    e=a['steadyDirectionEvidence'];Q,rate=base(e);rhs=rate/scale
    if e['scheme']=='zero-direction':
        assert np.all(e['directionU']==0) and np.all(e['value']==0) and e['actualOriginalResidualEvaluations']==0
        rebuilt=np.zeros_like(rhs)
    else:
        h=e['actualScalarStep'];assert h>0 and not cell_state(Q+h*e['directionU']*scale)
        if e['scheme']=='forward':np.testing.assert_array_equal(e['rateMinusOrBase'],rate);factor=1
        else:assert e['scheme']=='central' and not cell_state(Q-h*e['directionU']*scale);factor=2
        rebuilt=(e['ratePlus']-e['rateMinusOrBase'])/scale/(factor*h)
    error=relative(rebuilt,e['value']);assert error<3e-15
    eta=float(np.linalg.norm(-rebuilt-rhs)/np.linalg.norm(rhs));slope=float(np.sum(rhs*rebuilt)/np.sum(rhs*rhs))
    assert abs(eta-a['fullSystemRelativeDefect'])<3e-12 and abs(slope-a['normalizedMeritDerivative'])<3e-12
    assert abs(a['inexactNewtonUpperBound']-(-1+eta))<3e-12 and slope<=-1+eta+3e-12
    descending=slope<0 and slope<=-1+eta+3e-12
    assert a['meritDescentVerified']==descending and a['fullSystemTargetMet']==(eta<=a['target'] and descending)
    assert a['target']==r['options']['linearRelativeTolerance'] and a['preconditionerTimeStepOnly']
    linear_checks.append(dict(iteration=int(e['iteration']),retry=int(e['retry']),actualSteadyLinearDefect=eta,
                              normalizedMeritDerivative=slope,inexactNewtonUpperBound=-1+eta,
                              targetMet=bool(a['fullSystemTargetMet']),originalJvpArithmeticRelativeDifference=error))
candidate_calls=0
coordinates=r['options'].get('trialStateCoordinates','conservative')
assert coordinates in ('conservative','primitive-tangent')
def trial_state(Q,direction,alpha):
    if coordinates=='conservative':return Q+alpha*direction*scale
    # 独立链式回代原始变量切向路径；原NS Jv检查仍使用守恒扰动，不能混用。
    rho=Q[:,0];u=Q[:,1:4]/rho[:,None];T=(Q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
    dq=direction*scale;dr=dq[:,0];du=(dq[:,1:4]-u*dr[:,None])/rho[:,None]
    dt=(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*dr)/(rho*cv)
    new_rho=rho+alpha*dr;new_u=u+alpha*du;new_t=T+alpha*dt
    return np.column_stack([new_rho,new_rho[:,None]*new_u,new_rho*(cv*new_t+.5*np.sum(new_u*new_u,axis=1))])
for row in trials:
    matches=[a for a in logs if a['steadyDirectionEvidence']['iteration']==row['iteration'] and a['steadyDirectionEvidence']['retry']==row['retry']]
    assert len(matches)==1 and matches[0]['fullSystemTargetMet'];e=matches[0]['steadyDirectionEvidence'];Q,rate=base(e)
    q=trial_state(Q,e['directionU'],row['alpha']);state_error=cell_state(q);actual_error=row['errorIdentifier'] if isinstance(row['errorIdentifier'],str) else ''
    if state_error:assert state_error==actual_error
    else:candidate_calls+=1
    if not actual_error:
        sufficient=row['candidateRmsMerit']<=(1-r['options']['sufficientDecrease']*row['alpha'])*row['baseRmsMerit']
        assert row['accepted']==sufficient
    if row['accepted']:
        t=trace[int(row['iteration'])-1]
        if coordinates=='conservative':np.testing.assert_array_equal(q,t['Q'])
        else:
            assert t['trialStateCoordinates']==r['nonlinearTrialStateCoordinates']==coordinates
            assert relative(q,t['Q'])<5e-15
        np.testing.assert_array_equal(e['directionU'],t['linearUpdateU'])
        assert t['linearReport']['meritDescentVerified'] and 'steadyDirectionEvidence' not in t['linearReport']
jvp_rows=rows(r['jvpLog']);jvp_calls=sum(x['actualOriginalResidualEvaluations'] for x in jvp_rows)
assert jvp_calls==r['jacobianResidualEvaluations'] and 1+jvp_calls+candidate_calls==r['residualEvaluations']
assert len(trials)==r['candidateTrials'] and sum(x['accepted'] for x in trials)==r['acceptedNonlinearIterations']
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-all-steady-Newton-direction-evidence-descent-and-candidate-work-audit',case=case,kind=kind,variant=variant,processID=os.getpid(),
         linearDirections=linear_checks,allDirectionsIncludingFailuresSaved=True,candidateTrials=len(trials),candidateNSEvaluations=candidate_calls,
         jacobianNSEvaluations=int(jvp_calls),totalSolverNSEvaluations=int(r['residualEvaluations']),fullIndependentNSSolver=False,
         nonlinearTrialCoordinates=coordinates,physicalSolution=False,newProductionSolvedAreaM2=0,inputSHA256=bindings,scriptSHA256=sha(__file__))
if r['options'].get('jacobianAction')=='selected-analytic':
    models=[x for x in jvp_rows if x['scheme']=='selected-analytic'];checks=[x for x in jvp_rows if x['scheme']!='selected-analytic']
    assert len(models)==r['selectedAnalyticApplications']>0 and len(jvp_rows)==r['matrixVectorProducts']
    assert r['originalNSCombinedDirectionCheckPreserved'] and all(a['actualCombinedDirectionCheckedWithOriginalNS'] for a in logs)
    for row in models:
        assert row['actualOriginalResidualEvaluations']==0 and np.isnan(row['actualScalarStep'])
        assert row['boundaryStateSensitivityIncluded'] and row['fullWLSViscousAndWorkChainIncluded'] and not row['nonlinearAcceptanceClaimed']
    assert len(checks)==len(logs) and all(x['scheme'] in ('forward','central') for x in checks)
    out['selectedAnalyticApplications']=len(models);out['originalCombinedDirectionChecks']=len(checks)
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
