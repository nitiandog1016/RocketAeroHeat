# -*- coding: utf-8 -*-
"""独立核对局部粗空间小域的新Q重建、原NS方向及候选拒绝/缩步，不复写证据。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation';case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_targeted_integration.json');assert not target.exists()
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
file=VAL/case/'steady_newton_unit_evidence.mat';digest=sha(file);s=loadmat(file,simplify_cells=True);r=s['rebuilt'];gas=s['gas'];scale=r['residualScale']
rows=lambda x:[x] if isinstance(x,dict) else list(x)
assert r['coarseOperatorReuses']==1 and r['coarseOperatorBuilds']==1 and r['analyticCacheBuilds']==2 and r['coarseConstructionModelApplications']==140
assert r['acceptedNonlinearIterations']==2 and r['coarseBasis']=='conservative-targeted-spatial'
assert r['coarsePreconditionerModelApplications']==r['preconditionerApplications'] and r['coarseReuseVerificationModelApplications']==1
trace=rows(r['trace']);logs=rows(r['linearLog']);trials=rows(r['lineLog']);calls=1;invalid=0;directions=[]
for a in logs:
    e=a['steadyDirectionEvidence'];index=int(e['baseAcceptedIteration']);Q=r['initialQ'] if index==0 else trace[index-1]['Q'];rate=r['initialRate'] if index==0 else trace[index-1]['rate']
    h=e['actualScalarStep'];derivative=(e['ratePlus']-e['rateMinusOrBase'])/scale/h;np.testing.assert_array_equal(rate,e['rateMinusOrBase']);np.testing.assert_array_equal(derivative,e['value'])
    rhs=rate/scale;eta=np.linalg.norm(rhs+derivative)/np.linalg.norm(rhs);assert abs(eta-a['fullSystemRelativeDefect'])<3e-12 and eta<=.1
    slope=np.sum(rhs*derivative)/np.sum(rhs*rhs);assert slope<0 and slope<=-1+eta+3e-12;calls+=int(e['actualOriginalResidualEvaluations']);directions.append(float(eta))
    for trial in trials:
        if trial['iteration']!=e['iteration'] or trial['retry']!=e['retry']:continue
        q=Q+trial['alpha']*e['directionU']*scale
        error=trial['errorIdentifier'] if isinstance(trial['errorIdentifier'],str) else ''
        if error=='CutNS:State':assert np.any(q[:,0]<=0) or not np.isfinite(q).all();assert not trial['accepted'];invalid+=1;continue
        assert not error;calls+=1;rho=q[:,0];u=q[:,1:4]/rho[:,None];T=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*gas['R']/(gas['gamma']-1))
        assert np.all((T>=gas['temperatureRange'][0])&(T<=gas['temperatureRange'][1]))
        assert trial['accepted']==(trial['candidateRmsMerit']<=(1-r['options']['sufficientDecrease']*trial['alpha'])*trial['baseRmsMerit'])
        if trial['accepted']:np.testing.assert_array_equal(q,trace[int(trial['iteration'])-1]['Q'])
assert calls==r['residualEvaluations']==5 and invalid==8
assert s['result']['candidateTrials']==3 and s['result']['acceptedNonlinearIterations']==1 and s['limited']['candidateTrials']==0
assert sha(file)==digest
out=dict(status='passed-independent-targeted-new-Q-rebuild-original-direction-and-state-rejection-accounting',case=case,processID=os.getpid(),
    actualOriginalNSEvaluations=int(calls),acceptedIterations=2,invalidDensityCandidatesRejected=invalid,actualDirectionDefects=directions,
    newQCoarseBuilds=1,reusedInitialCoarseOperators=1,newCoarseConstructionModelApplications=140,physicalSolution=False,newProductionSolvedAreaM2=0,
    inputSHA256={str(file.relative_to(ROOT)):digest},scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
