# -*- coding: utf-8 -*-
"""独立核对真实三维同方向候选坐标、原残差范数及全部五量通量账。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_independent_audit.json');assert not target.exists()
file=VAL/case/'wall_trial_path_evidence.mat';domain_file=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()
bindings={str(p.relative_to(ROOT)):sha(p) for p in [file,domain_file]}
s=loadmat(file,simplify_cells=True);d=loadmat(domain_file,simplify_cells=True,variable_names=['domain'])['domain']
Q=s['baseQ'];scale=s['scale'];dQ=s['directionU']*scale;gas=s['gas'];cv=gas['R']/(gas['gamma']-1)
rho=Q[:,0];u=Q[:,1:4]/rho[:,None];T=(Q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
drho=dQ[:,0];du=(dQ[:,1:4]-u*drho[:,None])/rho[:,None]
dT=(dQ[:,4]-np.sum(u*dQ[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*drho)/(rho*cv)
base_merit=float(np.sqrt(np.mean((s['baseRate']/scale)**2)))
assert abs(base_merit-s['report']['baseMetrics']['rmsMerit'])<3e-15*base_merit
relative=lambda a,b:float(np.max(abs(a-b))/max(1,np.max(abs(b))))
f=d['openFaces'];w=d['walls'];volume=d['cells']['volume'];m=f['cellMinus'].astype(int)-1;p=f['cellPlus'].astype(int)-1
inside=p>=0;owner=w['ownerCell'].astype(int)-1;checks=[]
for c in s['cases']:
    report=c['summary'];alpha=report['alpha'];coordinates=report['coordinates']
    if coordinates=='conservative':expected=Q+alpha*dQ
    else:
        assert coordinates=='primitive-tangent'
        rr=rho+alpha*drho;uu=u+alpha*du;tt=T+alpha*dT
        expected=np.column_stack([rr,rr[:,None]*uu,rr*(cv*tt+.5*np.sum(uu*uu,axis=1))])
    update_error=relative(c['Q'],expected);assert update_error<3e-15
    item=dict(coordinates=coordinates,alpha=float(alpha),candidateArithmeticRelativeDifference=update_error,
              completeOriginalNSPassed=bool(report['completeOriginalNSPassed']),originalRmsDecreasePassed=bool(report['originalRmsDecreasePassed']))
    if report['completeOriginalNSPassed']:
        q=c['Q'];info=c['info'];rate=c['rate'];rr=q[:,0];uu=q[:,1:4]/rr[:,None]
        tt=(q[:,4]-.5*rr*np.sum(uu*uu,axis=1))/(rr*cv);assert np.all(rr>0) and np.all((tt>=80)&(tt<=2000))
        assert relative(np.column_stack([rr,uu,tt]),info['primitive'])<3e-15
        summed=np.zeros_like(q);flux=f['area'][:,None]*info['faceFlux']
        np.add.at(summed,m,-flux);np.add.at(summed,p[inside],flux[inside]);np.add.at(summed,owner,-w['area'][:,None]*info['wallFlux'])
        divergence=relative(summed/volume[:,None],rate);assert divergence<2e-12
        boundary=np.sum(flux[~inside].astype(np.longdouble),axis=0)+np.sum((w['area'][:,None]*info['wallFlux']).astype(np.longdouble),axis=0)
        denominator=np.maximum(1,np.sum(abs(flux[~inside]),axis=0)+np.sum(w['area'][:,None]*abs(info['wallFlux']),axis=0))
        conservation=float(np.max(abs(summed.astype(np.longdouble).sum(axis=0)+boundary)/denominator));assert conservation<1e-9
        merit=float(np.sqrt(np.mean((rate/scale)**2)));maximum=float(np.max(abs(rate/scale)))
        assert report['originalRmsDecreasePassed']==(merit<=(1-s['options']['sufficientDecrease']*alpha)*base_merit)
        assert abs(merit-report['metrics']['rmsMerit'])<3e-15*merit and abs(maximum-report['metrics']['maximumResidual'])<3e-15*maximum
        assert np.all(info['wallTemperatureInput']==300) and np.all(info['wallMechanicalPowerIntoFluid']==0)
        item.update(rmsMerit=merit,maximumResidual=maximum,fiveQuantityDivergenceRelativeDifference=divergence,
                    independentSpatialConservationRelative=conservation)
    assert not report['candidatePromotedToSolution'];checks.append(item)
assert s['report']['cases'][0]['fullConservativeControlExactlyReproduced']
assert not s['report']['spatialEquationChanged'] and not s['report']['physicalBoundaryChanged']
assert all(sha(ROOT/p)==h for p,h in bindings.items())
out=dict(status='passed-independent-real-same-direction-trial-path-and-full-flux-audit',case=case,processID=os.getpid(),
         cases=checks,fullIndependentNSSolver=False,candidatePromotedToSolution=False,physicalSolution=False,newProductionSolvedAreaM2=0,
         inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(out,indent=2)+'\n').encode('utf8'));print(json.dumps(out))
