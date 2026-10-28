# -*- coding: utf-8 -*-
"""比较同一真实NS切向方向的两种候选坐标可行半径，不重算NS或更新已接受状态。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1];assert case.replace('_','').isalnum();target=VAL/(case+'_primitive_trial_radius.json');assert not target.exists()
source=VAL/case/'wall_velocity_krylov_evidence.mat';domain_file=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


bindings={p.relative_to(ROOT).as_posix():sha(p) for p in [source,domain_file]}
s=loadmat(source,simplify_cells=True,variable_names=['result','gas']);r=s['result'];gas=s['gas']
tr=r['trace'];assert isinstance(tr,dict)
q=tr['initialQ'];dq=tr['linearUpdateU']*r['residualScale'];rho=q[:,0];u=q[:,1:4]/rho[:,None];cv=gas['R']/(gas['gamma']-1)
T=(q[:,4]-.5*rho*np.sum(u*u,axis=1))/(rho*cv)
dr=dq[:,0];du=(dq[:,1:4]-u*dr[:,None])/rho[:,None]
dt=(dq[:,4]-np.sum(u*dq[:,1:4],axis=1)+(.5*np.sum(u*u,axis=1)-cv*T)*dr)/(rho*cv)
d=loadmat(domain_file,simplify_cells=True,variable_names=['domain'])['domain'];nf=int(d['quality']['sourceFaces'])


def location(i):
    return dict(cellID=int(i+1),layer=int(i//nf+1),sourceFacet=int(i%nf+1),centroidM=d['cells']['centroid'][i].tolist())


constraints=[]
for name,delta,margin,condition in [('positiveDensity',dr,rho,dr<0),
                                    ('minimumTemperature',dt,T-gas['temperatureRange'][0],dt<0),
                                    ('maximumTemperature',-dt,gas['temperatureRange'][1]-T,dt>0)]:
    alpha=np.full(rho.shape,np.inf);alpha[condition]=-margin[condition]/delta[condition]
    i=int(np.argmin(alpha));constraints.append(dict(constraint=name,firstExitAlpha=float(alpha[i]),location=location(i)))
radius=min(a['firstExitAlpha'] for a in constraints)
trials=[]
for alpha in [1/512,1/1024,1/2048,1/4096]:
    rr=rho+alpha*dr;uu=u+alpha*du;tt=T+alpha*dt
    candidate=np.column_stack([rr,rr[:,None]*uu,rr*(cv*tt+.5*np.sum(uu*uu,axis=1))])
    legal=bool(np.all(rr>0)&np.all(tt>=gas['temperatureRange'][0])&np.all(tt<=gas['temperatureRange'][1]))
    # 恢复五量后二次独立EOS回代，不以构造时的T冒充真实候选温度。
    direct=(candidate[:,4]-.5*np.sum(candidate[:,1:4]**2,axis=1)/candidate[:,0])/(candidate[:,0]*cv)
    if legal:np.testing.assert_allclose(direct,tt,rtol=2e-13,atol=2e-10)
    trials.append(dict(alpha=alpha,cellStateLegal=legal,minDensity=float(rr.min()),minTemperature=float(tt.min()),
                       maxTemperature=float(tt.max()),fullNSOrInterfaceRecomputed=False))
out=dict(status='verified-same-direction-primitive-tangent-cell-state-radius',case=case,processID=os.getpid(),
         constraints=constraints,firstCellStateExitAlpha=radius,candidates=trials,flowStateChanged=False,
         newActualNSEvaluations=0,newProductionSolvedAreaM2=0,physicalSolution=False,inputSHA256=bindings,scriptSHA256=sha(__file__))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8');print(json.dumps({k:out[k] for k in ['status','constraints','firstCellStateExitAlpha','candidates']}))
