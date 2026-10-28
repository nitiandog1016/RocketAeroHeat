"""保留真实面状态，通过公开压力曲线定位采样转换；不冻结真实求解器波系。"""
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq
from riemannReservoirReference import namespace

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v84_validation'
p=loadmat(OUT/'diagnosis01/finite_step_patch.mat',simplify_cells=True)['patch']
idx=int(np.flatnonzero(p['fluxGlobalFaceID']==238752)[0]);left=int(p['fluxMinus'][idx])-1;right=int(p['fluxPlus'][idx])-1
n=p['fluxNormal'][idx];g=p['gas']['gamma'];R=p['gas']['R'];cv=R/(g-1)
def evaluate(fraction):
    q=p['Q0']+fraction*(p['Q1']-p['Q0']);states=[]
    for row in [q[left],q[right]]:
        rho=row[0];v=row[1:4]/rho;pressure=(g-1)*(row[4]-.5*rho*np.dot(v,v));states.append((rho,float(v@n),pressure))
    a=namespace['exact'](*states,g);rho,u,pressure=states[1]
    assert a['pressureStar']>pressure
    speed=u+np.sqrt(g*pressure/rho)*np.sqrt((g+1)/(2*g)*a['pressureStar']/pressure+(g-1)/(2*g))
    return {'fractionOfQ22Q24Step':float(fraction),'rightShockSpeedMs':float(speed),'sampleRegion':a['sampleRegion'],'pressureStarPa':a['pressureStar']}
cross=brentq(lambda x:evaluate(x)['rightShockSpeedMs'],.125,.25,xtol=1e-14)
d=loadmat(OUT/'diagnosis01/finite_step_decomposition.mat',simplify_cells=True)
c=loadmat(OUT/'conditioning01/conditioning.mat',simplify_cells=True)
report={'status':'diagnosed-wave-crossing-smooth-curvature-and-state-metric-conditioning',
    'faceID':238752,'crossingFractionOfAcceptedStep':cross,'crossingAlpha':cross/512,
    'waveSamples':[evaluate(f) for f in [0,.0625,.125,cross,.25,.5,1]],
    'conservativeScaledJacobianCondition':float(np.linalg.cond(c['J'])),
    'primitiveMetricJacobianCondition':float(np.linalg.cond(c['A'])),
    'nonInteriorConvectionRemainderToTotal':float(np.linalg.norm(d['total']-d['conv'])/np.linalg.norm(d['total'])),
    'unchangedRegionAndChangedRegionContributionsAreNotOrthogonal':True,
    'physicalSolution':False,'newProductionSolvedAreaM2':0}
(OUT/'mechanism.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8');print(json.dumps(report))
