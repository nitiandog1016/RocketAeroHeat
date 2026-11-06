"""模块用途：用独立解析剖面核验外缘选择、坐标身份、拒绝和采样敏感性。"""
import argparse
from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from nacaBoundaryLayerEdge import extract_edge, FIELDS
from su2ThermalRuntime import sha
from validateIdealAir import save


PROVENANCE = dict(stateKind='synthetic-test',flowStateSource='analytic-polynomial-v170',
                  normalColumnSource='exact-straight-normal-coordinates',
                  thermodynamicModel='synthetic-calorically-perfect-air-cp1004.5-R287')
NORMAL = np.array([0.,1.,0.])
TANGENT = np.array([-1.,0.,0.])


def profile(count=401, shock=False):
    """已知有限厚度三次速度/温度剖面；真实场资格保持为假。"""
    y = np.r_[0.,np.geomspace(1e-6,.03,count-1)]
    eta = np.minimum(y/.002,1.)
    blend = eta*eta*(3-2*eta)
    u = 400*blend
    t = 330-80*blend
    p = np.full(count,60000.)
    if shock:
        outside = y >= .012
        u[outside], t[outside], p[outside] = 527.,178.,26000.
    cp = np.full(count,1004.5)
    mu = 1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)
    return dict(distanceM=y, coordinatesM=y[:,None]*NORMAL,
                velocityMS=u[:,None]*TANGENT, pressurePa=p,temperatureK=t,
                densityKgM3=p/(287*t),enthalpyJkg=cp*t,specificHeatJkgK=cp,
                viscosityPas=mu,conductivityWmK=mu*cp/.72)


def call(data, **kwargs):
    args=dict(normal=NORMAL,tangent=TANGENT,station_s_m=.05,provenance=PROVENANCE)
    args.update(kwargs)
    return extract_edge(data,**args)


def main(out):
    out.mkdir(parents=True,exist_ok=False)
    files=[Path(__file__).resolve(),ROOT/'scripts/nacaBoundaryLayerEdge.py']
    locks={str(p):sha(p) for p in files}
    for p in files:(out/p.name).write_bytes(p.read_bytes())
    # 独立用解析多项式求根，不调用生产选择器建立真值。
    exact=.002*brentq(lambda eta:3*eta*eta-2*eta**3-.99,0,1,xtol=1e-15)
    accepted=[]
    for name,data in [('finite-polynomial',profile()),('postshock-local-state',profile(shock=True))]:
        r=call(data);save(out/(name+'.json'),r)
        assert r['accepted'],r['status']
        error=abs(r['delta99M']/exact-1)
        assert error<.002
        assert r['edgeSample']['pressurePa']==60000.
        assert r['edge']['temperatureK']==250. and r['edge']['speedMS']==400.
        assert r['edgeSample']['distanceM']>r['delta99M']
        assert r['samplingSensitivity']['accepted']
        assert not r['experimentalQualification'] and not r['physicalQualification']
        k=r['edgeSample']['index']
        assert all(r['edgeSample'][f]==data[f][k] for f in FIELDS)
        assert abs(r['edge']['densityKgM3']*287*r['edge']['temperatureK']/60000-1)<1e-15
        if name=='postshock-local-state':
            assert r['strongOuterChangeFirstIndex'] is not None
            assert k<r['usableLastIndex'] and data['pressurePa'][-1]==26000.
        accepted.append(dict(name=name,delta99RelativeError=error,exactDelta99M=exact))
    data=profile();data['enthalpyJkg']+=1e9
    shifted=call(data);base=call(profile())
    assert shifted['accepted'] and shifted['edgeSample']['index']==base['edgeSample']['index']
    assert shifted['delta99M']==base['delta99M']
    accepted.append(dict(name='enthalpy-reference-shift',sameSelectedNode=True))

    cases=[]
    d=profile();d['velocityMS'][:]=0;cases.append(('stagnant',d,{},'stagnant-or-near-stagnation-column'))
    cases.append(('exact-tip',profile(),dict(station_s_m=0),'stagnation-limit-not-defined'))
    d=profile();d['velocityMS'][40:100]=[10,0,0]
    cases.append(('reverse-flow',d,{},'reverse-flow-outside-method-scope'))
    d=profile();d['velocityMS'][:,2]=.1*np.linalg.norm(d['velocityMS'],axis=1)
    cases.append(('crossflow',d,{},'crossflow-outside-method-scope'))
    d=profile();d['velocityMS'][0]=[-1,0,0]
    cases.append(('slipping-wall',d,{},'wall-no-slip-not-satisfied'))
    d=profile();d['coordinatesM'][70,2]=1e-5
    cases.append(('not-normal',d,{},'not-a-straight-normal-column'))
    d=profile();d['distanceM'][10]=d['distanceM'][9]
    cases.append(('duplicate-distance',d,{},'missing-wall-or-unordered-normal-distance'))
    d=profile();d['temperatureK'][0]=np.nan
    cases.append(('nonfinite',d,{},'nonfinite-or-malformed-profile'))
    cases.append(('no-provenance',profile(),dict(provenance={}), 'missing-local-solution-provenance'))
    provenance=dict(PROVENANCE,flowStateSource='freestream-substituted-for-local-edge')
    cases.append(('freestream-fallback',profile(),dict(provenance=provenance),'freestream-substitution-forbidden'))
    d=profile();d['velocityMS']=-400*(d['distanceM']/.03)[:,None]*np.array([1.,0.,0.])
    cases.append(('no-platform',d,{},'no-resolved-local-platform'))
    d=profile();d['velocityMS'][(d['distanceM']>=.0018)&(d['distanceM']<=.003)]*=1.3
    cases.append(('overshoot',d,{},None))
    d=profile();y=d['distanceM'];u=-d['velocityMS'][:,0]
    # 两个真实不同的平台，中间以连续坡度连接；禁止挑方便的平台。
    u=np.where(y<=.003,u*.5,np.where(y<.006,200+200*(y-.003)/.003,400))
    d['velocityMS']=u[:,None]*TANGENT
    d['temperatureK'][:]=250;d['enthalpyJkg'][:]=1004.5*250
    d['densityKgM3'][:]=60000/(287*250)
    d['viscosityPas'][:]=d['viscosityPas'][-1];d['conductivityWmK'][:]=d['conductivityWmK'][-1]
    cases.append(('two-platforms',d,{},'multiple-platforms-ambiguous'))
    d=profile(shock=True);y=d['distanceM'];mask=y>=.0008
    d['pressurePa'][mask]=26000;d['densityKgM3'][mask]=26000/(287*d['temperatureK'][mask])
    cases.append(('shock-inside-layer',d,{},None))
    d=profile(31)
    cases.append(('coarse-normal-sampling',d,{},None))
    rejected=[]
    for name,data,kwargs,status in cases:
        r=call(data,**kwargs);save(out/(name+'.json'),r)
        assert not r['accepted'] and r['edge'] is None,(name,r['status'])
        if status is not None:assert r['status']==status,(name,r['status'],status)
        rejected.append(dict(name=name,status=r['status']))
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    report=dict(passed=True,acceptedTests=accepted,rejectedTests=rejected,
                inputSHA256=locks,inputsStable=stable,realFlowColumnsTested=0,
                syntheticTestsAreCFDValidation=False,experimentalQualification=False,
                newQualifiedArea_m2=0)
    save(out/'report.json',report)
    print(dict(passed=True,accepted=len(accepted),rejected=len(rejected)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True)
    args=p.parse_args();main(Path(args.output).resolve())
