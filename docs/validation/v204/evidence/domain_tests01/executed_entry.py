"""独立审核联合域拒绝、能量/熵身份、源模型保持及新旧边界区别。"""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from standardDryAir import StandardDryAir,pressure_ceiling,MODEL
from nistIntegratedAir import AIR as OLD
from scipy.optimize import brentq

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    start=time.monotonic();out=ROOT/'outputs/v204_validation/domain_tests01';out.mkdir(exist_ok=False)
    model=StandardDryAir();rejected=[]
    samples=[('T-low',1000,99.99),('T-high',1000,2000.01),('T-nan',1000,np.nan),('T-inf',1000,np.inf),
        ('p-zero',0,300),('p-negative',-1,300),('p-nan',np.nan,300),('p-inf',np.inf,300),
        ('cold-compressed',100000,100),('cap',1000000.01,1000),('curve',float(pressure_ceiling(250))+.01,250),
        ('empty',[],[]),('mixed-batch',[1000,100000],[300,100])]
    for name,p,T in samples:
        try:model.evaluate(p,T)
        except ValueError as exc:rejected.append(dict(case=name,message=str(exc)))
        else:raise AssertionError(name)
    for T in [100.,132.53,200.,250.,464.1588833612779,1000.,2000.]:
        p=float(pressure_ceiling(T));model.evaluate(p,T)
        try:model.evaluate(np.nextafter(p,np.inf),T)
        except ValueError:rejected.append(dict(case='next-float-above-%.8g'%T,message='strict boundary rejected'))
        else:raise AssertionError('strict pressure boundary')
    a=model.evaluate(100000,298.15);assert float(a['enthalpyJkg'])==0 and float(a['entropyJkgK'])==0
    assert float(a['internalEnergyJkg'])==-model.R*298.15
    t=np.array([150.,226.50908361133003,400.,1274.,1900.]);dt=.001
    lo=model.evaluate(1000,t-dt);hi=model.evaluate(1000,t+dt);mid=model.evaluate(1000,t)
    dH=float(np.max(abs((hi['enthalpyJkg']-lo['enthalpyJkg'])/(2*dt)/mid['cpJkgK']-1)))
    dE=float(np.max(abs((hi['internalEnergyJkg']-lo['internalEnergyJkg'])/(2*dt)/mid['cvJkgK']-1)))
    dS=float(np.max(abs((hi['entropyJkgK']-lo['entropyJkgK'])/(2*dt)/(mid['cpJkgK']/t)-1)))
    assert max(dH,dE,dS)<2e-8
    # 改变参考压力不会扩域；旧生产合同仍拒绝原v203压力。
    try:OLD.evaluate(100292.29488236419,1231.3523721818547)
    except ValueError:oldClosed=True
    else:raise AssertionError('旧模型被意外放宽')
    for p,h in json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8')).items():assert sha(ROOT/p)==h
    # 同一Mach/高度的新候选总态重新推导，不能继承旧气体的速度/总焓。
    pinf=1197.0279007383137;Tinf=226.50908361133003;g=model.evaluate(pinf,Tinf);speed=5*float(g['soundSpeedMps'])
    h0=float(g['enthalpyJkg'])+speed*speed/2
    T0=brentq(lambda temp:float(model.evaluate(1000,temp)['enthalpyJkg'])-h0,100,2000,xtol=1e-10)
    s0=float(model.evaluate(1000,T0)['entropyJkgK']);si=float(model.evaluate(1000,Tinf)['entropyJkgK']);p0=pinf*np.exp((s0-si)/model.R)
    model.evaluate(p0,T0)
    r=dict(passed=True,rejections=rejected,derivativeRelativeErrors=dict(h=dH,e=dE,s=dS),oldPressureContractUnchanged=oldClosed,
        protectedUnchanged=29,model=MODEL,designedFreestreamSpeedMps=speed,designedVelocityMps=[-speed*np.cos(np.deg2rad(20)),0.,speed*np.sin(np.deg2rad(20))],
        isentropicDesign=dict(T0=T0,p0=p0,h0=h0,actualObserved=False),nativeTested=False,CFDRuns=0,builds=0,
        elapsedSeconds=time.monotonic()-start,sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),ROOT/'scripts/standardDryAir.py',ROOT/'data/thermo/standard_dry_air_cp_v1.json']})
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');print(json.dumps(r))

if __name__=='__main__':main()
