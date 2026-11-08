"""新低温模型/相态/恢复与参考焓接口，独立量热和差分检验。"""
from pathlib import Path
import sys,json,time,shutil
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from phaseAwareDryAir import PhaseAwareDryAir,pressure_ceiling,dew_pressure,SOURCE
from standardDryAir import StandardDryAir
from idealAirProperties import IdealAir
from prepareOpenFoamControl import save,sha

def main():
    start=time.monotonic();out=ROOT/'outputs/v206_validation/python_tests01';out.mkdir(exist_ok=False)
    for name in ['scripts/phaseAwareDryAir.py','tests/testPhaseAwareAir.py']:
        shutil.copyfile(ROOT/name,out/Path(name).name)
    g=PhaseAwareDryAir();old=StandardDryAir();helper=IdealAir();failed=[]
    def reject(name,func):
        try:func()
        except ValueError as e:failed.append(dict(case=name,error=str(e)));return
        raise AssertionError('未拒绝 '+name)
    # 模块一：完整低温端点/中间点与成熟理想项独立比对，不重跑旧物性网格。
    t=np.array([60.,60.01,61.75,67.3,75.5,88.2,99.994,100.,226.50908361133003,298.15,750.,1500.,2000.])
    p=pressure_ceiling(t)*.5;s=g.evaluate(p,t);tr=g.transport(p,t)
    e_inv=g.inverse(p,s['internalEnergyJkg'],'e');h_inv=g.inverse(p,s['enthalpyJkg'],'h')
    inverse_error=float(max(abs(e_inv-t).max(),abs(h_inv-t).max()));assert inverse_error<2e-8
    legacy=helper.state('PT',p,t);target=legacy['h']-g.data['oldOTHEnthalpyAt29815Jkg']
    independent_error=float(np.max(abs(s['enthalpyJkg']-target)/(s['cpJkgK']*t)));assert independent_error<2e-6
    for label,pp,tt in [('low',1,59.99),('high',1,2000.01),('p-zero',0,300),('nan',1,np.nan),('inf',np.inf,300),
        ('cold-phase',float(dew_pressure(60))*.5,60),('cap',1000001,1000),('empty',[],[])]:
        reject(label,lambda pp=pp,tt=tt:g.evaluate(pp,tt))
    for tt in [60,62,70,100,132.6312,1000,2000]:
        reject('next-p-'+str(tt),lambda tt=tt:g.evaluate(np.nextafter(float(pressure_ceiling(tt)),np.inf),tt))
    reject('outer-broadcast',lambda:g.evaluate(np.ones((2,1)),np.full((1,2),300)))
    reject('e-below',lambda:g.inverse(1,float(g.thermal(60)[1])-g.R*60-1,'e'))
    reject('e-above',lambda:g.inverse(1,float(g.thermal(2000)[1])-g.R*2000+1,'e'))
    reject('inverse-final-phase',lambda:g.inverse(100000,float(g.thermal(60)[1]),'h'))
    reject('recovery-invalid',lambda:g.reference_state(1000,300,500,300,1.01))
    reject('recovery-energy-outside',lambda:g.reference_state(1000,300,5000,300,.9))
    reject('old-model-stays-closed',lambda:old.evaluate(1,80))
    # 模块二：恢复/参考焓由独立Cp积分审核，壁温导数用中心差分验证。
    cases=np.array([[10,70,1200,300,.85],[1000,180,1500,300,.9],[50000,700,900,600,.88],[1000,1200,0,1800,0.]])
    p2,te,u,tw,r=cases.T;ref=g.reference_state(p2,te,u,tw,r);dh=[]
    def integral(a,b):
        sign=1 if b>=a else -1;lo,hi=sorted((a,b));points=g.T[(g.T>lo)&(g.T<hi)]
        return sign*quad(lambda x:float(g.thermal(x)[0]),lo,hi,points=points,limit=len(points)+50,epsabs=1e-7,epsrel=1e-11)[0]
    for i in range(len(cases)):
        recovered=integral(te[i],float(ref['recoveryTemperatureK'][i]));href=integral(te[i],float(ref['referenceTemperatureK'][i]))
        dh.extend([abs(recovered-r[i]*u[i]**2/2),abs(href-(.5*integral(te[i],tw[i])+.22*r[i]*u[i]**2/2))])
    step=.01;plus=g.reference_state(p2,te,u,tw+step,r);minus=g.reference_state(p2,te,u,tw-step,r)
    derivative=float(np.max(abs((plus['referenceTemperatureK']-minus['referenceTemperatureK'])/(2*step)-ref['referenceTemperatureDerivativeByWallTemperature'])))
    assert max(dh)<1e-6 and derivative<1e-7
    assert ref['thermalDrivingEnthalpyJkg'][-1]<0,'热驱动必须保留负值'
    # 跨语言样本写一次，MATLAB读取同模型状态及独立差分/积分结果。
    fixture=dict(pressurePa=p.tolist(),temperatureK=t.tolist(),states={k:np.asarray(v).tolist() for k,v in dict(s,**tr).items()},
        referenceInputs=cases.tolist(),referenceOutputs={k:np.asarray(v).tolist() for k,v in ref.items()},dataSHA256=sha(SOURCE))
    save(out/'matlab_fixture.json',fixture)
    report=dict(passed=True,rejections=failed,positiveStates=len(t),inverseMaxErrorK=inverse_error,independentIdealMaxScaled=independent_error,
        referenceIntegralMaxAbsJkg=max(dh),wallDerivativeMaxAbs=derivative,oldModelUnchanged=True,physicalQualification=False,
        nativeRuns=0,sourceSHA256={n:sha(ROOT/n) for n in ['scripts/phaseAwareDryAir.py','tests/testPhaseAwareAir.py','data/thermo/phase_aware_dry_air_v1.json']},
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','positiveStates','inverseMaxErrorK','referenceIntegralMaxAbsJkg','wallDerivativeMaxAbs','elapsedSeconds']}))

if __name__=='__main__':main()
