"""当前合法温域内的膨胀能力参照；不反算表外状态、不定位真实失败单元。"""
import math,time,shutil,json
import numpy as np
from scipy.integrate import quad
from prepareOpenFoamControl import ROOT,sha,save
from standardDryAir import StandardDryAir,MODEL

def main():
    start=time.monotonic();out=ROOT/'outputs/v205_validation/cold_scope01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');gas=StandardDryAir()
    admission=json.loads((ROOT/'outputs/v205_validation/rocket_run01/admission.json').read_text(encoding='utf-8'))
    condition=admission['condition'];p1=condition['pressurePa'];t1=condition['temperatureK'];t2=100.
    first=gas.evaluate(p1,t1);speed=condition['Mach']*float(first['soundSpeedMps']);h0=float(first['enthalpyJkg'])+speed**2/2
    # 模块一：仅在合法区间积分已有简单波关系，节点处分段避免把表插值误当光滑多项式。
    def kernel(t):
        state=gas.evaluate(p1,t);u2=2*(h0-state['enthalpyJkg']);m2=u2/state['soundSpeedMps']**2
        if m2<=1:raise ValueError('超声速简单波范围失效')
        return float(state['cpJkgK']*math.sqrt(m2-1)/u2)
    points=gas.T[(gas.T>t2)&(gas.T<t1)]
    angle,error=quad(kernel,t2,t1,points=points,limit=len(points)+50,epsabs=1e-12,epsrel=1e-11)
    endpoint_at_p1=gas.evaluate(p1,t2)
    p2=p1*math.exp((float(endpoint_at_p1['entropyJkgK'])-float(first['entropyJkgK']))/gas.R)
    last=gas.evaluate(p2,t2);u2=math.sqrt(2*(h0-float(last['enthalpyJkg'])))
    integral=quad(lambda t:float(gas.evaluate(p1,t)['cpJkgK']),t2,t1,points=points,limit=len(points)+50,epsabs=1e-8,epsrel=1e-12)[0]
    dh_error=abs(integral-(float(first['enthalpyJkg'])-float(last['enthalpyJkg'])))
    s_error=abs(float(last['entropyJkgK'])-float(first['entropyJkgK']))
    if dh_error>1e-7 or s_error>1e-9:raise ValueError('合法域积分/等熵闭合未过')
    result=dict(passed=True,model=MODEL,upstream=condition,evaluatedTemperatureRangeK=[t2,t1],
        sameEOSExpansionTurnDeg=math.degrees(angle),limitPressurePa=p2,limitMach=u2/float(last['soundSpeedMps']),
        cpIntegralEnthalpyErrorJkg=dh_error,entropyClosureErrorJkgK=s_error,quadratureErrorEstimateRad=error,
        relation='dtheta=sqrt(M^2-1)*dU/U; cp=dhs/dT; hs+U^2/2=constant',
        sourceURL='https://www.grc.nasa.gov/WWW/BGH/pranmyer.html',
        reusedVariableCpDerivation='scripts/hybridEulerReference.py; scripts/checkRocketExpansionTemperatureLimit.py',
        meaning='2D isentropic scope reference only; not actual 3D cell, transient temperature or unique failure cause',
        changedThermoRange=False,CFDRuns=0,builds=0,physicalQualification=False,
        sourceSHA256={p:sha(ROOT/p) for p in ['scripts/checkStandardAirColdScope.py','scripts/standardDryAir.py',
            'scripts/hybridEulerReference.py','scripts/checkRocketExpansionTemperatureLimit.py','docs/validation/v205/COLD_SCOPE_PLAN.md']},
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','sameEOSExpansionTurnDeg','limitPressurePa','limitMach','cpIntegralEnthalpyErrorJkg','elapsedSeconds']}))

if __name__=='__main__':main()
