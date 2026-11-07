"""只在原NASA7有效温域内检查膨胀到200K所需转角；不推断真实异常单元。"""
from pathlib import Path
import json,time,math,shutil
import numpy as np
from scipy.integrate import quad
from prepareOpenFoamControl import ROOT,sha,save
from openFoamAirThermo import evaluate,R


def main():
    start=time.monotonic();out=ROOT/'outputs/v201_validation/temperature_scope01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    p1=1197.0279007383137;T1=226.50908361133003;T2=200.;gas=evaluate(p1,T1)
    u1=5*float(gas['soundSpeedMps']);h1=float(gas['enthalpyJkg']);h0=h1+u1*u1/2
    # 模块一：复用已有同EOS简单波关系dtheta=sqrt(M²−1)dU/U，cp=dh/dT。
    def kernel(t):
        g=evaluate(p1,t);u2=2*(h0-g['enthalpyJkg']);m2=u2/g['soundSpeedMps']**2
        return float(g['cpJkgK']*np.sqrt(m2-1)/u2)
    angle,error=quad(kernel,T2,T1,epsabs=1e-13,epsrel=1e-12)
    p2=101325*math.exp((float(evaluate(101325,T2)['entropyJkgK'])-float(gas['entropyJkgK']))/R)
    g2=evaluate(p2,T2);u2=math.sqrt(2*(h0-float(g2['enthalpyJkg'])))
    integral=quad(lambda t:float(evaluate(p1,t)['cpJkgK']),T2,T1,epsabs=1e-9)[0]
    # 模块二：定gamma闭式和同gamma积分互验公式，保持其对照模型身份。
    gamma=1.4;cp=gamma*R/(gamma-1);speed=5*math.sqrt(gamma*R*T1);target=math.sqrt((speed**2+2*cp*(T1-T2))/(gamma*R*T2))
    def nu(m):return math.sqrt((gamma+1)/(gamma-1))*math.atan(math.sqrt((gamma-1)/(gamma+1)*(m*m-1)))-math.atan(math.sqrt(m*m-1))
    exact=nu(target)-nu(5)
    integ=quad(lambda t:cp*math.sqrt((speed**2+2*cp*(T1-t))/(gamma*R*t)-1)/(speed**2+2*cp*(T1-t)),T2,T1,epsabs=1e-13)[0]
    dh_error=abs(integral-(h1-float(g2['enthalpyJkg'])))
    assert dh_error<1e-8 and abs(exact-integ)<1e-13
    r=dict(passed=True,upstream=dict(Mach=5,pressurePa=p1,temperatureK=T1),limitTemperatureK=T2,
        sameEOSExpansionTurnDeg=math.degrees(angle),limitPressurePa=p2,limitMach=u2/float(g2['soundSpeedMps']),
        constantGammaExpansionTurnDeg=math.degrees(exact),constantGammaIntegralErrorRad=abs(exact-integ),
        cpIntegralEnthalpyErrorJkg=dh_error,quadratureErrorEstimateRad=error,
        sourceURL='https://www.grc.nasa.gov/WWW/BGH/pranmyer.html',
        existingImplementationSHA256=sha(ROOT/'scripts/hybridEulerReference.py'),
        evaluatedTemperatureRangeK=[T2,T1],changedThermoRange=False,
        meaning='2D isentropic feasibility check only; not 3D transient diagnosis or actual minimum temperature',
        physicalQualification=False,CFDRuns=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',r);print(json.dumps(r))


if __name__=='__main__':main()
