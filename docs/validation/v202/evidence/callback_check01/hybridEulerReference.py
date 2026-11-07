# -*- coding: utf-8 -*-
"""低温同EOS二维弱斜激波及等熵膨胀参照；不提供实验或三维箭体资格。"""
import math
import numpy as np
from scipy.optimize import root, brentq
from scipy.integrate import quad
from eilmerAirThermo import evaluate, R


def control_reference(eos=None,gas_constant=None,model=None,upstream_pressure=5000.,upstream_temperature=300.,temperature_min=200.,entropy_pressure=101325.):
    """模块一：Rankine-Hugoniot守恒关系求10度弱支，不用定gamma激波表。"""
    # 原CEA默认保持；新的具名原生入口必须显式传入自身EOS和气体常数。
    evaluate=eos if eos is not None else globals()['evaluate']
    R=gas_constant if gas_constant is not None else globals()['R']
    theta=math.radians(10);p1=float(upstream_pressure);T1=float(upstream_temperature)
    gas1=evaluate(p1,T1);rho1=float(gas1['densityKgm3'])
    U1=3*float(gas1['soundSpeedMps']);h1=float(gas1['enthalpyJkg']);h0=h1+U1**2/2
    def equations(x):
        beta,r=x;un1=U1*math.sin(beta);un2=un1/r
        p2=p1+rho1*un1**2*(1-1/r);T2=p2/(rho1*r*R)
        h2=float(evaluate(p2,T2)['enthalpyJkg'])
        return [(h2+un2**2/2-h1-un1**2/2)/U1**2,
                beta-math.atan(math.tan(beta)/r)-theta]
    solved=root(equations,[math.radians(27),1.7],tol=1e-11)
    if not solved.success or max(abs(np.asarray(equations(solved.x))))>1e-11:
        raise ValueError('同EOS激波参照未闭合')
    beta,r=solved.x;un1=U1*math.sin(beta);ut=U1*math.cos(beta);un2=un1/r
    p2=p1+rho1*un1**2*(1-1/r);T2=p2/(rho1*r*R);g2=evaluate(p2,T2)
    U2=math.hypot(un2,ut);s2=float(g2['entropyJkgK'])
    if not (math.asin(1/3)<beta<math.radians(40) and r>1 and temperature_min<T2<800
            and U2>g2['soundSpeedMps'] and s2>gas1['entropyJkgK']):
        raise ValueError('不是目标低温弱压缩支')
    # 模块二：简单波dtheta=sqrt(M^2-1)*dU/U；该温域cp=dh/dT。
    def integrand(T):
        g=evaluate(p2,T);u2=2*(h0-g['enthalpyJkg']);m2=u2/g['soundSpeedMps']**2
        if not temperature_min<=T<800 or m2<=1: raise ValueError('膨胀参照超出低温超声速域')
        return float(g['cpJkgK']*np.sqrt(m2-1)/u2)
    angle=lambda T:quad(integrand,T,T2,epsabs=1e-13,epsrel=1e-12)[0]
    T3=brentq(lambda T:angle(T)-theta,temperature_min,T2,xtol=1e-11)
    s_atm=float(evaluate(entropy_pressure,T3)['entropyJkgK'])
    p3=entropy_pressure*math.exp((s_atm-s2)/R);g3=evaluate(p3,T3)
    U3=math.sqrt(2*(h0-g3['enthalpyJkg']))
    states=[]
    for name,p,T,U,turn in [('upstream',p1,T1,U1,0),('compressed',p2,T2,U2,theta),
                           ('expanded',p3,T3,U3,0)]:
        g=evaluate(p,T);d={k:float(v) for k,v in g.items()}
        d.update(name=name,velocityMps=[-U*math.cos(turn),0.,U*math.sin(turn)],
                 totalEnthalpyJkg=float(g['enthalpyJkg']+U**2/2),Mach=float(U/g['soundSpeedMps']))
        states.append(d)
    residuals=dict(massRelative=abs(g2['densityKgm3']*un2/(rho1*un1)-1),
        normalMomentumRelative=abs((p2+g2['densityKgm3']*un2**2)/(p1+rho1*un1**2)-1),
        totalEnthalpyRelative=max(abs(d['totalEnthalpyJkg']/h0-1) for d in states),
        entropyExpansionAbsolute=abs(float(g3['entropyJkgK'])-s2),
        expansionAngleRad=abs(angle(T3)-theta),
        cpIntegralEnthalpyAbsolute=abs(quad(lambda t:float(evaluate(p1,t)['cpJkgK']),T1,T2,
                                             epsabs=1e-8,epsrel=1e-12)[0]-(g2['enthalpyJkg']-h1)))
    result=dict(schema='RAH.HybridEulerReference.v1',states=states,
                shockAngleDeg=math.degrees(beta),densityRatio=float(r),
                residuals={k:float(v) for k,v in residuals.items()},
                temperatureRangeK=[temperature_min,800],physicalQualification=False,
                scope='same-EOS 2D inviscid reference, not experimental heat or rocket validation')
    if model is not None:result['gasModel']=model
    return result
