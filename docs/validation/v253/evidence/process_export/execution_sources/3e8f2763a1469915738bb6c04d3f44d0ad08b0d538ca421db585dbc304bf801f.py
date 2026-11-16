"""模块职责：一般流管尺度的CR-112180质量表与实际剖面入口匹配；不授物理资格。"""
import numpy as np
from scipy.integrate import solve_ivp,simpson
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from phaseAwareViscousAir import PhaseAwareViscousAir


class ShockMassModel:
    """模块一：同一相态气体的压缩根；保持法向质量/动量和完整总焓。"""
    def __init__(self,p,t,speed):
        self.air=PhaseAwareViscousAir();self.p=float(p);self.t=float(t);self.speed=float(speed)
        a=self.air.evaluate(p,t);self.rho=float(a['densityKgm3']);self.a=float(a['soundSpeedMps'])
        self.h=float(a['enthalpyJkg']);self.H=self.h+.5*speed**2;self.Q=self.rho*speed**2
        if speed<=self.a:raise ValueError('INLET_SUBSONIC_UPSTREAM')
        self.normal=self.shock(1.)
        self.t0=float(self.air.inverse(p,self.H))
        _,_,s0=self.air.thermal(self.t0);_,_,s2=self.air.thermal(self.normal['temperatureK'])
        self.p0=float(self.normal['pressurePa']*np.exp((s0-s2)/self.air.R))
        self.air.evaluate(self.p0,self.t0)

    def shock(self,sin_angle):
        if not np.isfinite(sin_angle) or not self.a/self.speed<=sin_angle<=1:raise ValueError('INLET_SHOCK_ANGLE_DOMAIN')
        v=self.speed*sin_angle;R=self.air.R
        if v==self.a:return dict(pressurePa=self.p,temperatureK=self.t,normalSpeed=v,entropyJkgK=float(self.air.evaluate(self.p,self.t)['entropyJkgK']))
        A=v*v/R;B=A-self.t;upper=(B+np.sqrt(B*B+4*A*(self.t-60)))/(2*A)
        cuts=[0.,upper];disc=B*B-4*A*(2000-self.t)
        if disc>0:cuts.extend([(B-np.sqrt(disc))/(2*A),(B+np.sqrt(disc))/(2*A)])
        cuts=sorted(set(x for x in cuts if 0<=x<=upper))
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            if self.t+B*(lo+hi)/2-A*((lo+hi)/2)**2>2000:continue
            def residual(z):
                temp=self.t+B*z-A*z*z
                # 只对解析区间端点的舍入还原准确表端，不裁未知物理态。
                if z in (lo,hi):
                    if abs(temp-60)<1e-9:temp=60.
                    if abs(temp-2000)<1e-9:temp=2000.
                cp,h,_=self.air.thermal(temp)
                mean=float(cp) if temp==self.t else float((h-self.h)/(temp-self.t))
                if abs(temp-self.t)<1e-3:mean=float(self.air.thermal((temp+self.t)/2)[0])
                return mean*(B-A*z)+.5*v*v*(z-2)
            if residual(lo)>0 and residual(hi)<=0:
                z=brentq(residual,lo,hi,xtol=1e-14);temp=self.t+B*z-A*z*z;un=v*(1-z)
                p=self.p+self.rho*v*(v-un);state=self.air.evaluate(p,temp)
                return dict(pressurePa=float(p),temperatureK=float(temp),normalSpeed=float(un),entropyJkgK=float(state['entropyJkgK']))
        raise ValueError('INLET_NO_PHASE_SHOCK_ROOT')

    def pressure_table(self):
        wmax=1-(self.a/self.speed)**2
        w=np.unique(np.r_[np.linspace(0,wmax,257),1e-6,1e-5,1e-4])
        p=np.array([self.shock(float(np.sqrt(1-x)))['pressurePa'] if x<wmax else self.p for x in w])
        return w,p,PchipInterpolator(w,p,extrapolate=False)

    def recover(self,pressure,sin_angle):
        """按匹配到的激波熵与物理压力恢复；该位置不是原物面流线同熵延拓。"""
        jump=self.shock(sin_angle);entropy=jump['entropyJkgK']
        target=entropy+self.air.R*np.log(pressure/1e5)
        lo=float(self.air.thermal(60.)[2]);hi=float(self.air.thermal(2000.)[2])
        if not lo<=target<=hi:raise ValueError('INLET_ENTROPY_OUTSIDE_TABLE')
        t=brentq(lambda T:float(self.air.thermal(T)[2])-target,60.,2000.,xtol=1e-10)
        state=self.air.evaluate(pressure,t);kinetic=2*(self.H-float(state['enthalpyJkg']))
        if kinetic<=0:raise ValueError('INLET_NONPOSITIVE_EDGE_KINETIC')
        return dict(pressurePa=float(pressure),temperatureK=float(t),speedMps=float(np.sqrt(kinetic)),
            densityKgm3=float(state['densityKgm3']),entropyJkgK=entropy,totalEnthalpyJkg=self.H,
            shockAngleDeg=float(np.rad2deg(np.arcsin(sin_angle))),shock=jump)


def integrate_metric_shock(model,stations,metric,physical_pressure,inlet):
    """模块二：式65—67的一般h推进；有限上游入口不能被隐式置零。"""
    s=np.asarray(stations,float);h=np.asarray(metric,float);p=np.asarray(physical_pressure,float)
    if len(s)<2 or h.shape!=s.shape or p.shape!=s.shape or not np.isfinite(np.r_[s,h,p]).all() or np.any(np.diff(s)<=0) or np.any(h<=0):raise ValueError('INLET_METRIC_DATA')
    if not isinstance(inlet,dict) or set(('massCoordinate','sinShockAngle','sourceSHA256'))-inlet.keys():raise ValueError('INLET_FINITE_SOURCE_REQUIRED')
    if inlet['massCoordinate']<=0 or len(inlet['sourceSHA256'])!=64:raise ValueError('INLET_FINITE_SOURCE_REQUIRED')
    hm=PchipInterpolator(s,h);pm=PchipInterpolator(s,p);_,_,psh=model.pressure_table()
    def rhs(x,y):
        u,psi=y
        if not model.a/model.speed<u<=1 or psi<=0:raise ValueError('INLET_SHOCK_ODE_DOMAIN')
        pb=model.p+(pm(x)-model.p)*(model.normal['pressurePa']-model.p)/(model.p0-model.p)
        return [(pb-psh(1-u*u))*hm(x)/(model.speed*psi),model.rho*model.speed*u*hm(x)]
    sol=solve_ivp(rhs,[s[0],s[-1]],[inlet['sinShockAngle'],inlet['massCoordinate']],t_eval=s,rtol=1e-9,atol=[1e-11,1e-13],method='DOP853')
    if not sol.success or np.any(np.diff(sol.y[1])<=0):raise ValueError('INLET_SHOCK_MASS_INTEGRATION')
    return sol.y


def spherical_source_table(model,radius,max_angle):
    """模块三：仅生成已声明球帽的上游质量供体，不延拓到后体或假称真实激波面。"""
    w,p,psh=model.pressure_table();probe=np.array([1e-3,2.5e-4,6.25e-5])
    d=np.array([(model.normal['pressurePa']-model.shock(np.sqrt(1-x))['pressurePa'])/x for x in probe])
    d1=(4*d[1]-d[0])/3;D=(4*d[2]-d[1])/3
    if abs(D-d1)/abs(D)>1e-4:raise ValueError('INLET_STAGNATION_DERIVATIVE')
    A=(model.normal['pressurePa']-model.p)/(model.Q+2*D)
    q0=1e-6;init=[1-A*q0,q0/2-(A+1/6)*q0*q0/4]
    def rhs(q,y):
        theta=np.sqrt(q);u,psi=y
        if psi<=0 or not model.a/model.speed<u<=1:raise ValueError('INLET_CAP_ODE_DOMAIN')
        pb=model.p+(model.normal['pressurePa']-model.p)*np.cos(theta)**2
        return [(pb-psh(1-u*u))/model.Q*np.sin(theta)/(2*theta*psi),u*np.sin(theta)/(2*theta)]
    theta=np.r_[0,np.geomspace(.001,max_angle,321)];theta[-1]=max_angle
    sol=solve_ivp(rhs,[q0,max_angle**2],init,t_eval=theta[1:]**2,rtol=1e-9,atol=[1e-11,1e-13],method='DOP853')
    if not sol.success:raise ValueError('INLET_CAP_INTEGRATION')
    u=np.r_[1,sol.y[0]];mass=model.rho*model.speed*radius**2*np.r_[0,sol.y[1]]
    if np.any(np.diff(mass)<=0) or np.any(np.diff(u)>1e-12):raise ValueError('INLET_CAP_MONOTONICITY')
    return dict(thetaRad=theta.tolist(),sinShockAngle=u.tolist(),massCoordinate=mass.tolist(),radiusM=radius,
        A=float(A),pressureDerivative=float(D),derivativeRelativeIndicator=float(abs(D-d1)/abs(D)),
        shockPressureTableW=w.tolist(),shockPressureTablePa=p.tolist(),physicalQualification=False,
        sourceKind='CR112180-B1-rotated-spherical-cap-engineering-mass-source')


def match_inlet_mass(model,table,pressure,mass_per_width,metric,label_jacobian):
    """模块四：固定已积分剖面质量，回溯同标签激波入口；不修改剖面去强凑外缘。"""
    if not np.isfinite([pressure,mass_per_width,metric,label_jacobian]).all() or min(pressure,mass_per_width,metric,label_jacobian)<=0:raise ValueError('INLET_MASS_OR_LABEL_DOMAIN')
    mass=mass_per_width*metric;source_mass=mass/label_jacobian
    coordinate=np.array(table['massCoordinate']);u=np.array(table['sinShockAngle'])
    if source_mass>coordinate[-1]:raise ValueError('INLET_MASS_OUTSIDE_SOURCE')
    fun=PchipInterpolator(coordinate,u,extrapolate=False);sin_angle=float(fun(source_mass))
    state=model.recover(pressure,sin_angle)
    state.update(matchedMassPerLabel=mass,sourceMassCoordinate=source_mass,labelJacobian=label_jacobian,
        profileMassPerWidth=mass_per_width,projectedMetricM=metric,sourceSinShockAngle=sin_angle,
        acceptedEdge=False,profileReadyForCFD=False)
    return state
