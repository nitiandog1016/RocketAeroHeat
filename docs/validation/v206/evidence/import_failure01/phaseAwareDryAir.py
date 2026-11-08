"""相态有界标准干空气：完整来源温区、相容能量/熵与具名低密度输运。"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/thermo/phase_aware_dry_air_v1.json'
MODEL='rah-phase-aware-dry-air-of13-v1'
DOMAIN='rah-lemmon-gas-low-density-60-2000-v1'
TRANGE=(60.,2000.);TREF=298.15;PREF=100000.

def dew_pressure(t):
    """模块一：Lemmon2000表10露点压力，只在该式发表的温区使用。"""
    t=np.asarray(t,dtype=float)
    if np.any(~np.isfinite(t)) or np.any((t<59.75)|(t>132.6312)):raise ValueError('DEW_T_OUTSIDE_SOURCE')
    theta=1-t/132.6312
    value=-.1567266*theta**.5-5.539635*theta+.7567212*theta**2.5-3.514322*theta**4
    return 3785020.*np.exp(132.6312/t*value)

def pressure_ceiling(t):
    t=np.asarray(t,dtype=float)
    if np.any(~np.isfinite(t)) or np.any((t<60)|(t>2000)):raise ValueError('PHASE_AIR_T_OUTSIDE_DOMAIN')
    ratio=t/100.;upper=np.minimum(1e6,1e4*((ratio*ratio)*ratio))
    cold=t<=132.6312
    return np.where(cold,np.minimum(upper,.1*dew_pressure(np.where(cold,t,132.6312))),upper)

def validate(p,t):
    p,t=np.broadcast_arrays(np.asarray(p,dtype=float),np.asarray(t,dtype=float))
    if not p.size or np.any(~np.isfinite(p)) or np.any(~np.isfinite(t)):raise ValueError('PHASE_AIR_INVALID_INPUT')
    ceiling=pressure_ceiling(t)
    if np.any((p<=0)|(p>ceiling)):raise ValueError('PHASE_AIR_PT_OUTSIDE_DOMAIN')
    return p,t

class PhaseAwareDryAir:
    def __init__(self):
        d=json.loads(SOURCE.read_text(encoding='utf-8'));self.data=d
        if d['model']!=MODEL or d['domain']!=DOMAIN:raise ValueError('PHASE_AIR_IDENTITY')
        self.T=np.array(d['temperatureK']);self.Cp=np.array(d['cpJkgK']);self.R=d['gasConstantJkgK']
        if self.R!=8.31451/.02896546 or self.T.shape!=self.Cp.shape or self.T.ndim!=1:raise ValueError('PHASE_AIR_TABLE')
        if not np.isfinite(self.T).all() or not np.isfinite(self.Cp).all() or self.T[0]!=60 or self.T[-1]!=2000:
            raise ValueError('PHASE_AIR_TABLE')
        if not np.all(np.diff(self.T)>0) or not np.all(self.Cp>self.R):raise ValueError('PHASE_AIR_TABLE')
        dt=np.diff(self.T);self.slopes=np.diff(self.Cp)/dt
        self.H=np.r_[0,np.cumsum((self.Cp[:-1]+self.Cp[1:])*dt/2)]
        self.S=np.r_[0,np.cumsum((self.Cp[:-1]-self.slopes*self.T[:-1])*np.log(self.T[1:]/self.T[:-1])+self.slopes*dt)]
        ids=np.flatnonzero(self.T==TREF)
        if len(ids)!=1:raise ValueError('PHASE_AIR_REFERENCE')
        self.H-=self.H[ids[0]];self.S-=self.S[ids[0]]

    def thermal(self,t):
        """模块二：单独的量热数学函数，不能据此授予p/T物理资格。"""
        t=np.asarray(t,dtype=float)
        if not t.size or np.any(~np.isfinite(t)) or np.any((t<60)|(t>2000)):raise ValueError('PHASE_AIR_T_OUTSIDE_DOMAIN')
        i=np.clip(np.searchsorted(self.T,t,side='right')-1,0,len(self.T)-2);dt=t-self.T[i];grad=self.slopes[i]
        cp=self.Cp[i]+grad*dt;h=self.H[i]+self.Cp[i]*dt+grad*dt*dt/2
        s=self.S[i]+(self.Cp[i]-grad*self.T[i])*np.log(t/self.T[i])+grad*dt
        return cp,h,s

    def evaluate(self,p,t):
        p,t=validate(p,t);cp,h,s=self.thermal(t);rho=p/(self.R*t)
        s=s-self.R*(np.log(p)-np.log(PREF))
        if np.any(rho<=0) or np.any(~np.isfinite(rho)) or np.any(~np.isfinite(s)):raise ValueError('PHASE_AIR_NUMERICAL_RANGE')
        return dict(pressurePa=p,temperatureK=t,cpJkgK=cp,cvJkgK=cp-self.R,enthalpyJkg=h,sensibleEnthalpyJkg=h,
            internalEnergyJkg=h-self.R*t,entropyJkgK=s,densityKgm3=rho,soundSpeedMps=np.sqrt(cp/(cp-self.R)*self.R*t))

    def transport(self,p,t):
        """模块三：Lemmon2004式2/5；不将低压气体自动判作连续介质。"""
        state=self.evaluate(p,t);t=state['temperatureK'];p=state['pressurePa'];x=np.log(t/103.3)
        omega=np.exp(.431-.4623*x+.08406*x*x+.005341*x**3-.00331*x**4)
        mu=2.66958e-8*np.sqrt(28.9586*t)/(.36**2*omega)
        tau=132.6312/t;k=.001308*(mu*1e6)+.001405*tau**(-1.1)-.001036*tau**(-.3)
        if np.any(~np.isfinite(mu)) or np.any(~np.isfinite(k)) or np.any(mu<=0) or np.any(k<=0):raise ValueError('PHASE_AIR_TRANSPORT')
        return dict(dynamicViscosityPas=mu,conductivityWmK=k,Prandtl=state['cpJkgK']*mu/k,
            meanFreePathViscosityDefinitionM=mu/p*np.sqrt(np.pi*self.R*t/2))

    def inverse(self,p,value,kind='h'):
        """模块四：严格有界焓/内能反解，最终p/T再核验，不裁目标或状态。"""
        if kind not in ('h','e'):raise ValueError('PHASE_AIR_INVERSE_KIND')
        p,value=np.broadcast_arrays(np.asarray(p,dtype=float),np.asarray(value,dtype=float))
        if not p.size or not np.isfinite(p).all() or not np.isfinite(value).all():raise ValueError('PHASE_AIR_INVALID_INPUT')
        def caloric(t):return self.thermal(t)[1]-(self.R*t if kind=='e' else 0.)
        lo,hi=map(float,[caloric(60.),caloric(2000.)])
        if np.any((value<lo)|(value>hi)):raise ValueError('PHASE_AIR_ENERGY_OUTSIDE_TABLE')
        ts=np.array([60. if v==lo else 2000. if v==hi else brentq(lambda t:float(caloric(t))-v,60.,2000.,xtol=1e-10)
            for v in value.flat]).reshape(value.shape)
        self.evaluate(p,ts)
        return ts

    def reference_state(self,p,te,speed,tw,recovery):
        """模块五：相同气体的恢复焓/参考焓，r由有来源的区域方法显式提供。"""
        p,te,speed,tw,recovery=np.broadcast_arrays(*[np.asarray(v,dtype=float) for v in (p,te,speed,tw,recovery)])
        if np.any(~np.isfinite(speed)) or np.any(speed<0) or np.any(~np.isfinite(recovery)) or np.any((recovery<0)|(recovery>1)):
            raise ValueError('PHASE_AIR_RECOVERY_INPUT')
        edge=self.evaluate(p,te);wall=self.evaluate(p,tw)
        haw=edge['enthalpyJkg']+recovery*speed**2/2
        taw=self.inverse(p,haw);href=.5*(edge['enthalpyJkg']+wall['enthalpyJkg'])+.22*(haw-edge['enthalpyJkg'])
        tref=self.inverse(p,href);ref=self.evaluate(p,tref)
        return dict(recoveryTemperatureK=taw,referenceTemperatureK=tref,recoveryEnthalpyJkg=haw,referenceEnthalpyJkg=href,
            referenceTemperatureDerivativeByWallTemperature=.5*wall['cpJkgK']/ref['cpJkgK'],
            thermalDrivingEnthalpyJkg=haw-wall['enthalpyJkg'],**self.transport(p,tref))

    def dictionary(self):
        native_ru=6.02214e23*1.38065e-23
        values='\n'.join('(%.17g %.17g)'%x for x in zip(self.T,self.Cp))
        return 'specie { molWeight %.17g; }\n'%(1000*native_ru/self.R)+'thermodynamics { hf 0; sf 0; Cp { values (\n'+values+'\n); } }\ntransport { mu 0; Pr 0.72; }\n'

    def contract(self):
        return dict(model=MODEL,domain=DOMAIN,temperatureRangeK=list(TRANGE),gasConstantJkgK=self.R,moleFractions=[.7812,.2096,.0092],
            pressureDomain='0<p<=min(1e6,1e4*(T/100)^3,0.1*p_dew(T) when T<=132.6312)',
            energyReference='hf=0; hs(298.15K)=0; es=hs-RT',entropyReference='relative s(298.15K,100000Pa)=0',
            sourceSHA256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            transportModel='Lemmon2004-dilute-mu-k-v1',transportRequiresContinuumCheck=True,
            nativeFinalStateGuardVerified=False,physicalQualification=False,engineeringHeatReady=False)
