"""具名标准干空气候选：表积分、联合压力温度范围、明确相对能量/熵身份。"""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/thermo/standard_dry_air_cp_v1.json'
MODEL='rah-standard-dry-air-cp-table-of13-v1'
DOMAIN='rah-standard-dry-air-low-density-v1'
TREF=298.15;PREF=100000.;TRANGE=(100.,2000.)

def pressure_ceiling(t):
    """模块一：低温收紧压力的已登记工程包络，不是文献的普适临界公式。"""
    ratio=np.asarray(t,dtype=float)/100.
    return np.minimum(1e6,1e4*((ratio*ratio)*ratio))

def validate(p,t):
    p,t=np.broadcast_arrays(np.asarray(p,dtype=float),np.asarray(t,dtype=float))
    if not p.size or np.any(~np.isfinite(p)) or np.any(~np.isfinite(t)):
        raise ValueError('STANDARD_AIR_INVALID_INPUT')
    if np.any((t<TRANGE[0])|(t>TRANGE[1])):raise ValueError('STANDARD_AIR_TEMPERATURE_OUTSIDE_DOMAIN')
    if np.any((p<=0)|(p>pressure_ceiling(t))):raise ValueError('STANDARD_AIR_PRESSURE_TEMPERATURE_OUTSIDE_DOMAIN')
    return p,t

class StandardDryAir:
    def __init__(self):
        d=json.loads(SOURCE.read_text(encoding='utf-8'))
        if d['model']!=MODEL or d['domain']!=DOMAIN:raise ValueError('STANDARD_AIR_MODEL_IDENTITY')
        self.data=d;self.T=np.array(d['temperatureK']);self.Cp=np.array(d['cpJkgK']);self.R=d['gasConstantJkgK']
        if (not np.isfinite(self.R) or self.R!=8.31451/.02896546 or self.T.shape!=self.Cp.shape
                or self.T.ndim!=1 or len(self.T)<2 or not np.isfinite(self.T).all() or not np.isfinite(self.Cp).all()):
            raise ValueError('STANDARD_AIR_INVALID_TABLE')
        if self.T[0]!=100 or self.T[-1]!=2000 or not np.all(np.diff(self.T)>0) or not np.all(self.Cp>self.R):
            raise ValueError('STANDARD_AIR_INVALID_TABLE')
        self.slopes=np.diff(self.Cp)/np.diff(self.T)
        self.H=np.r_[0,np.cumsum((self.Cp[:-1]+self.Cp[1:])*np.diff(self.T)/2)]
        self.S=np.r_[0,np.cumsum((self.Cp[:-1]-self.slopes*self.T[:-1])*np.log(self.T[1:]/self.T[:-1])+self.slopes*np.diff(self.T))]
        i=np.flatnonzero(self.T==TREF)
        if len(i)!=1:raise ValueError('STANDARD_AIR_MISSING_REFERENCE_NODE')
        self.H-=self.H[i[0]];self.S-=self.S[i[0]]

    def evaluate(self,p,t):
        """模块二：积分同一Cp表保证dh/dT=Cp、de/dT=Cv；域外明确拒绝。"""
        p,t=validate(p,t);i=np.clip(np.searchsorted(self.T,t,side='right')-1,0,len(self.T)-2)
        dt=t-self.T[i];grad=self.slopes[i];cp=self.Cp[i]+grad*dt
        h=self.H[i]+self.Cp[i]*dt+grad*dt*dt/2
        s=self.S[i]+(self.Cp[i]-grad*self.T[i])*np.log(t/self.T[i])+grad*dt-self.R*(np.log(p)-np.log(PREF))
        cv=cp-self.R
        rho=p/(self.R*t)
        if np.any(rho<=0) or np.any(~np.isfinite(rho)) or np.any(~np.isfinite(s)):
            raise ValueError('STANDARD_AIR_NUMERICAL_RANGE')
        return dict(pressurePa=p,temperatureK=t,cpJkgK=cp,cvJkgK=cv,enthalpyJkg=h,sensibleEnthalpyJkg=h,
            internalEnergyJkg=h-self.R*t,entropyJkgK=s,densityKgm3=rho,soundSpeedMps=np.sqrt(cp/cv*self.R*t))

    def dictionary(self):
        """模块三：准备原生表字典，实际原生联合域守卫尚须独立集成验证。"""
        native_ru=6.02214e23*1.38065e-23
        values='\n'.join('(%.17g %.17g)'%x for x in zip(self.T,self.Cp))
        return 'specie { molWeight %.17g; }\n'%(1000*native_ru/self.R)+\
            'thermodynamics { hf 0; sf 0; Cp { values (\n'+values+'\n); } }\ntransport { mu 0; Pr 0.72; }\n'

    def contract(self):
        """模块四：独立模型身份和能量基准；不复用旧NIST或OTH状态资格。"""
        return dict(model=MODEL,domain=DOMAIN,temperatureRangeK=list(TRANGE),pressureUpperFormula='min(1000000,10000*(T/100)^3)',
            pressureLowerExclusivePa=0.,moleFractions=self.data['moleFractions'],gasConstantJkgK=self.R,
            sensibleEnthalpyReference='hs(298.15K)=0; hf=0; es=hs-R*T',entropyReference='relative s(298.15K,100000Pa)=0; not absolute chemical entropy',
            sourceSHA256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),oldNISTModelUnchanged=True,nativeDomainDispatchVerified=False,
            physicalQualification=False,engineeringHeatReady=False)
