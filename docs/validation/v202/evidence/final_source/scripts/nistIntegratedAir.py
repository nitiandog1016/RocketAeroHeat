"""新具名NIST固定组分积分比热空气；原NASA7与工程默认不变。"""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/thermo/nist_air_shomate_v1.json'
MODEL='rah-nist-integrated-air-of13-v1'
RU=8.31446261815324
TREF=298.15
PREF=100000.
TRANGE=(100.,2500.)
DATA=json.loads(SOURCE.read_text(encoding='utf-8'))
Y=np.array([DATA[k]['massFraction'] for k in ['N2','O2','Ar']])
M=np.array([DATA[k]['molarMassKgMol'] for k in ['N2','O2','Ar']])
R=float(np.sum(Y*RU/M));X=(Y/M)/np.sum(Y/M)
SMIX=float(-np.sum(Y*RU/M*np.log(X)))
NATIVE_RU=6.02214e23*1.38065e-23
EOS_MOLAR_PARAMETER=NATIVE_RU*1000/R


def validate_temperature(t):
    t=np.asarray(t,dtype=float)
    if np.any(~np.isfinite(t)) or np.any((t<TRANGE[0])|(t>TRANGE[1])):raise ValueError('新NIST模型只支持100—2500K，不裁温/外推')
    return t


def segment_properties(t,a):
    """模块一：原NIST Shomate公式，H原始kJ/mol在此显式转J/mol。"""
    z=np.asarray(t)/1000.;A,B,C,D,E,F,G,H=a
    cp=A+B*z+C*z*z+D*z**3+E/z**2
    h=1000*(A*z+B*z*z/2+C*z**3/3+D*z**4/4-E/z+F-H)
    s=A*np.log(z)+B*z+C*z*z/2+D*z**3/3-E/(2*z*z)+G
    return cp,h,s


def species_properties(name,t,continuous=True):
    """模块二：同一分段Cp积分得到连续H/S；原始拟合接缝仍可单独读取。"""
    t=validate_temperature(t)
    if name=='Ar':
        cp=np.full(t.shape,2.5*RU);return cp,cp*(t-TREF),DATA['Ar']['standardEntropy29815JmolK']+cp*np.log(t/TREF)
    segments=DATA[name]['segments'];cp=np.zeros(t.shape);rawh=np.zeros(t.shape);raws=np.zeros(t.shape)
    h=np.zeros(t.shape);s=np.zeros(t.shape)
    href=0.;sref=0.
    for j,row in enumerate(segments):
        lo,hi=row['rangeK'];a=row['coefficients'];at=np.clip(t,lo,hi)
        c,hh,ss=segment_properties(at,a);_,hlo,slo=segment_properties(lo,a)
        mask=(t>=lo)&((t<hi) if j<len(segments)-1 else (t<=hi))
        cp=np.where(mask,c,cp);rawh=np.where(mask,hh,rawh);raws=np.where(mask,ss,raws)
        h+=hh-hlo;s+=ss-slo
        _,hr,sr=segment_properties(np.clip(TREF,lo,hi),a);href+=hr-hlo;sref+=sr-slo
    if not continuous:return cp,rawh,raws
    anchor=segment_properties(TREF,segments[0]['coefficients'])[2]
    return cp,h-href,s-sref+anchor


def mixture_reference(t,continuous=True):
    cp,h,s=[np.zeros(np.asarray(t).shape) for _ in range(3)]
    for i,name in enumerate(['N2','O2','Ar']):
        a,b,c=species_properties(name,t,continuous);cp+=Y[i]/M[i]*a;h+=Y[i]/M[i]*b;s+=Y[i]/M[i]*c
    return cp,h,s+SMIX


class IntegratedAir:
    """模块三：固定0.5K表，独立实现分段线性Cp及其H/S积分供原生对照。"""
    def __init__(self):
        self.T=np.unique(np.r_[np.arange(100.,2500.25,.5),TREF])
        self.Cp=mixture_reference(self.T)[0]
        self.slopes=np.diff(self.Cp)/np.diff(self.T)
        stepH=(self.Cp[:-1]+self.Cp[1:])*np.diff(self.T)/2
        stepS=(self.Cp[:-1]-self.slopes*self.T[:-1])*np.log(self.T[1:]/self.T[:-1])+self.slopes*np.diff(self.T)
        self.H=np.r_[0.,np.cumsum(stepH)];self.S=np.r_[0.,np.cumsum(stepS)]
        i=int(np.where(self.T==TREF)[0][0]);self.H-=self.H[i];self.S-=self.S[i]
        self.sf=float(mixture_reference(TREF)[2])

    def evaluate(self,p,t):
        p,t=np.broadcast_arrays(np.asarray(p,dtype=float),validate_temperature(t))
        if np.any(~np.isfinite(p)) or np.any((p<=0)|(p>PREF)):raise ValueError('当前低压气体合同要求0<p<=100000Pa')
        i=np.clip(np.searchsorted(self.T,t,side='right')-1,0,len(self.T)-2)
        dt=t-self.T[i];grad=self.slopes[i];cp=self.Cp[i]+grad*dt
        h=self.H[i]+self.Cp[i]*dt+grad*dt*dt/2
        s=self.S[i]+(self.Cp[i]-grad*self.T[i])*np.log(t/self.T[i])+grad*dt+self.sf-R*np.log(p/PREF)
        cv=cp-R;rho=p/(R*t);e=h-R*t
        return dict(temperatureK=t,pressurePa=p,cpJkgK=cp,cvJkgK=cv,enthalpyJkg=h,sensibleEnthalpyJkg=h,
            internalEnergyJkg=e,entropyJkgK=s,densityKgm3=rho,soundSpeedMps=np.sqrt(cp/cv*R*t))

    def dictionary(self):
        """模块四：只生成原生字典；等效molWeight保持实际R，不冒称实际分子量改变。"""
        values='\n'.join('(%.17g %.17g)'%(t,c) for t,c in zip(self.T,self.Cp))
        return ('specie { molWeight %.17g; }\n'%EOS_MOLAR_PARAMETER+
            'thermodynamics { hf 0; sf %.17g; Cp { values (\n'%self.sf+values+'\n); } }\ntransport { mu 0; Pr 0.72; }\n')


AIR=IntegratedAir()
evaluate=AIR.evaluate


def contract():
    """模块五：分别登记数据来源、相容表示、物理适用域与未取得的热资格。"""
    return dict(model=MODEL,source=SOURCE.relative_to(ROOT).as_posix(),sourceSHA256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        temperatureRangeK=list(TRANGE),pressureRangePa=[0,PREF],pressureLowerExclusive=True,
        massFractions=Y.tolist(),moleFractions=X.tolist(),molarMassKgMol=M.tolist(),gasConstantJkgK=R,
        nativeEOSMolarParameter=EOS_MOLAR_PARAMETER,nativeUniversalConstantJmolK=NATIVE_RU,
        physicalMixtureMolarMassGmol=float(1000/np.sum(Y/M)),hfJkg=0.,referenceTemperatureK=TREF,
        standardPressurePa=PREF,standardEntropyJkgK=AIR.sf,idealMixingEntropyJkgK=SMIX,tableNodes=len(AIR.T),
        energy='hs=integral(Cp,298.15,T); ha=hs (hf=0); es=hs-R*T; h0=ha+|U|^2/2',
        integration='piecewise-linear Cp integrated exactly; H and S continuous; source Shomate seams reported separately',
        nativeClass='hIcoTabulatedThermo<perfectGas<specie>>',nativeOutOfRange='NonUniformTable Fatal; no clamp or extrapolation',
        chemistry='fixed N2/O2/Ar ideal gas; no dissociation or condensation',transport='Euler only; viscosity and heat conduction not qualified',
        physicalQualification=False,engineeringHeatReady=False)
