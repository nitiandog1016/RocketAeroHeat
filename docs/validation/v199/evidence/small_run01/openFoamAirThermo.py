"""工程固定组分 NASA7 与 OpenFOAM 原生 janaf 的显式相容接口。"""
from pathlib import Path
import hashlib
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MODEL = 'rah-nasa7-air-of13-sensible-e-v1'
R_UNIVERSAL = 8.31446261815324
Y = np.array([.755, .232, .013])
M = np.array([.0280134, .0319988, .039948])
LOW = np.array([[3.53100528,-1.23660987e-4,-5.02999433e-7,2.43530612e-9,-1.40881235e-12,-1046.97628,2.96747468],
 [3.78245636,-2.99673416e-3,9.84730201e-6,-9.68129509e-9,3.24372837e-12,-1063.94356,3.65767573],
 [2.5,0,0,0,0,0,0]])
HIGH = np.array([[2.95257626,1.39690040e-3,-4.92631603e-7,7.86010367e-11,-4.60755321e-15,-923.948688,5.87188762],
 [3.28253784,1.48308754e-3,-7.57966669e-7,2.09470555e-10,-2.16717794e-14,-1088.45772,5.45323129],
 [2.5,0,0,0,0,0,0]])
FACTORS = Y*R_UNIVERSAL/M
R = float(FACTORS.sum())
MIX_MOLAR_MASS_GMOL = float(1000/(Y/M).sum())
# 原生发行包采用舍入的NA/k。用等效EOS参数保持目标R；它不是重新定义物理分子量。
FOAM_R_UNIVERSAL = 6.02214e23*1.38065e-23
EOS_MOLAR_PARAMETER = FOAM_R_UNIVERSAL*1000/R
LOW_MIX = FACTORS@LOW/R
HIGH_MIX = FACTORS@HIGH/R
TEMPERATURE_RANGE = (200., 2500.)


def bind_engineering_source():
    """模块一：逐项核对当前MATLAB权威系数，防止两个实现静默分叉。"""
    path=ROOT/'airThermalProperties.m'; text=path.read_text(encoding='utf-8')
    for name,expected in [('mass',Y),('molecularMass',M),('low',LOW),('high',HIGH)]:
        literal=re.search(r'\b'+name+r'=\[([^]]+)\]',text).group(1).replace('...','')
        parsed=np.array([[float(x) for x in row.split()] for row in literal.split(';')])
        if not np.array_equal(parsed.reshape(expected.shape),expected):
            raise ValueError('工程物性源系数已变化: '+name)
    if 'R=8.31446261815324;' not in text:raise ValueError('工程气体常数改变')
    return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def absolute_properties(T, coefficients=None):
    """模块二：按物种累加Cp、绝对参考焓和标准熵，独立于原生等效系数计算。"""
    t=np.asarray(T,dtype=float)
    if np.any(~np.isfinite(t)) or np.any((t<200)|(t>2500)):
        raise ValueError('NASA7原生入口仅支持200—2500 K，不裁温或外推')
    a=np.where((t[...,None,None]<1000),LOW,HIGH) if coefficients is None else coefficients
    z=t[...,None]
    cp=np.sum(FACTORS*(a[...,0]+a[...,1]*z+a[...,2]*z**2+a[...,3]*z**3+a[...,4]*z**4),axis=-1)
    ha=np.sum(FACTORS*(a[...,0]*z+a[...,1]*z**2/2+a[...,2]*z**3/3+a[...,3]*z**4/4+a[...,4]*z**5/5+a[...,5]),axis=-1)
    s=np.sum(FACTORS*(a[...,0]*np.log(z)+a[...,1]*z+a[...,2]*z**2/2+a[...,3]*z**3/3+a[...,4]*z**4/4+a[...,6]),axis=-1)
    return cp,ha,s


HF=float(absolute_properties(298.15)[1])


def evaluate(p,T):
    """模块三：与原生sensibleInternalEnergy对齐；负es允许，不误作热力学负温。"""
    p,t=np.broadcast_arrays(np.asarray(p,dtype=float),np.asarray(T,dtype=float))
    if np.any(~np.isfinite(p)) or np.any(p<=0):raise ValueError('压力必须正且有限')
    cp,ha,s=absolute_properties(t);cv=cp-R;rho=p/(R*t)
    return dict(temperatureK=t,pressurePa=p,densityKgm3=rho,cpJkgK=cp,cvJkgK=cv,
        enthalpyJkg=ha,sensibleEnthalpyJkg=ha-HF,internalEnergyJkg=ha-HF-R*t,
        entropyJkgK=s-R*np.log(p/101325.),soundSpeedMps=np.sqrt(cp/cv*R*t))


def mixture_dictionary():
    """模块四：不改原生内核的等效单组分配置；质量组分仍在独立合同中保存。"""
    coeff=lambda x:'('+' '.join(format(v,'.17g') for v in x)+')'
    return ('specie { molWeight %.17g; }\n'%EOS_MOLAR_PARAMETER+
        'thermodynamics { Tlow 200; Thigh 2500; Tcommon 1000;\n'+
        ' lowCpCoeffs '+coeff(LOW_MIX)+';\n highCpCoeffs '+coeff(HIGH_MIX)+'; }\n'+
        'transport { mu 0; Pr 0.72; }\n')


def contract():
    """模块五：明示能量零点/温域及分段跳跃，任何结果不能自动获得壁热资格。"""
    jump=float(absolute_properties(1000,HIGH)[1]-absolute_properties(1000,LOW)[1])
    return dict(model=MODEL,source=bind_engineering_source(),temperatureRangeK=list(TEMPERATURE_RANGE),
        gasConstantJkgK=R,molarMassGmol=MIX_MOLAR_MASS_GMOL,massFractions=Y.tolist(),
        nativeUniversalConstantJmolK=FOAM_R_UNIVERSAL,nativeEOSMolarParameter=EOS_MOLAR_PARAMETER,
        molarParameterMeaning='RR_native/R_target; effective EOS parameter, not physical molecular mass',
        hfJkg=HF,referenceTemperatureK=298.15,energy='es=ha-hf-R*T',
        engineTotalEnergy='E=es+0.5*|U|^2; engineering h0=E+hf+R*T; ha=es+hf+R*T',
        enthalpyJumpAt1000KJkg=jump,originalJumpRetained=True,
        chemistry='fixed composition, no dissociation',physicalQualification=False,
        nativeLimitBehavior='warn and clamp; any warning rejects run, not a pre-query guarantee')
