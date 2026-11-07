# -*- coding: utf-8 -*-
"""固定Eilmer5伪空气CEA9回调；保留原生能量/熵零点，不与NASA7混名。"""
from pathlib import Path
import hashlib
import numpy as np

R = 8.31451/.02896
COEFFICIENTS = np.array([
    [1.009950160e4,-1.968275610e2,5.009155110,-5.761013730e-3,1.066859930e-5,
     -7.940297970e-9,2.185231910e-12,-1.767967310e2,-3.921504225],
    [2.415214430e5,-1.257874600e3,5.144558670,-2.138541790e-4,7.065227840e-8,
     -1.071483490e-11,6.577800150e-16,6.462263190e3,-8.147411905]])


def evaluate(p, T):
    """模块一：原生系数线性混合。表外拒绝；不继承原生端点延拓的资格。"""
    p,T = np.broadcast_arrays(np.asarray(p,dtype=float),np.asarray(T,dtype=float))
    if not (np.isfinite(p).all() and np.isfinite(T).all() and np.all(p>0)
            and np.all((T>=200)&(T<=6000))):
        raise ValueError('CEA空气接口只接受正压、有限200—6000K状态')
    w=np.clip((T-800)/400,0,1)
    a=((1-w)[...,None]*COEFFICIENTS[0]+w[...,None]*COEFFICIENTS[1]).T
    # 只接标量或一维逐面序列，避免多维广播隐式改变原面顺序。
    if T.ndim>1: raise ValueError('物性输入必须为标量或逐面向量')
    cp=R*(a[0]/T**2+a[1]/T+a[2]+a[3]*T+a[4]*T**2+a[5]*T**3+a[6]*T**4)
    h=R*(-a[0]/T+a[1]*np.log(T)+a[2]*T+a[3]*T**2/2+a[4]*T**3/3+
         a[5]*T**4/4+a[6]*T**5/5+a[7])
    s=R*(-a[0]/(2*T**2)-a[1]/T+a[2]*np.log(T)+a[3]*T+a[4]*T**2/2+
         a[5]*T**3/3+a[6]*T**4/4+a[8]-np.log(p/101325.))
    return dict(pressurePa=p,temperatureK=T,densityKgm3=p/(R*T),internalEnergyJkg=h-R*T,
                enthalpyJkg=h,entropyJkgK=s,cpJkgK=cp,soundSpeedMps=np.sqrt(cp/(cp-R)*R*T),
                viscosityPas=1.716e-5*(T/273.)**1.5*(273+111)/(T+111),
                conductivityWmK=.0241*(T/273.)**1.5*(273+194)/(T+194))


def bind_model(gas_file, expected_sha):
    """模块二：绑定已原生核对的完整生成气体，不能只靠model名称接受别的EOS。"""
    path=Path(gas_file); raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected_sha: raise ValueError('气体文件SHA变化')
    canonical='\n'.join(s for s in raw.decode().splitlines() if not s.startswith('-- Auto-generated'))
    if hashlib.sha256(canonical.encode()).hexdigest()!=GAS_CANONICAL_SHA:
        raise ValueError('不是固定CompositeGas/air/CEA9生成文件')
    return dict(name='eilmer-air-cea-v1',nativeClass='CompositeGas/ThermPerfGasMixture',
                sourceCommit='f53f4609a0331d48efee69a4e4f3c3598378cc03',gasSHA256=expected_sha,
                temperatureRangeK=[200.,6000.],gasConstantJkgK=R,
                enthalpyReference='CEA9 official air polynomial, no added shift',
                entropyReference='CEA9 s(T) minus R*ln(p/101325 Pa); not s(298.15K)=0',
                blending='linear coefficients 800--1200 K; cp is not exact dh/dT there',
                physicalQualification=False)


# 生成器的日期注释除外，其余字节合同固定；来源由原生查询另行核对。
GAS_CANONICAL_SHA = '6544e040a47671f6b3a9477f7d73b25d573d722edbec950cdecc281ad71bd66a'
