"""模块用途：半球柱几何、实测壁温和局部热参数；测量量与计算量严格分开。"""
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.integrate import quad

from idealAirProperties import IdealAir
from su2ThermalRuntime import ROOT, read, sha

DATA = ROOT/'tests/data/naca_tn3344'
D, R = .1016, .0508
SHOULDER = np.pi*R/2
END = SHOULDER+.127-R


@lru_cache(maxsize=1)
def source_data():
    """模块一：复用图读数据及逐字清单，拒绝静默改用别的实验或温度。"""
    contract = read(DATA/'contract.json')
    manifest = read(DATA/'source_manifest.json')
    assert sha(DATA/'source_manifest.json') == contract['sourceManifestSHA256']
    for row in manifest['measurements']:
        assert sha(DATA/row['path']) == row['sha256']
    wall = np.genfromtxt(str(DATA/'wall_temperature.csv'), delimiter=',', names=True)
    heat = np.genfromtxt(str(DATA/'heat_parameter.csv'), delimiter=',', names=True)
    return contract, wall, heat


def meridian(s, offset=0.):
    """模块二：保留真实 r=0 鼻尖；offset 沿外法向，仅用于外层几何。"""
    s, offset = np.broadcast_arrays(np.asarray(s, float), np.asarray(offset, float))
    assert np.isfinite(s).all() and np.isfinite(offset).all()
    assert np.all((s >= 0) & (s <= END+1e-14)) and np.all(offset >= 0)
    theta = np.minimum(s/R, np.pi/2)
    x = (R+offset)*np.cos(theta)-R-np.maximum(s-SHOULDER, 0)
    radius = (R+offset)*np.sin(theta)
    return np.stack((x, radius), axis=-1)


def wall_temperature(s):
    """模块三：测站间线性插值，末测站后固定末值；显式返回延拓位置。"""
    _, wall, _ = source_data()
    s = np.asarray(s, float)
    assert np.isfinite(s).all() and np.all((s >= 0) & (s <= END+1e-14))
    last = wall['s_over_D'][-1]*D
    value = np.interp(s, wall['s_over_D']*D, wall['Tw_K'])
    return value, s > last


def strip_mean(function, lower, upper):
    """模块四：按真实回转壁面积积分，可用于温度和段热功率，非简单节点均值。"""
    assert 0 <= lower < upper <= END+1e-14
    cuts = np.r_[source_data()[1]['s_over_D']*D, SHOULDER]
    points = np.unique(cuts[(cuts>lower)&(cuts<upper)]).tolist()
    weight = lambda s: float(meridian(s)[1])*2*np.pi
    area = quad(weight, lower, upper, points=points, epsabs=1e-13)[0]
    value = quad(lambda s: float(function(s))*weight(s), lower, upper,
                 points=points, epsabs=1e-10, limit=150)[0]/area
    return float(value), float(area)


def heat_parameter(s, q_into_wall, Tw, edge, recovery, provenance):
    """模块五：原式1/10/11的局部 Nu/√Re，必须显式提供外缘和恢复温度来源。"""
    required = ('temperatureK','densityKgM3','speedMS','viscosityPas','conductivityWmK')
    assert set(edge) == set(required) and provenance.get('edge') and provenance.get('recovery')
    arrays = np.broadcast_arrays(*[np.asarray(v, float) for v in
        [s, q_into_wall, Tw, recovery]+[edge[k] for k in required]])
    assert all(np.isfinite(a).all() for a in arrays)
    s,q,Tw,Tr,Te,rho,u,mu,k = arrays
    assert np.all(s > 0), '驻点 s=Ue=0 需专门极限，不能代入0/0'
    assert np.all((rho>0)&(u>0)&(mu>0)&(k>0)&(Te>0))
    assert np.all(Tw>Tr), '本加热壁基准必须 Tw>Tr'
    assert np.all(q<0), '本实验正入壁热流应为负，不能取绝对值掩盖符号'
    Re = rho*u*s/mu
    h = -q/(Tw-Tr)
    Nu = h*s/k
    return dict(parameter=Nu/np.sqrt(Re), Reynolds=Re, Nusselt=Nu, heatTransferCoefficient=h,
                qOut=-q, experimentalQualification=False)


def recovery_temperature(Te, T0, coefficient, source):
    """原式10；Cr须是有来源的局部量，不能把近鼻部未加热表面T当Tr。"""
    assert source and source != 'unheated-wall-near-nose'
    Te,T0,coefficient = np.broadcast_arrays(Te,T0,coefficient)
    assert np.isfinite(Te+T0+coefficient).all()
    assert np.all((Te>0)&(T0>=Te)&(coefficient>0)&(coefficient<=1))
    return Te+coefficient*(T0-Te)


def transport(temperature, pressure):
    """模块六：求解器同一物性只用于计算边缘量，不称作原实验物性表。"""
    thermo = IdealAir().state('PT', pressure, temperature)
    t = thermo['T']
    mu = 1.716e-5*(t/273.15)**1.5*(273.15+110.4)/(t+110.4)
    return thermo, mu, mu*thermo['cp']/.72
