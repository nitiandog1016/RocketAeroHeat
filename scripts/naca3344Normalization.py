"""模块用途：原文恢复因子曲线和局部热参数对照；尚未闭合的测量段界显式拒绝。"""
import numpy as np

from naca3344Model import DATA,D,source_data,heat_parameter,recovery_temperature
from su2ThermalRuntime import read,sha


def published_recovery(s):
    """采用Figure3(b)原式12实线的图读Cr，不用散点平衡壁温或任意sqrt(Pr)代替。"""
    folder=DATA/'normalization';manifest=read(folder/'manifest.json')
    for name,h in manifest['files'].items():assert sha(folder/name)==h
    data=np.genfromtxt(str(folder/'recovery_factor.csv'),delimiter=',',names=True)
    s=np.asarray(s,float)/D
    assert np.isfinite(s).all() and np.all((s>=min(data['s_over_D']))&(s<=max(data['s_over_D'])))
    return np.interp(s,data['s_over_D'],data['Cr'])


def compare_local(s,q_into_wall,Tw,edge,edge_method):
    """模块二：可追溯的点值对照，不把缺失的原实验有限段平均伪装成已完成验证。"""
    assert edge_method and edge_method!='freestream-substituted-for-local-edge'
    contract,_,heat=source_data()
    coefficient=published_recovery(s)
    Tr=recovery_temperature(edge['temperatureK'],contract['freeStream']['totalTemperature_K'],coefficient,
        'NACA-TN3344-equation12-Figure3b-published-reduction-curve')
    result=heat_parameter(s,q_into_wall,Tw,edge,Tr,dict(edge=edge_method,
        recovery='NACA-TN3344-equation10/12-Figure3b-digitized'))
    sd=np.asarray(s)/D
    assert np.all((sd>=min(heat['s_over_D']))&(sd<=max(heat['s_over_D'])))
    result.update(recoveryTemperatureK=Tr,recoveryFactor=coefficient,
        measuredParameterPointInterpolation=np.interp(sd,heat['s_over_D'],heat['Nu_s_over_sqrt_Re_s']),
        definition='pointwise diagnostic; electrode-segment average still required',
        originalTransportTableEquivalent=False,measurementOperatorQualified=False,
        physicalQualification=False)
    return result


def require_experimental_qualification(segment_bounds,edge_qualified,converged,grid_and_domain_qualified):
    """模块三：正式实验对照必须有原测量段界及数值证据，不能用中点近似获得资格。"""
    assert segment_bounds is not None,'Figure1(c)各经线电极段界尚未逐点匹配，拒绝实验资格'
    bounds=np.asarray(segment_bounds,float)
    assert bounds.shape==(20,2) and np.isfinite(bounds).all() and np.all(bounds[:,1]>bounds[:,0])
    assert edge_qualified and converged and grid_and_domain_qualified
