"""模块职责：标准Air与SU2具名内部能量坐标之间的显式交换。

这里只转换量的表示，不求流场、不裁状态；完整气体/原生分派另行接入。
流动数组顺序固定为[rho, rho*u, rho*v, rho*w, rho*E]；通量和残差
使用对应的五分量，SA的非守恒nu_tilde必须单独传递。
"""
import numpy as np
from phaseAwareDryAir import MODEL, DOMAIN
from phaseAwareViscousAir import TRANSPORT

PHYSICAL_REFERENCE='rah-hs298.15-zero-v1'
SOLVER_REFERENCE='rah-su2-phase-air-eplus300k-v1'
SHIFT_JKG=300000.0
SCHEMA='rah-flow-energy-packet-v1'
QUANTITIES=('state','flux','residual','increment')


def contract():
    """模块一：对外声明气体、输运和两个零点，不允许默认为历史OTH。"""
    return dict(schema=SCHEMA,gasModel=MODEL,domain=DOMAIN,transport=TRANSPORT,
                physicalReference=PHYSICAL_REFERENCE,solverReference=SOLVER_REFERENCE,
                solverMinusPhysicalJkg=SHIFT_JKG,components=['mass','momentumX','momentumY','momentumZ','totalEnergy'],
                units='SI',nativeIntegrationVerified=False,physicalQualification=False)


def _offset(reference):
    if reference==PHYSICAL_REFERENCE:return 0.0
    if reference==SOLVER_REFERENCE:return SHIFT_JKG
    raise ValueError('ENERGY_REFERENCE_UNKNOWN')


def specific_energy(values,source_reference,target_reference):
    """模块二：转换比内能或比焓；不用于Cp、温度、熵或单位体积能量。"""
    delta=_offset(target_reference)-_offset(source_reference)
    data=np.asarray(values,dtype=float)
    if not data.size or not np.isfinite(data).all():raise ValueError('ENERGY_VALUES_INVALID')
    with np.errstate(over='ignore',invalid='ignore'):out=data+delta
    if not np.isfinite(out).all():raise ValueError('ENERGY_NUMERICAL_RANGE')
    return out


def flow_packet(values,reference,quantity):
    """模块三：创建显式身份的数据包；调用者负责数据确实采用所声明的零点。"""
    packet=dict(schema=SCHEMA,gasModel=MODEL,domain=DOMAIN,transport=TRANSPORT,
                reference=reference,quantity=quantity,units='SI',values=np.array(values,dtype=float,copy=True))
    _validate_packet(packet)
    return packet


def _validate_packet(packet):
    required={'schema','gasModel','domain','transport','reference','quantity','units','values'}
    if not isinstance(packet,dict) or set(packet)!=required:raise ValueError('ENERGY_PACKET_SCHEMA')
    for key,value in [('schema',SCHEMA),('gasModel',MODEL),('domain',DOMAIN),('transport',TRANSPORT),('units','SI')]:
        if packet[key]!=value:raise ValueError('ENERGY_PACKET_IDENTITY_'+key)
    _offset(packet['reference'])
    if packet['quantity'] not in QUANTITIES:raise ValueError('ENERGY_PACKET_QUANTITY')
    values=np.asarray(packet['values'],dtype=float)
    if values.ndim<1 or values.shape[-1]!=5 or not values.size:raise ValueError('ENERGY_PACKET_SHAPE')
    if not np.isfinite(values).all():raise ValueError('ENERGY_VALUES_INVALID')
    # 密度正性仅约束状态；残差、增量及通量的质量分量允许零或负。
    if packet['quantity']=='state' and np.any(values[...,0]<=0):raise ValueError('ENERGY_STATE_DENSITY')
    return values


def convert_flow(packet,target_reference):
    """模块四：同一线性坐标变换用于守恒量、通量和逐单元有符号残差。

    必须在残差求范数之前调用；本函数不接受RMS/log残差，也不能从这些量恢复
    物理能量残差。返回独立数组且更新零点身份，相同目标不会重复加C。
    """
    data=_validate_packet(packet)
    delta=_offset(target_reference)-_offset(packet['reference'])
    out=data.copy()
    with np.errstate(over='ignore',invalid='ignore'):out[...,4]+=delta*data[...,0]
    if not np.isfinite(out).all():raise ValueError('ENERGY_NUMERICAL_RANGE')
    return flow_packet(out,target_reference,packet['quantity'])
