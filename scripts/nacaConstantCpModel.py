"""模块用途：NACA低焓参考的唯一定比热定义；正状态计算、开发范围与物理资格分开。"""
from pathlib import Path
import numpy as np

from naca3344Model import DATA, source_data
from su2ThermalRuntime import sha

MODEL = 'SU2-IDEAL_GAS-gamma1.4-Sutherland-Pr0.72-NACA-v1'
GAMMA = 1.4
R = 8.31451/.02896546
PR = .72
MU_REF = 1.716e-5
MU_T_REF = 273.15
MU_S = 110.4
CP = GAMMA/(GAMMA-1)*R
CV = CP-R
DEVELOPMENT_TEMPERATURE_RANGE_K = (150., 360.)


def qualification_range():
    """仅返回预登记的比较/开发检查范围；此名称不表示已取得物理资格。"""
    contract, wall, _ = source_data()
    return dict(comparisonTemperatureRangeK=[float(contract['freeStream']['staticTemperature_K']),
                float(np.max(wall['Tw_K']))],
        developmentTemperatureRangeK=list(DEVELOPMENT_TEMPERATURE_RANGE_K),
        physicalQualified=False, experimentalQualified=False,
        meaning='predeclared comparison and development applicability checks; not a solved-field envelope')


def _inputs(first, second):
    a, b = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    if (a.ndim > 1 or b.ndim > 1 or a.size == 0 or b.size == 0
            or (a.size != 1 and b.size != 1 and a.shape != b.shape)):
        raise ValueError('状态输入须为标量或同尺寸一维数组')
    a, b = np.broadcast_arrays(a, b)
    a, b = a.reshape(-1), b.reshape(-1)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('状态输入必须全部有限')
    if np.any(a <= 0) or np.any(b <= 0):
        raise ValueError('PT或rhoe输入必须严格为正；不裁温或补能量')
    return a, b


def in_development_range(temperature):
    """只返回150—360K开发检查掩码；超界不改温度，也不说明正状态物理公式非法。"""
    values = np.asarray(temperature, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('开发范围检查不接受非有限温度')
    return ((values >= DEVELOPMENT_TEMPERATURE_RANGE_K[0])
            & (values <= DEVELOPMENT_TEMPERATURE_RANGE_K[1]))


def state(mode, first, second):
    """SI状态：PT或rhoe；rhoe中的e是单位质量静内能，不是rhoE或含动能总能量。"""
    if mode not in ('PT', 'rhoe'):
        raise ValueError('本模型仅接受PT或rhoe；不得混入OTH焓/内能定义')
    a, b = _inputs(first, second)
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            if mode == 'PT':
                p, t = a, b
                rho = p/(R*t)
            else:
                rho, energy = a, b
                t = energy/CV
                p = rho*R*t
            mu = MU_REF*(t/MU_T_REF)**1.5*(MU_T_REF+MU_S)/(t+MU_S)
            values = dict(T=t.copy(), rho=rho, p=p, e=CV*t, h=CP*t,
                cp=np.full(t.shape, CP), cv=np.full(t.shape, CV),
                gamma=np.full(t.shape, GAMMA), sound=np.sqrt(GAMMA*R*t),
                mu=mu, conductivity=mu*CP/PR)
    except FloatingPointError as error:
        raise ValueError('定比热状态算术溢出或非法') from error
    if not all(np.isfinite(value).all() and np.all(value > 0) for value in values.values()):
        raise ValueError('定比热输出不是有限正状态')
    return values


def freestream(mach, reynolds, total_temperature, diameter):
    """由M/Re_D/T0/D重建相容静态条件；全程使用本模型能量参考。"""
    anchors = np.asarray([mach, reynolds, total_temperature, diameter], dtype=float)
    if not np.isfinite(anchors).all() or np.any(anchors <= 0):
        raise ValueError('M、Re_D、T0、D须全部有限且严格为正')
    t = float(total_temperature/(1+.5*(GAMMA-1)*mach*mach))
    trial = state('PT', 101325., t)
    speed = float(mach*trial['sound'][0])
    density = float(reynolds*trial['mu'][0]/(speed*diameter))
    pressure = density*R*t
    actual = state('PT', pressure, t)
    return dict(Mach=float(mach), ReynoldsBasedOnDiameter=float(reynolds),
        totalTemperature_K=float(total_temperature), diameter_m=float(diameter),
        staticTemperature_K=t, pressure_Pa=float(pressure), density_kg_m3=density,
        speed_m_s=speed, soundSpeed_m_s=float(actual['sound'][0]),
        viscosity_Pas=float(actual['mu'][0]), conductivity_WmK=float(actual['conductivity'][0]),
        staticEnergy_Jkg=float(actual['e'][0]), staticEnthalpy_Jkg=float(actual['h'][0]),
        totalEnergy_Jkg=float(actual['e'][0]+.5*speed*speed), totalEnthalpy_Jkg=float(CP*total_temperature),
        derivedStateIsMeasured=False)


def conditions():
    """固定原实验锚点及Tw，不改旧合同、不读取任何CFD场。"""
    contract, wall, heat = source_data()
    c = contract['freeStream']
    assert c['Mach'] == 1.97 and c['ReynoldsBasedOnDiameter'] == 2.28e6
    assert c['totalTemperature_K'] == 570.*5/9 and contract['geometry']['diameter_m'] == .1016
    assert len(wall) == 22 and len(heat) == 20
    answer = freestream(c['Mach'], c['ReynoldsBasedOnDiameter'], c['totalTemperature_K'],
                        contract['geometry']['diameter_m'])
    paths = [Path(__file__).resolve(), DATA/'contract.json', DATA/'source_manifest.json',
             DATA/'wall_temperature.csv', DATA/'heat_parameter.csv']
    answer.update(model=MODEL, gamma=GAMMA, gasConstant_JkgK=R, cp_JkgK=CP, cv_JkgK=CV,
        Prandtl=PR, sutherland=dict(muRef_Pas=MU_REF, temperatureRef_K=MU_T_REF, S_K=MU_S),
        energyReference='e=cv*T and h=cp*T; zero at T=0 by caloric model, no OTH offset',
        qualificationRange=qualification_range(), physicalQualified=False,
        wallTemperatureRangeK=[float(np.min(wall['Tw_K'])), float(np.max(wall['Tw_K']))],
        sourceSHA256={str(path): sha(path) for path in paths})
    return answer
