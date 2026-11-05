"""固定理想空气核的平均流总态转换，供公开工况准备复用。"""
import math

from scipy.integrate import quad
from scipy.optimize import brentq

from su2IdealAirInterface import R


def total_to_static(model, mach, total_temperature, total_pressure):
    """模块一。保持Mach和平均流总温/总压，用总焓与等熵关系求静态。

    总态不包含SST湍动能k，求解器守恒能量由调用者另加k。
    60—2000K为固定核数值温区，任何无相容根的输入立即拒绝。
    """
    if not all(math.isfinite(v) for v in (mach, total_temperature, total_pressure)) or mach < 0:
        raise ValueError('Mach and total state must be finite; Mach must be nonnegative')
    total = model.state('PT', total_pressure, total_temperature)
    h0 = float(total['h'][0])

    def energy_balance(temperature):
        value = model.state('PT', total_pressure, temperature)
        return float(value['h'][0] + .5*(mach*value['sound'][0])**2 - h0)

    temperature = total_temperature
    if mach > 0:
        if energy_balance(60) > 0:
            raise ValueError('Compatible static temperature lies below 60 K')
        temperature = brentq(energy_balance, 60, total_temperature, xtol=1e-11, rtol=1e-13)
    entropy_integral, error = quad(
        lambda t: float(model.state('PT', total_pressure, t)['cp'][0])/t,
        temperature, total_temperature, epsabs=1e-8, epsrel=2e-12)
    pressure = total_pressure*math.exp(-entropy_integral/R)
    value = {key: float(column[0]) for key, column in model.state('PT', pressure, temperature).items()}
    value.update(mach=mach, speed=mach*value['sound'], totalTemperature=total_temperature,
                 totalPressure=total_pressure, totalEnthalpy=h0,
                 entropyIntegral=entropy_integral, entropyQuadratureError=error)
    return value


def static_to_total(model, pressure, temperature, speed):
    """模块二。由静态和平均流速度反解总温/总压，用于往返校验和推导参考。"""
    if not math.isfinite(speed) or speed < 0:
        raise ValueError('Speed must be finite and nonnegative')
    value = model.state('PT', pressure, temperature)
    enthalpy = float(value['h'][0]) + .5*speed**2
    total_temperature = float(model.state('Ph', pressure, enthalpy)['T'][0])
    entropy_integral = quad(
        lambda t: float(model.state('PT', pressure, t)['cp'][0])/t,
        temperature, total_temperature, epsabs=1e-8, epsrel=2e-12)[0]
    return dict(totalTemperature=total_temperature,
                totalPressure=pressure*math.exp(entropy_integral/R), totalEnthalpy=enthalpy)
