// 模块用途。集中记录SST总能量与热状态之间的换算，供原生调用点和差分测试共用。
#pragma once
#include "CFluidModel.hpp"

namespace rah_air {
inline double thermal_energy(double total_energy,double velocity_squared,double turbulent_k) {
  return total_energy-.5*velocity_squared-turbulent_k;
}

// 冻结比湍动能k，密度/动量/总能量密度作为守恒自变量。
inline double temperature_density_derivative(const CFluidModel& fluid,
    double total_energy,double velocity_squared) {
  return fluid.GetdTdrho_e()+fluid.GetdTde_rho()*(velocity_squared-total_energy)/fluid.GetDensity();
}

inline double temperature_energy_derivative(const CFluidModel& fluid) {
  return fluid.GetdTde_rho()/fluid.GetDensity();
}

inline double conservative_chi(const CFluidModel& fluid,double turbulent_k) {
  const double kappa=fluid.GetdPde_rho()/fluid.GetDensity();
  return fluid.GetdPdrho_e()-kappa*(fluid.GetStaticEnergy()+turbulent_k);
}
} // namespace rah_air
