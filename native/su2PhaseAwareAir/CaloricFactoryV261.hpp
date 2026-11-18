// 模块职责：有量纲初始化和各线程使用同一工厂，原生气体枚举继续标识同一种相态空气。
#pragma once
#include "CSU2PhaseAwareAirStableV260.hpp"

namespace rah_caloric_factory {
inline CFluidModel* create(const CConfig& config) {
  if (config.OptionIsSet("RAH_CALORIC_INVERSE")) return new CSU2PhaseAwareAirStableV260();
  return new CSU2PhaseAwareAir();
}

inline const char* method(const CFluidModel* gas) {
  if (dynamic_cast<const CSU2PhaseAwareAirStableV260*>(gas)) return "STABLE_PIECEWISE";
  if (dynamic_cast<const CSU2PhaseAwareAir*>(gas)) return "TOMS748";
  return "OTHER_GAS";
}

inline bool isPhaseAir(const CFluidModel* gas) {
  return dynamic_cast<const CSU2PhaseAwareAir*>(gas) || dynamic_cast<const CSU2PhaseAwareAirStableV260*>(gas);
}
}
