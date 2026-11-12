// 模块职责：把具名气体的非法试探接到SU2原生旧状态回退；不放宽域或吞数值错误。
#pragma once
#include "fluid/CSU2PhaseAwareAir.hpp"
#include <string>
namespace rah_phase_recovery {
inline bool input_domain_error(const std::domain_error& error) {
  const std::string name=error.what();
  return name=="RAH_PHASE_NONFINITE" || name=="RAH_PHASE_DENSITY" ||
         name=="RAH_PHASE_INVERSE_INPUT" || name=="RAH_PHASE_ENERGY_ENTROPY_DOMAIN" ||
         name=="RAH_PHASE_T_DOMAIN" || name=="RAH_PHASE_PT_DOMAIN";
}
inline bool set_trial_state(CFluidModel* model,su2double density,su2double energy) {
  try {model->SetTDState_rhoe(density,energy);return true;}
  catch(const std::domain_error& error) {
    if(!dynamic_cast<CSU2PhaseAwareAir*>(model)||!input_domain_error(error))throw;
    return false; // 调用方进入既有Solution_Old回退，并返回RightVol=false。
  }
}
}
