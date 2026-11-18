// 模块职责：通过原CFluidModel虚接口接稳定E/H反算，其余气体/输运和拒绝沿用。
#pragma once
#include "CSU2PhaseAwareAir.hpp"
#include "StableCaloricInverseV260.hpp"

class CSU2PhaseAwareAirStableV260 : public CSU2PhaseAwareAir {
 public:
  rah_phase_air_stable::InverseInfo lastCaloricInverse;
  void SetTDState_rhoe(su2double rho,su2double e) override {
    rah_phase_air::require(std::isfinite(rho)&&std::isfinite(e),"RAH_PHASE_NONFINITE");
    const double t=rah_phase_air_stable::inverse(e,rah_phase_air::Quantity::E,rah_phase_air::Pref,&lastCaloricInverse);
    CSU2PhaseAwareAir::SetTDState_rhoT(rho,t);
  }
  void SetTDState_hs(su2double h,su2double entropy) override {
    rah_phase_air::require(std::isfinite(h)&&std::isfinite(entropy),"RAH_PHASE_NONFINITE");
    const double t=rah_phase_air_stable::inverse(h,rah_phase_air::Quantity::H,rah_phase_air::Pref,&lastCaloricInverse);
    const double p=rah_phase_air::Pref*std::exp((rah_phase_air::thermal(t).s-entropy)/rah_phase_air::R);
    CSU2PhaseAwareAir::SetTDState_PT(p,t);
  }
};
