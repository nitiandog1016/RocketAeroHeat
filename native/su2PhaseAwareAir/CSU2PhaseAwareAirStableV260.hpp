// 模块职责：原CFluidModel完整接口的同模型适配；只对E/H选择稳定反算。
// 原CSU2PhaseAwareAir为final；本类直接实现CFluidModel，不改旧类或历史对象。
#pragma once
#include "CSU2PhaseAwareAir.hpp"
#include "StableCaloricInverseV260.hpp"

template<bool StableCaloric> class CPhaseAirCaloricV260 : public CFluidModel {
  // 此字段提交与旧类一致；只有完整联合域接受后才写入，禁止部分更新。
  void assign(const rah_phase_air::State& q) {
    Density=q.rho;Temperature=q.t;Pressure=q.p;StaticEnergy=q.e;Enthalpy=q.h;Entropy=q.s;
    Cp=q.cp;Cv=q.cv;SoundSpeed2=q.a2;dPdrho_e=q.dpdrho;dPde_rho=q.dpde;
    dTdrho_e=q.dtdrho;dTde_rho=q.dtde;dhdrho_P=q.dhdrho;dhdP_rho=q.dhdp;
    dsdrho_P=q.dsdrho;dsdP_rho=q.dsdp;Mu=q.mu;Kt=q.k;dmudrho_T=0;dktdrho_T=0;
    dmudT_rho=q.dmu;dktdT_rho=q.dk;
    static_cast<CPhaseAirViscosity*>(LaminarViscosity.get())->commit(q);
  }
 public:
  rah_phase_air_stable::InverseInfo lastCaloricInverse;
  double receivedRho=0,receivedE=0;unsigned calls=0;
  CPhaseAirCaloricV260() {
    static_assert(std::is_same<su2double,double>::value,"PhaseAir requires direct real double");
    LaminarViscosity.reset(new CPhaseAirViscosity());assign(rah_phase_air::PT(1000.,300.));
  }
  void SetTDState_PT(su2double p,su2double t) override {assign(rah_phase_air::state("PT",p,t));}
  void SetTDState_Prho(su2double p,su2double r) override {assign(rah_phase_air::state("Prho",p,r));}
  void SetTDState_rhoT(su2double r,su2double t) override {assign(rah_phase_air::state("rhoT",r,t));}
  void SetTDState_rhoe(su2double r,su2double e) override {
    receivedRho=r;receivedE=e;++calls;
    if constexpr(!StableCaloric) {assign(rah_phase_air::state("rhoe",r,e));}
    else {
      rah_phase_air::require(std::isfinite(r)&&std::isfinite(e),"RAH_PHASE_NONFINITE");
      const double t=rah_phase_air_stable::inverse(e,rah_phase_air::Quantity::E,rah_phase_air::Pref,&lastCaloricInverse);
      assign(rah_phase_air::rhoT(r,t));
    }
  }
  void SetTDState_hs(su2double h,su2double entropy) override {
    if constexpr(!StableCaloric) {assign(rah_phase_air::state("hs",h,entropy));}
    else {
      rah_phase_air::require(std::isfinite(h)&&std::isfinite(entropy),"RAH_PHASE_NONFINITE");
      const double t=rah_phase_air_stable::inverse(h,rah_phase_air::Quantity::H,rah_phase_air::Pref,&lastCaloricInverse);
      const double p=rah_phase_air::Pref*std::exp((rah_phase_air::thermal(t).s-entropy)/rah_phase_air::R);
      assign(rah_phase_air::PT(p,t));
    }
  }
  void SetTDState_Ps(su2double p,su2double s) override {assign(rah_phase_air::state("Ps",p,s));}
  void SetEnergy_Prho(su2double p,su2double r) override {SetTDState_Prho(p,r);}
  void ComputeDerivativeNRBC_Prho(su2double p,su2double r) override {SetTDState_Prho(p,r);}
  su2double GetLaminarViscosity() override {return Mu;}
  su2double GetThermalConductivity() override {return Kt;}
  void SetLaminarViscosityModel(const CConfig* c) override {
    if(c->GetKind_ViscosityModel()!=VISCOSITYMODEL::RAH_PHASE_LINEAR)throw std::invalid_argument("RAH_PHASE_REQUIRES_LINEAR_VISCOSITY");
  }
  void SetThermalConductivityModel(const CConfig* c) override {
    if(c->GetKind_ConductivityModel()!=CONDUCTIVITYMODEL::RAH_PHASE_LINEAR)throw std::invalid_argument("RAH_PHASE_REQUIRES_LINEAR_CONDUCTIVITY");
  }
};
// true为新复用类；false仅供本轮对照并须与未改原final类逐状态核对。
using CSU2PhaseAwareAirStableV260=CPhaseAirCaloricV260<true>;
