// 模块职责：将标准Air相态物性与显式能量坐标接入真实CFluidModel虚接口。
#pragma once
#include "CFluidModel.hpp"
#include "PhaseAwareAir.hpp"
#include <type_traits>

class CSU2PhaseAwareAir final : public CFluidModel {
  void assign(const rah_phase_air::State& q) {
    Density=q.rho;Temperature=q.t;Pressure=q.p;StaticEnergy=q.e;Enthalpy=q.h;Entropy=q.s;
    Cp=q.cp;Cv=q.cv;SoundSpeed2=q.a2;dPdrho_e=q.dpdrho;dPde_rho=q.dpde;
    dTdrho_e=q.dtdrho;dTde_rho=q.dtde;dhdrho_P=q.dhdrho;dhdP_rho=q.dhdp;
    dsdrho_P=q.dsdrho;dsdP_rho=q.dsdp;Mu=q.mu;Kt=q.k;dmudrho_T=0;dktdrho_T=0;
    dmudT_rho=q.dmu;dktdT_rho=q.dk;
  }
 public:
  CSU2PhaseAwareAir() {
    static_assert(std::is_same<su2double,double>::value,"PhaseAir requires direct real double");
    assign(rah_phase_air::PT(1000.,300.));
  }
  // 模块二：完整状态/导数入口，异常不提交字段，不使用旧OTH参考。
  void SetTDState_PT(su2double p,su2double t) override {assign(rah_phase_air::state("PT",p,t));}
  void SetTDState_Prho(su2double p,su2double rho) override {assign(rah_phase_air::state("Prho",p,rho));}
  void SetTDState_rhoT(su2double rho,su2double t) override {assign(rah_phase_air::state("rhoT",rho,t));}
  void SetTDState_rhoe(su2double rho,su2double e) override {assign(rah_phase_air::state("rhoe",rho,e));}
  void SetTDState_hs(su2double h,su2double s) override {assign(rah_phase_air::state("hs",h,s));}
  void SetTDState_Ps(su2double p,su2double s) override {assign(rah_phase_air::state("Ps",p,s));}
  void SetEnergy_Prho(su2double p,su2double rho) override {SetTDState_Prho(p,rho);}
  void ComputeDerivativeNRBC_Prho(su2double p,su2double rho) override {SetTDState_Prho(p,rho);}
  // 模块三：返回分子输运；原生数值层单独加入Cp*mu_t/Pr_t，禁止重复加入。
  su2double GetLaminarViscosity() override {return Mu;}
  su2double GetThermalConductivity() override {return Kt;}
  void SetLaminarViscosityModel(const CConfig* config) override {
    if(config->GetKind_ViscosityModel()!=VISCOSITYMODEL::RAH_PHASE_LINEAR)
      throw std::invalid_argument("RAH_PHASE_REQUIRES_LINEAR_VISCOSITY");
  }
  void SetThermalConductivityModel(const CConfig* config) override {
    if(config->GetKind_ConductivityModel()!=CONDUCTIVITYMODEL::RAH_PHASE_LINEAR)
      throw std::invalid_argument("RAH_PHASE_REQUIRES_LINEAR_CONDUCTIVITY");
  }
};
