// 模块用途。把已验证理想空气核接入真实SU2 CFluidModel，全部输入使用SI。
// 编译包含目录须提供固定SU2的fluid目录。本类只允许直接双精度求解。
#pragma once
#include "CFluidModel.hpp"
#include "CSutherland.hpp"
#include "IdealAir.hpp"
#include "PropertyQueryContext.hpp"
#include <type_traits>

class CSU2IdealAir final : public CFluidModel {
 private:
  double prandtl_ = .72;
  double dcp_dt_ = 0;

  // 模块一。先得到完整合法状态，随后一次提交到基类字段。
  void assign(const rah_air::State& value) {
    Density=value.rho; Temperature=value.temperature; Pressure=value.pressure;
    StaticEnergy=value.energy; Enthalpy=value.enthalpy; Entropy=value.entropy;
    Cp=value.cp; Cv=value.cv; SoundSpeed2=value.sound_speed*value.sound_speed;
    dPdrho_e=value.dp_drho_e; dPde_rho=value.dp_de_rho;
    dTdrho_e=value.dt_drho_e; dTde_rho=value.dt_de_rho;
    dhdrho_P=value.dh_drho_p; dhdP_rho=value.dh_dp_rho;
    dsdrho_P=value.ds_drho_p; dsdP_rho=value.ds_dp_rho;
    dcp_dt_=value.dcp_dt;
  }

 public:
  CSU2IdealAir() {
    static_assert(std::is_same<su2double,double>::value, "IdealAir requires direct double precision");
    LaminarViscosity.reset(new CSutherland(1.716e-5,273.15,110.4));
    assign(rah_air::from_pressure_temperature(101325,300));
  }

  // 模块二。保留原生虚接口，避免调用者读到未赋值的导数或旧Cp。
  void SetTDState_rhoe(su2double rho,su2double e) override { assign(rah_query::evaluate("rhoe",rho,e)); }
  void SetTDState_PT(su2double p,su2double t) override { assign(rah_query::evaluate("PT",p,t)); }
  void SetTDState_Prho(su2double p,su2double rho) override { assign(rah_query::evaluate("Prho",p,rho)); }
  void SetTDState_rhoT(su2double rho,su2double t) override { assign(rah_query::evaluate("rhoT",rho,t)); }
  void SetTDState_hs(su2double h,su2double s) override { assign(rah_query::evaluate("hs",h,s)); }
  void SetTDState_Ps(su2double p,su2double s) override { assign(rah_query::evaluate("Ps",p,s)); }
  void SetEnergy_Prho(su2double p,su2double rho) override { SetTDState_Prho(p,rho); }
  void ComputeDerivativeNRBC_Prho(su2double p,su2double rho) override { SetTDState_Prho(p,rho); }

  // 模块三。复用SU2的Sutherland黏度，常Pr导热偏导增加cp'(T)项。
  // 原生可压缩数值层另外加入湍流导热，本接口只返回分子导热。
  void SetLaminarViscosityModel(const CConfig* config) override {
    if (config->GetKind_ViscosityModel()!=VISCOSITYMODEL::SUTHERLAND)
      throw std::invalid_argument("IDEAL_AIR_VARIABLE_CP requires SUTHERLAND");
    CFluidModel::SetLaminarViscosityModel(config);
  }
  void SetThermalConductivityModel(const CConfig* config) override {
    if (config->GetKind_ConductivityModel()!=CONDUCTIVITYMODEL::CONSTANT_PRANDTL)
      throw std::invalid_argument("IDEAL_AIR_VARIABLE_CP requires CONSTANT_PRANDTL");
    prandtl_=config->GetPrandtl_Lam();
    rah_air::require(std::isfinite(prandtl_) && prandtl_>0,"invalid laminar Prandtl number");
  }
  su2double GetThermalConductivity() override {
    GetLaminarViscosity();
    Kt=Mu*Cp/prandtl_;
    dktdrho_T=dmudrho_T*Cp/prandtl_;
    dktdT_rho=(dmudT_rho*Cp+Mu*dcp_dt_)/prandtl_;
    return Kt;
  }
};
