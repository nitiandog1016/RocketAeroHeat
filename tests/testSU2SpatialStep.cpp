// 模块用途：真实CNS变量/真实物性下检验层流同向限步和只读查询上下文；不运行CFD。
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2IdealAir.hpp"
#include "fluid/PropertyQueryContext.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

using State = std::array<double, 5>;
void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
double scaled(double actual, double expected) { return std::abs(actual-expected)/std::max(1., std::abs(expected)); }

State conservative(double rho, double temperature, double u=500., double v=-20., double w=75.) {
  const auto air = rah_air::from_rho_temperature(rho, temperature);
  return {{rho, rho*u, rho*v, rho*w, rho*(air.energy+.5*(u*u+v*v+w*w))}};
}
double thermal(const State& state) {
  return state[4]/state[0]-.5*(state[1]*state[1]+state[2]*state[2]+state[3]*state[3])/(state[0]*state[0]);
}
State conserved(CNSVariable& node) {
  State result{}; for (unsigned j=0; j<5; ++j) result[j]=node.GetSolution(0,j); return result;
}
bool same_context(const rah_query::Context& a, const rah_query::Context& b) {
  return std::string(a.phase)==std::string(b.phase) && a.localKnown==b.localKnown && a.globalKnown==b.globalKnown &&
    a.markerKnown==b.markerKnown && a.innerKnown==b.innerKnown && a.local==b.local && a.global==b.global &&
    a.marker==b.marker && a.inner==b.inner;
}

double audit_state(CNSVariable& node, CSU2IdealAir& air, const State& old, const State& trial, bool unchanged) {
  // 模块一：五个守恒量必须共用theta；从接收态重算全部热物性，不读旧缓存当真值。
  const State accepted=conserved(node);
  if (unchanged) require(accepted==trial, "legal laminar conservative state changed");
  unsigned major=0;
  for (unsigned j=1; j<5; ++j)
    if (std::abs(trial[j]-old[j])>std::abs(trial[major]-old[major])) major=j;
  const double theta=unchanged?1.:(accepted[major]-old[major])/(trial[major]-old[major]);
  if (!unchanged) {
    require(theta>=-1e-14 && theta<=.5+1e-14, "accepted state outside declared backtracking segment");
    for (unsigned j=0; j<5; ++j)
      require(scaled(accepted[j],old[j]+theta*(trial[j]-old[j]))<2e-13, "not one common conservative direction");
  }
  const auto exact=rah_air::state("rhoe",accepted[0],thermal(accepted));
  require(exact.temperature>=60. && exact.temperature<=2000., "accepted state outside original thermal domain");
  require(scaled(node.GetTemperature(0),exact.temperature)<2e-12, "accepted temperature mismatch");
  require(scaled(node.GetPressure(0),exact.pressure)<2e-12, "accepted pressure mismatch");
  require(scaled(node.GetSoundSpeed(0),exact.sound_speed)<2e-12, "accepted sound speed mismatch");
  require(scaled(node.GetEnthalpy(0),accepted[4]/accepted[0]+exact.pressure/accepted[0])<2e-12, "accepted enthalpy mismatch");
  require(scaled(node.GetSpecificHeatCp(0),exact.cp)<2e-12, "accepted cp mismatch");
  require(scaled(node.GetThermalConductivity(0),air.GetLaminarViscosity()*exact.cp/.72)<2e-12, "accepted conductivity mismatch");
  require(node.GetFlowPrimitiveK(0)==0., "laminar state acquired turbulent kinetic energy");
  require(node.GetEddyViscosity(0)==0., "laminar state acquired eddy viscosity");
  return theta;
}

class OtherDomainFluid: public CFluidModel {
 public: void SetTDState_rhoe(double,double) override { throw std::domain_error("other-model-domain-error"); }
};
class OtherRuntimeFluid: public CFluidModel {
 public: void SetTDState_rhoe(double,double) override { throw std::runtime_error("other-model-runtime-error"); }
};

void audit_property_context() {
  // 模块二：合法查询所有18字段逐值保持；异常只增加身份并继续抛出。
  const auto original=rah_query::context();
  rah_query::Context outer; outer.phase="fixture.outer"; outer.localKnown=true; outer.local=7; outer.innerKnown=true; outer.inner=11;
  {
    rah_query::Guard guard(outer);
    auto expected=rah_air::state("PT",90000.,300.);
    auto actual=rah_query::evaluate("PT",90000.,300.);
    for (auto field:{&rah_air::State::rho,&rah_air::State::temperature,&rah_air::State::pressure,
        &rah_air::State::energy,&rah_air::State::enthalpy,&rah_air::State::entropy,&rah_air::State::cp,
        &rah_air::State::cv,&rah_air::State::sound_speed,&rah_air::State::dp_drho_e,&rah_air::State::dp_de_rho,
        &rah_air::State::dt_drho_e,&rah_air::State::dt_de_rho,&rah_air::State::dh_drho_p,&rah_air::State::dh_dp_rho,
        &rah_air::State::ds_drho_p,&rah_air::State::ds_dp_rho,&rah_air::State::dcp_dt})
      require(actual.*field==expected.*field, "context wrapper changed a legal property value");
    rah_query::Context inner=outer; inner.phase="fixture.riemann"; inner.globalKnown=true; inner.global=101;
    inner.markerKnown=true; inner.marker=3;
    bool caught=false;
    try {
      rah_query::Guard nested(inner); nested.phase("riemann.exterior");
      rah_query::evaluate("Prho",-100.,.02);
    } catch (const std::domain_error& error) {
      const std::string message=error.what();
      require(message.find("RAH_PROPERTY_QUERY")!=std::string::npos && message.find("mode=Prho")!=std::string::npos &&
          message.find("phase=riemann.exterior")!=std::string::npos && message.find("global=101")!=std::string::npos &&
          message.find("marker=3")!=std::string::npos && message.find("input1=")!=std::string::npos,
          "boundary query lost its explicit failure context");
      caught=true;
    }
    require(caught, "invalid boundary pressure query was recovered or clipped");
    require(same_context(rah_query::context(),outer), "nested query context leaked after exception");
  }
  require(same_context(rah_query::context(),original), "outer query context was not restored");
}

int main(int argc,char** argv) {
  try {
    require(argc==3, "expected spatial-laminar config and report path");
    std::ifstream stream(argv[1]); require(bool(stream), "configuration unavailable");
    const std::string config_text((std::istreambuf_iterator<char>(stream)),std::istreambuf_iterator<char>());
    require(config_text.find("MARKER_PYTHON_CUSTOM")!=std::string::npos &&
        config_text.find("MARKER_ISOTHERMAL")!=std::string::npos, "fixture must retain spatial-temperature wall configuration");
    CConfig config(argv[1],SU2_COMPONENT::SU2_CFD);
    require(config.GetKind_Solver()==MAIN_SOLVER::NAVIER_STOKES && config.GetKind_Turb_Model()==TURB_MODEL::NONE,
        "fixture is not the declared laminar Navier-Stokes configuration");
    require(config.GetnMGLevels()==0 && config.GetPrandtl_Lam()==.72, "fixture changed wall/grid/transport contract");
    // 这里只核对空间壁温配置相容；全壁节点实际施加和热流须由随后唯一小CFD另验。
    double velocity[3]={500.,-20.,75.}; const State old=conservative(.02,300.);
    CNSVariable node(old[0],velocity,old[4]/old[0],1,3,5,&config); CSU2IdealAir air;
    unsigned legal=0,limited=0,zero_steps=0,rejected_old=0,other_exceptions=0;
    for (const auto trial:{old,conservative(.03,700.),conservative(.02,60.+1e-6),
        conservative(.02,2000.-1e-6),conservative(1.,60.,0.,0.,0.),conservative(1.,2000.,0.,0.,0.)}) {
      air.SetTDState_PT(90000.,1500.);
      node.SetSolution(0,trial.data()); node.SetSolution_Old(0,old.data());
      require(node.SetPrimVar(0,0.,0.,&air), "legal laminar trial rejected");
      audit_state(node,air,old,trial,true); ++legal;
    }
    State low=conservative(.03,60.),high=conservative(.03,2000.),zero_density=old,negative_density=old;
    low[4]-=300.; high[4]+=300.; zero_density[0]=0.; negative_density[0]=-.02;
    for (const auto trial:{low,high,zero_density,negative_density}) {
      air.SetTDState_PT(90000.,1500.);
      node.SetSolution(0,trial.data()); node.SetSolution_Old(0,old.data());
      require(!node.SetPrimVar(0,0.,0.,&air), "inadmissible laminar trial reported unchanged");
      require(audit_state(node,air,old,trial,false)>0., "valid old state lost available positive step"); ++limited;
    }
    // 模块三：大反向速度使中间混合态越过上界，不能把下界凸性套到整个温区。
    const auto hot_old=conservative(.02,1900.,3000.); auto hot_trial=conservative(.03,2000.,-3000.); hot_trial[4]+=300.;
    node.SetSolution(0,hot_trial.data()); node.SetSolution_Old(0,hot_old.data());
    require(!node.SetPrimVar(0,0.,0.,&air), "upper-bound crossing was not limited");
    const double hot_theta=audit_state(node,air,hot_old,hot_trial,false);
    require(hot_theta>0. && hot_theta<.5, "nonconvex upper-bound example did not reject the half-step"); ++limited;
    for (double temperature:{60.,2000.}) {
      const auto boundary=conservative(1.,temperature,0.,0.,0.); auto outward=boundary;
      outward[4]+=(temperature==60.?-10000.:10000.);
      node.SetSolution(0,outward.data()); node.SetSolution_Old(0,boundary.data());
      require(!node.SetPrimVar(0,0.,0.,&air), "outward boundary direction not rejected");
      require(audit_state(node,air,boundary,outward,false)==0. && conserved(node)==boundary,
          "zero-step case invented progress or changed old state"); ++zero_steps;
    }
    for (const auto invalid_old:{low,negative_density}) {
      node.SetSolution(0,low.data()); node.SetSolution_Old(0,invalid_old.data());
      bool caught=false;
      try { node.SetPrimVar(0,0.,0.,&air); }
      catch (const std::domain_error& error) {
        const std::string message=error.what();
        require(message.find("mode=rhoe")!=std::string::npos && message.find("phase=node.old")!=std::string::npos,
            "invalid old state lost its property query phase"); caught=true;
      }
      require(caught, "invalid old state silently accepted"); ++rejected_old;
    }
    // 其他模型/非domain_error必须原样透传，不能扩成吞异常的通用恢复。
    OtherDomainFluid other_domain; OtherRuntimeFluid other_runtime;
    node.SetSolution(0,old.data()); node.SetSolution_Old(0,old.data());
    try { node.SetPrimVar(0,0.,0.,&other_domain); }
    catch (const std::domain_error& error) {
      require(std::string(error.what())=="other-model-domain-error", "other model domain exception changed"); ++other_exceptions;
    }
    require(other_exceptions==1 && conserved(node)==old, "other model exception swallowed or conservative state changed");
    try { node.SetPrimVar(0,0.,0.,&other_runtime); }
    catch (const std::runtime_error& error) {
      require(std::string(error.what())=="other-model-runtime-error", "non-domain exception changed"); ++other_exceptions;
    }
    require(other_exceptions==2 && conserved(node)==old, "non-domain exception swallowed or conservative state changed");
    audit_property_context();

    // 模块四：独立变量/物性对象，两个线程分别安装上下文并确认退出后恢复。
    std::vector<std::unique_ptr<CNSVariable>> nodes; std::vector<std::unique_ptr<CSU2IdealAir>> fluids;
    for (unsigned i=0;i<16;++i) {
      nodes.emplace_back(new CNSVariable(old[0],velocity,old[4]/old[0],1,3,5,&config));
      fluids.emplace_back(new CSU2IdealAir());
    }
    int parallel_failures=0;
    #pragma omp parallel for num_threads(2) reduction(+:parallel_failures)
    for (int i=0;i<16;++i) {
      try {
        const auto before=rah_query::context();
        {
          rah_query::Context context; context.phase="fixture.parallel"; context.globalKnown=true; context.global=1000+i;
          rah_query::Guard guard(context);
          auto& n=*nodes[i]; auto& fluid=*fluids[i]; n.SetSolution(0,high.data()); n.SetSolution_Old(0,old.data());
          require(!n.SetPrimVar(0,0.,0.,&fluid), "parallel limit absent");
          require(audit_state(n,fluid,old,high,false)>0., "parallel positive step lost");
          require(same_context(rah_query::context(),context), "CNS query context leaked into caller");
          bool caught=false;
          try { rah_query::evaluate("Prho",-10.,.02); }
          catch (const std::domain_error& error) {
            require(std::string(error.what()).find("global="+std::to_string(1000+i))!=std::string::npos,
                "thread-local query context mixed global identities"); caught=true;
          }
          require(caught, "parallel invalid boundary query accepted");
        }
        require(same_context(rah_query::context(),before), "parallel guard failed to restore prior context");
      } catch (...) { ++parallel_failures; }
    }
    require(parallel_failures==0, "parallel state/context test failed");
    std::ofstream report(argv[2]); report<<std::setprecision(17);
    report<<"{\"passed\":true,\"legalCases\":"<<legal<<",\"limitedCases\":"<<limited
      <<",\"zeroStepCases\":"<<zero_steps<<",\"invalidOldCases\":"<<rejected_old
      <<",\"otherExceptionCases\":"<<other_exceptions<<",\"parallelCases\":16,\"propertyLegalFieldsExactlyEqual\":18"
      <<",\"upperNonconvexAcceptedTheta\":"<<hot_theta
      <<",\"k\":0,\"eddyViscosity\":0,\"sourceSSTEventPairsUsed\":0,\"CFDLaunches\":0"
      <<",\"spatialWallConfigurationChecked\":true,\"spatialWallActualApplicationValidated\":false"
      <<",\"boundaryQueryRecoveryIntroduced\":false,\"physicalQualified\":false}\n";
    require(bool(report), "cannot write spatial-step fixture report");
    return 0;
  } catch (const std::exception& error) {
    std::cerr<<error.what()<<'\n'; return 1;
  }
}
