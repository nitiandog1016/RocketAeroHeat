// 模块用途。直接调用真实CNSVariable和真实空气物性，验证回退后的完整状态与异常边界。
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2IdealAir.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <omp.h>

using State = std::array<double,5>;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

State conservative(double rho, double temperature, double k) {
  const auto air = rah_air::from_rho_temperature(rho, temperature);
  return {{rho, rho*500, -rho*20, rho*75, rho*(air.energy+0.5*(500*500+20*20+75*75)+k)}};
}

double error(double a, double b) { return std::abs(a-b)/std::max(1.0,std::abs(b)); }

double check(CNSVariable& node, CSU2IdealAir& air, const State& expected, double k, double mu_t) {
  double maximum = 0;
  for (unsigned i=0;i<5;++i) require(node.GetSolution(0,i)==expected[i], "conservative state changed");
  const double rho = expected[0];
  const double v2 = (expected[1]*expected[1]+expected[2]*expected[2]+expected[3]*expected[3])/(rho*rho);
  const auto state = rah_air::state("rhoe",rho,expected[4]/rho-0.5*v2-k);
  maximum = std::max({error(node.GetDensity(0),rho), error(node.GetPressure(0),state.pressure),
    error(node.GetTemperature(0),state.temperature), error(node.GetSoundSpeed(0),state.sound_speed),
    error(node.GetEnthalpy(0),expected[4]/rho+state.pressure/rho),
    error(node.GetLaminarViscosity(0),air.GetLaminarViscosity()),
    error(node.GetThermalConductivity(0),air.GetLaminarViscosity()*state.cp/.72),
    error(node.GetSpecificHeatCp(0),state.cp),error(node.GetEddyViscosity(0),mu_t)});
  require(node.GetFlowPrimitiveK(0)==k, "captured stage k differs");
  require(maximum < 2e-12, "recovered thermodynamics differs");
  return maximum;
}

class OtherFluid : public CFluidModel {
 public:
  void SetTDState_rhoe(double, double) override { throw std::domain_error("other-model-error"); }
};

int main(int argc,char** argv) {
  try {
    require(argc==3, "expected config and report paths");
    CConfig cfg(argv[1],SU2_COMPONENT::SU2_CFD);
    const double velocity[3] = {500,-20,75};
    const double k = 200, mu_t = 1e-4;
    const State old = conservative(.02,300,k), valid = conservative(.03,700,k);
    CNSVariable node(old[0],velocity,old[4]/old[0],1,3,5,&cfg);
    CSU2IdealAir air;
    double maximum = 0;
    unsigned recovered = 0, rejected = 0;
    for (const State& state : {old,valid}) {
      node.SetSolution(0,state.data()); node.SetSolution_Old(0,old.data());
      require(node.SetPrimVar(0,mu_t,k,&air), "legal state reported invalid");
      maximum=std::max(maximum,check(node,air,state,k,mu_t));
    }
    State low=valid, high=valid, zero=valid, negative=valid;
    const double kinetic=0.5*(500*500+20*20+75*75);
    low[4]=low[0]*(rah_air::from_rho_temperature(low[0],60).energy-10000+kinetic+k);
    high[4]=high[0]*(rah_air::from_rho_temperature(high[0],2000).energy+10000+kinetic+k);
    zero[0]=0; negative[0]=-.02;
    for (const State& trial : {low,high,zero,negative}) {
      air.SetTDState_PT(90000,1500); // 故意保留不同物性缓存，检验恢复确实重算。
      node.SetSolution(0,trial.data()); node.SetSolution_Old(0,old.data());
      require(!node.SetPrimVar(0,mu_t,k,&air), "recovery did not report invalid candidate");
      maximum=std::max(maximum,check(node,air,old,k,mu_t)); ++recovered;
    }
    // 旧能量原本对应300K，过大的当前k使恢复后的热内能低于60K，必须抛出。
    const double incompatible_k = k+rah_air::from_rho_temperature(.02,300).energy-
                                  rah_air::from_rho_temperature(.02,60).energy+1000;
    node.SetSolution(0,low.data());node.SetSolution_Old(0,old.data());
    try { node.SetPrimVar(0,mu_t,incompatible_k,&air); }
    catch (const std::domain_error&) { ++rejected; }
    require(rejected==1, "incompatible old state accepted");
    OtherFluid other;
    node.SetSolution(0,valid.data());node.SetSolution_Old(0,old.data());
    try { node.SetPrimVar(0,mu_t,k,&other); }
    catch (const std::domain_error& ex) { require(std::string(ex.what())=="other-model-error","wrong exception");++rejected; }
    require(rejected==2, "other model exception swallowed");
    for(unsigned i=0;i<5;++i) require(node.GetSolution(0,i)==valid[i],"other model rolled back");
    // 每个线程使用独立节点和物性实例，记录器共享锁保证行完整。
    int parallel_failed=0;
    #pragma omp parallel for num_threads(2) reduction(+:parallel_failed)
    for (int i=0;i<16;++i) {
      try {
        CNSVariable local(old[0],velocity,old[4]/old[0],1,3,5,&cfg);CSU2IdealAir fluid;
        local.SetSolution(0,high.data());local.SetSolution_Old(0,old.data());
        require(!local.SetPrimVar(0,mu_t,k,&fluid),"parallel recovery missing");
        check(local,fluid,old,k,mu_t);
      } catch (...) { ++parallel_failed; }
    }
    require(parallel_failed==0,"parallel state recovery failed");
    std::ofstream report(argv[2]);report<<std::setprecision(17);
    report<<"{\"passed\":true,\"legalCases\":2,\"recoveredCases\":"<<recovered
          <<",\"rejectedCases\":"<<rejected<<",\"parallelCases\":16,\"maximumError\":"<<maximum<<"}\n";
    return report ? 0 : 2;
  } catch (const std::exception& ex) { std::cerr<<ex.what()<<'\n'; return 1; }
}
