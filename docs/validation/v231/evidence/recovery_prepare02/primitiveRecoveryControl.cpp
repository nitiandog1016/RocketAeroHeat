// 模块一：在真实Euler变量上比较非法试探处理，不推进流场。
#include "drivers/CSinglezoneDriver.hpp"
#include "variables/CEulerVariable.hpp"
#include "PhaseAirPrimitiveRecovery.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <array>
#include <limits>

class PrimitiveRecoveryControl final : public CSinglezoneDriver {
 public:
  using CSinglezoneDriver::CSinglezoneDriver;
  void check(const std::string& folder,const std::string& input) {
    auto* geo=geometry_container[0][0][0];auto* flow=solver_container[0][0][0][FLOW_SOL];
    auto* node=dynamic_cast<CEulerVariable*>(flow->GetNodes());auto* gas=flow->GetFluidModel();
    if(geo->GetnPoint()!=24||!node||!dynamic_cast<CSU2PhaseAwareAir*>(gas)||flow->GetnVar()!=5)
      throw std::runtime_error("RECOVERY_CONTROL_IDENTITY");
    std::ifstream in(input);std::string line;std::getline(in,line);std::array<double,5> old{},bad{};
    for(auto* values:{&old,&bad}){std::getline(in,line);std::istringstream row(line);
      for(int j=0;j<5;++j){std::string value;std::getline(row,value,',');(*values)[j]=std::stod(value);}}
    if(!in)throw std::runtime_error("RECOVERY_INPUT_READ");
    auto state=[](double p,double t){const auto q=rah_phase_air::PT(p,t);const double u=600.,v=-20.,w=90.;
      return std::array<double,5>{q.rho,q.rho*u,q.rho*v,q.rho*w,q.rho*(q.e+.5*(u*u+v*v+w*w))};};
    const auto cold=state(100.,80.); // 物理内能可为负，内部坐标仍相容。
    auto high=old;high[4]=old[0]*1e8;auto density=old;density[0]=-1.;
    auto nan=old;nan[4]=std::numeric_limits<double>::quiet_NaN();
    auto phase=cold;for(auto& value:phase)value*=1e5; // T合法但压力超联合域。
    std::ofstream out(folder+"/primitive.csv");out<<std::setprecision(17)<<"case,returned,exception";
    for(int j=0;j<5;++j)out<<",U"<<j;out<<",p,T,dpdrho,dpde\n";
    auto run=[&](const char* name,const std::array<double,5>& prior,const std::array<double,5>& trial) {
      node->SetSolution(0,old.data());node->SetSolution_Old(0,old.data());node->SetPrimVar(0,gas);
      node->SetSolution_Old(0,prior.data());node->SetSolution(0,trial.data());
      int physical=-1;std::string error="none";
      try{physical=node->SetPrimVar(0,gas);node->SetSecondaryVar(0,gas);}
      catch(const std::domain_error& e){error=e.what();}
      out<<name<<','<<physical<<','<<error;
      for(int j=0;j<5;++j)out<<','<<node->GetSolution(0,j);
      if(error=="none")out<<','<<node->GetPressure(0)<<','<<node->GetTemperature(0)<<','<<node->GetSecondary(0)[0]<<','<<node->GetSecondary(0)[1];
      else out<<",unknown,unknown,unknown,unknown";out<<'\n';
    };
    run("valid",old,old);run("cold_negative_physical_e",old,cold);
    run("actual_raw_invalid",old,bad);run("high_energy",old,high);run("negative_density",old,density);
    run("nonfinite",old,nan);run("joint_domain",old,phase);run("invalid_old",bad,bad);
    // 模块二：其他气体与非域数值错误不得变成可恢复的物性拒绝。
    class OtherGas final : public CFluidModel {void SetTDState_rhoe(su2double,su2double) override {throw std::domain_error("RAH_PHASE_ENERGY_ENTROPY_DOMAIN");}} other;
    bool propagated=false;try{rah_phase_recovery::set_trial_state(&other,1.,1.);}catch(const std::domain_error&){propagated=true;}
    if(!propagated||rah_phase_recovery::input_domain_error(std::domain_error("RAH_PHASE_INVERSE_CONVERGENCE"))||
       rah_phase_recovery::input_domain_error(std::domain_error("RAH_PHASE_NUMERICAL_RANGE")))throw std::runtime_error("ERROR_CLASSIFICATION");
    std::ofstream result(folder+"/control.json");result<<"{\"passed\":true,\"cases\":8,\"otherGasPropagated\":true,\"numericalErrorsNotRecovered\":true,\"flowUpdates\":0}\n";
    if(!out||!result)throw std::runtime_error("RECOVERY_WRITE");
  }
};
int main(int argc,char** argv){try{
  SU2_MPI::Init(&argc,&argv);if(argc!=4)throw std::runtime_error("CONFIG_OUTPUT_INPUT_REQUIRED");
  PrimitiveRecoveryControl driver(argv[1],1,SU2_MPI::GetComm());driver.check(argv[2],argv[3]);
  driver.Finalize();SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
