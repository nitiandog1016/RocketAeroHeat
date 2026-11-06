// 模块用途。直接测试真实CNSVariable、真实物性和有限限步；覆盖温界、k阶段、二维及并行。
#include "variables/CNSVariable.hpp"
#include "variables/ConservativeStepTrace.hpp"
#include "fluid/CSU2IdealAir.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

using State=std::array<double,5>;
void require(bool ok,const char* message) { if(!ok) throw std::runtime_error(message); }
double relative(double a,double b) { return std::abs(a-b)/std::max(1.,std::abs(b)); }
State conservative(double rho,double temperature,double k,double u=500,double v=-20,double w=75) {
  auto air=rah_air::from_rho_temperature(rho,temperature);
  return {{rho,rho*u,rho*v,rho*w,rho*(air.energy+.5*(u*u+v*v+w*w)+k)}};
}
double energy(const State& s,double k) {
  return s[4]/s[0]-.5*(s[1]*s[1]+s[2]*s[2]+s[3]*s[3])/(s[0]*s[0])-k;
}
double check(CNSVariable& node,CSU2IdealAir& air,const State& old,const State& trial,double k,bool unchanged) {
  State result{};
  for(unsigned j=0;j<5;++j) result[j]=node.GetSolution(0,j);
  if(unchanged) require(result==trial,"legal conservative state changed");
  unsigned major=0;
  for(unsigned j=1;j<5;++j)
    if(std::abs(trial[j]-old[j])>std::abs(trial[major]-old[major])) major=j;
  double theta=unchanged?1.:(result[major]-old[major])/(trial[major]-old[major]);
  if(!unchanged) {
    require(theta>=-1e-14 && theta<=.5+1e-14,"step outside bounded segment");
    for(unsigned j=0;j<5;++j)
      require(relative(result[j],old[j]+theta*(trial[j]-old[j]))<2e-13,"not one conservative direction");
  }
  auto exact=rah_air::state("rhoe",result[0],energy(result,k));
  require(relative(node.GetTemperature(0),exact.temperature)<2e-12,"temperature mismatch");
  require(relative(node.GetPressure(0),exact.pressure)<2e-12,"pressure mismatch");
  require(relative(node.GetSoundSpeed(0),exact.sound_speed)<2e-12,"sound mismatch");
  require(relative(node.GetEnthalpy(0),result[4]/result[0]+exact.pressure/result[0])<2e-12,"enthalpy mismatch");
  require(relative(node.GetSpecificHeatCp(0),exact.cp)<2e-12,"cp mismatch");
  require(relative(node.GetThermalConductivity(0),air.GetLaminarViscosity()*exact.cp/.72)<2e-12,"conductivity mismatch");
  require(node.GetFlowPrimitiveK(0)==k,"actual stage k changed");
  require(node.GetEddyViscosity(0)==1e-4,"eddy viscosity changed");
  return theta;
}
class OtherFluid:public CFluidModel {
 public:void SetTDState_rhoe(double,double) override { throw std::domain_error("other-model-error"); }
};

int main(int argc,char** argv) {
  try {
    require(argc==4,"expected config, original events, and report");
    CConfig cfg(argv[1],SU2_COMPONENT::SU2_CFD);
    double velocity[3]={500,-20,75};
    const double k=200;
    const State old=conservative(.02,300,k);
    CNSVariable node(old[0],velocity,old[4]/old[0],1,3,5,&cfg);
    CSU2IdealAir air;
    // 夹具明确是合成上下文，不能作为真实网格定位证据。
    rah_step::ContextGuard context([](unsigned long point) {
      rah_step::Context c;c.known=true;c.global=100+point;c.domain=true;c.volume=.125;c.inner=7;
      c.coordinates={{1,2,3}};c.alpha=.75;return c;
    });
    unsigned legal=0, limited=0, rejected=0;
    for(auto trial:{old,conservative(.03,700,k),conservative(.02,60+1e-7,k),
                    conservative(.02,2000-1e-7,k),conservative(1,60,0,0,0,0),conservative(1,2000,0,0,0,0)}) {
      double used_k=(trial[0]==1)?0:k;
      node.SetSolution(0,trial.data());node.SetSolution_Old(0,old.data());
      require(node.SetPrimVar(0,1e-4,used_k,&air),"legal state rejected");
      check(node,air,old,trial,used_k,true);++legal;
    }
    State low=conservative(.03,60,k),high=conservative(.03,2000,k),zero=old,negative=old;
    low[4]-=300;high[4]+=300;zero[0]=0;negative[0]=-.02;
    for(auto trial:{low,high,zero,negative}) {
      air.SetTDState_PT(90000,1500);
      node.SetSolution(0,trial.data());node.SetSolution_Old(0,old.data());
      require(!node.SetPrimVar(0,1e-4,k,&air),"rejected trial reported legal");
      require(check(node,air,old,trial,k,false)>0,"positive segment not retained");++limited;
    }
    // 大速度差下上界不能用下界凸性代替；完整物性逐次检查。
    auto hot_old=conservative(.02,1900,k,3000),hot_trial=conservative(.03,2000,k,-3000);
    hot_trial[4]+=300;
    node.SetSolution(0,hot_trial.data());node.SetSolution_Old(0,hot_old.data());
    require(!node.SetPrimVar(0,1e-4,k,&air),"upper crossing not limited");
    require(check(node,air,hot_old,hot_trial,k,false)>0,"upper crossing lost progress");++limited;
    // 精确下界向外的方向不存在正步，32次后明确theta=0，不能伪造进展。
    auto boundary=conservative(1,60,0,0,0,0),outside=boundary;outside[4]-=10000;
    node.SetSolution(0,outside.data());node.SetSolution_Old(0,boundary.data());
    require(!node.SetPrimVar(0,1e-4,0,&air),"boundary no-step not reported");
    require(check(node,air,boundary,outside,0,false)==0,"boundary accepted outward step");
    // 当前k不同于构造旧态时的k，旧态不合法必须停止。
    const double bad_k=k+rah_air::from_rho_temperature(.02,300).energy-rah_air::from_rho_temperature(.02,60).energy+1000;
    node.SetSolution(0,low.data());node.SetSolution_Old(0,old.data());
    try { node.SetPrimVar(0,1e-4,bad_k,&air); }
    catch(const std::domain_error&) { ++rejected; }
    require(rejected==1,"invalid old state accepted");
    OtherFluid other;node.SetSolution(0,old.data());
    try {node.SetPrimVar(0,1e-4,k,&other);}
    catch(const std::domain_error& ex) {require(std::string(ex.what())=="other-model-error","wrong exception");++rejected;}
    require(rejected==2,"other model exception swallowed");
    for(unsigned j=0;j<5;++j) require(node.GetSolution(0,j)==old[j],"other model state mutated");
    // 逐条读取v164实测旧态/候选/当前k，独立原始证据不是新流场。
    std::ifstream events(argv[2]);std::string line;std::getline(events,line);unsigned real_pairs=0;
    while(std::getline(events,line)) {
      std::stringstream stream(line);std::vector<std::string> values;std::string token;
      while(std::getline(stream,token,',')) values.push_back(token);
      require(values.size()==18,"event CSV schema changed");
      const double actual_k=std::stod(values[4]);State trial{},prior{};
      for(unsigned j=0;j<5;++j) {trial[j]=std::stod(values[5+j]);prior[j]=std::stod(values[11+j]);}
      node.SetSolution(0,trial.data());node.SetSolution_Old(0,prior.data());
      require(!node.SetPrimVar(0,1e-4,actual_k,&air),"actual rejected pair accepted directly");
      require(check(node,air,prior,trial,actual_k,false)>0,"actual pair cannot retain positive step");++real_pairs;
    }
    require(real_pairs==14,"missing actual event pairs");
    // 二维能量列和三维记录中的零w必须一致。
    auto old2=conservative(.02,300,k,500,-20,0),low2=conservative(.03,60,k,500,-20,0);low2[4]-=300;
    CNSVariable node2(old2[0],velocity,old2[4]/old2[0],1,2,4,&cfg);
    double u2[4]={low2[0],low2[1],low2[2],low2[4]},o2[4]={old2[0],old2[1],old2[2],old2[4]};
    node2.SetSolution(0,u2);node2.SetSolution_Old(0,o2);
    require(!node2.SetPrimVar(0,1e-4,k,&air),"2D limit missing");
    require(node2.GetTemperature(0)>=60 && node2.GetTemperature(0)<=2000,"2D temperature invalid");
    std::vector<std::unique_ptr<CNSVariable>> nodes;std::vector<std::unique_ptr<CSU2IdealAir>> fluids;
    for(int i=0;i<16;++i) {nodes.emplace_back(new CNSVariable(old[0],velocity,old[4]/old[0],1,3,5,&cfg));fluids.emplace_back(new CSU2IdealAir());}
    int parallel_failed=0;
    #pragma omp parallel for num_threads(2) reduction(+:parallel_failed)
    for(int i=0;i<16;++i) {
      try {
        auto& n=*nodes[i];auto& f=*fluids[i];
        n.SetSolution(0,high.data());n.SetSolution_Old(0,old.data());
        require(!n.SetPrimVar(0,1e-4,k,&f),"parallel missing limit");
        require(check(n,f,old,high,k,false)>0,"parallel lost direction");
      } catch(...) {++parallel_failed;}
    }
    require(parallel_failed==0,"parallel test failed");
    std::ofstream report(argv[3]);report<<"{\"passed\":true,\"legalCases\":"<<legal<<",\"limitedCases\":"<<limited
      <<",\"zeroStepCases\":1,\"rejectedCases\":"<<rejected<<",\"realPairs\":"<<real_pairs<<",\"twoDimensionalCases\":1,\"parallelCases\":16}\n";
    return report?0:2;
  } catch(const std::exception& ex) {std::cerr<<ex.what()<<'\n';return 1;}
}
