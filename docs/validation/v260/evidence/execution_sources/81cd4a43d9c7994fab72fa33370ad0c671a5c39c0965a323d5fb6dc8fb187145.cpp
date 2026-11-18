// 模块职责：实际CNSVariable→气体虚接口的小状态重放；保存旧/新真正收到的相同输入。
#include "CSU2PhaseAwareAirStableV260.hpp"
#include "variables/CNSVariable.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstring>

using OldTrace260=CPhaseAirCaloricV260<false>;
using NewTrace260=CSU2PhaseAwareAirStableV260;
std::array<double,22> gasValues260(CFluidModel& q){return {q.GetTemperature(),q.GetPressure(),q.GetDensity(),q.GetStaticEnergy(),q.GetEnthalpy(),q.GetEntropy(),
  q.GetCp(),q.GetCv(),q.GetSoundSpeed(),q.GetLaminarViscosity(),q.GetThermalConductivity(),q.GetdPdrho_e(),q.GetdPde_rho(),q.GetdTde_rho(),
  q.Getdhdrho_P(),q.GetdhdP_rho(),q.Getdsdrho_P(),q.GetdsdP_rho(),q.GetdmudT_rho(),q.GetdktdT_rho(),q.Getdmudrho_T(),q.Getdktdrho_T()};}
void stateRow(std::ofstream& out,unsigned sample,const std::string& mode,CFluidModel& gas,double r,double e,const char* method){
  out<<sample<<','<<mode<<','<<r<<','<<e<<','<<gas.GetTemperature()<<','<<gas.GetPressure()<<','<<gas.GetDensity()<<','
    <<gas.GetStaticEnergy()<<','<<gas.GetEnthalpy()<<','<<gas.GetEntropy()<<','<<gas.GetCp()<<','<<gas.GetCv()<<','
    <<gas.GetLaminarViscosity()<<','<<gas.GetThermalConductivity()<<','<<method<<'\n';
}
int main(int argc,char** argv){
  std::string folder=argc>1?argv[1]:".";
  try{
    if(argc!=2)throw std::runtime_error("OUTPUT_REQUIRED");
    std::istringstream configText("SOLVER= RANS\nKIND_TURB_MODEL= SA\nSA_OPTIONS= (NEGATIVE, WITHFT2)\n"
      "SYSTEM_MEASUREMENTS= SI\nREF_DIMENSIONALIZATION= DIMENSIONAL\nFLUID_MODEL= RAH_PHASE_AIR\n"
      "VISCOSITY_MODEL= RAH_PHASE_LINEAR\nCONDUCTIVITY_MODEL= RAH_PHASE_LINEAR\n"
      "CONV_NUM_METHOD_FLOW= ROE\nTIME_DISCRE_FLOW= EULER_IMPLICIT\nTIME_DISCRE_TURB= EULER_IMPLICIT\n"
      "INIT_OPTION= TD_CONDITIONS\nMACH_NUMBER= 5\nFREESTREAM_PRESSURE= 1197.0279007383137\nFREESTREAM_TEMPERATURE= 226.50908361133003\n"
      "OUTPUT_FILES= (RESTART_ASCII)\n");
    CConfig cfg(configText,SU2_COMPONENT::SU2_CFD,false);OldTrace260 oldGas;NewTrace260 newGas;CSU2PhaseAwareAir originalGas;
    oldGas.SetLaminarViscosityModel(&cfg);newGas.SetLaminarViscosityModel(&cfg);oldGas.SetThermalConductivityModel(&cfg);newGas.SetThermalConductivityModel(&cfg);
    std::ifstream input(folder+"/states_input.csv");std::string line;std::getline(input,line);
    std::ofstream states(folder+"/states.csv");states<<std::setprecision(17)<<"sample,mode,receivedRho,receivedE,T,p,rho,e,h,s,Cp,Cv,mu,k,method\n";
    unsigned samples=0;std::vector<double> actualTargets;
    while(std::getline(input,line)){
      std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);unsigned sample,globalID,stage;std::array<double,6>a{};
      row>>sample>>globalID>>stage;for(auto& v:a)row>>v;if(!row||sample!=samples)throw std::runtime_error("INPUT_ROW");
      // 原保存物理能量→原生内部能量；两节点的Solution与OldSolution逐位同输入。
      a[4]+=rah_phase_air::shift*a[0];double velocity[3]={a[1]/a[0],a[2]/a[0],a[3]/a[0]};const double E=a[4]/a[0];
      CNSVariable oldNode(a[0],velocity,E,1,3,5,&cfg),newNode(a[0],velocity,E,1,3,5,&cfg),originalNode(a[0],velocity,E,1,3,5,&cfg);
      for(unsigned j=0;j<5;++j){oldNode.SetSolution(0,j,a[j]);oldNode.SetSolution_Old(0,j,a[j]);newNode.SetSolution(0,j,a[j]);newNode.SetSolution_Old(0,j,a[j]);originalNode.SetSolution(0,j,a[j]);originalNode.SetSolution_Old(0,j,a[j]);}
      oldGas.calls=0;newGas.calls=0;
      if(!oldNode.SetPrimVar(0,0.,0.,&oldGas)||!newNode.SetPrimVar(0,0.,0.,&newGas))throw std::runtime_error("VALID_REPLAY_PRIMITIVE_REFUSED");
      if(!originalNode.SetPrimVar(0,0.,0.,&originalGas))throw std::runtime_error("ORIGINAL_CLASS_REFUSED");
      const auto originalValues=gasValues260(originalGas),tracedValues=gasValues260(oldGas);
      if(std::memcmp(originalValues.data(),tracedValues.data(),sizeof(originalValues)))throw std::runtime_error("OLD_TRACE_DIFFERS_FROM_ORIGINAL_FINAL_CLASS");
      if(oldGas.calls!=1||newGas.calls!=1||std::memcmp(&oldGas.receivedRho,&newGas.receivedRho,8)||std::memcmp(&oldGas.receivedE,&newGas.receivedE,8))throw std::runtime_error("PRIMITIVE_TARGET_NOT_IDENTICAL");
      stateRow(states,sample,"old",oldGas,oldGas.receivedRho,oldGas.receivedE,"legacy-TOMS748");
      stateRow(states,sample,"new",newGas,newGas.receivedRho,newGas.receivedE,newGas.lastCaloricInverse.method);
      actualTargets.push_back(oldGas.receivedE);++samples;
    }
    states.close();if(!states||samples==0||samples>195)throw std::runtime_error("REPLAY_WRITE_OR_COUNT");
    // 模块二：节点/相邻可表示目标和非法输入。未知返回NaN，不填零。
    std::ofstream controls(folder+"/controls.csv");controls<<std::setprecision(17)<<"id,kind,label,sourceT,target,expectedValid,oldAccepted,newAccepted,oldT,newT,newMethod,newInterval\n";
    unsigned count=0,mismatches=0;
    auto control=[&](rah_phase_air::Quantity kind,const std::string& label,double sourceT,double target){
      const double lo=rah_phase_air_stable::caloricValue(rah_phase_air::Tmin,kind),hi=rah_phase_air_stable::caloricValue(rah_phase_air::Tmax,kind);
      const bool expected=std::isfinite(target)&&target>=lo&&target<=hi;double a=std::numeric_limits<double>::quiet_NaN(),b=a;bool okA=false,okB=false;rah_phase_air_stable::InverseInfo info;
      try{a=rah_phase_air::inverse(target,kind);okA=true;}catch(const std::domain_error&){}
      try{b=rah_phase_air_stable::inverse(target,kind,rah_phase_air::Pref,&info);okB=true;}catch(const std::domain_error&){}
      mismatches+=(okA!=expected)||(okB!=expected);
      controls<<count++<<','<<(kind==rah_phase_air::Quantity::E?"E":"H")<<','<<label<<','<<sourceT<<','<<target<<','<<expected<<','<<okA<<','<<okB<<','<<a<<','<<b<<','<<info.method<<','<<info.interval<<'\n';
    };
    for(auto kind:{rah_phase_air::Quantity::E,rah_phase_air::Quantity::H}){
      for(double t:{60.,80.,132.5,226.5,298.15,300.,1000.,1600.,2000.}){
        const double value=rah_phase_air_stable::caloricValue(t,kind);control(kind,"node",t,value);
        control(kind,"next-below",t,std::nextafter(value,-std::numeric_limits<double>::infinity()));
        control(kind,"next-above",t,std::nextafter(value,std::numeric_limits<double>::infinity()));
      }
      for(double target:{std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()})control(kind,"nonfinite",0,target);
    }
    controls.close();if(!controls)throw std::runtime_error("CONTROL_WRITE");
    // 模块三：真实H/S虚接口，S通路保持原实现；只取四个内点，不重跑七入口大回归。
    std::ofstream hs(folder+"/hs.csv");hs<<std::setprecision(17)<<"id,mode,targetH,targetS,inputP,T,p,h,s,method\n";unsigned index=0,sBitMismatches=0;
    for(double t:{226.50908361133003,300.123,1000.123,1600.123}){
      const double p=.05*rah_phase_air::ceiling(t);oldGas.SetTDState_PT(p,t);const double h=oldGas.GetEnthalpy(),s=oldGas.GetEntropy();
      oldGas.SetTDState_hs(h,s);newGas.SetTDState_hs(h,s);
      hs<<index<<",oldH,"<<h<<','<<s<<','<<p<<','<<oldGas.GetTemperature()<<','<<oldGas.GetPressure()<<','<<oldGas.GetEnthalpy()<<','<<oldGas.GetEntropy()<<",legacy-TOMS748\n";
      hs<<index<<",newH,"<<h<<','<<s<<','<<p<<','<<newGas.GetTemperature()<<','<<newGas.GetPressure()<<','<<newGas.GetEnthalpy()<<','<<newGas.GetEntropy()<<','<<newGas.lastCaloricInverse.method<<'\n';
      oldGas.SetTDState_Ps(p,s);newGas.SetTDState_Ps(p,s);double a=oldGas.GetTemperature(),b=newGas.GetTemperature();sBitMismatches+=std::memcmp(&a,&b,8)!=0;
      hs<<index<<",oldS,"<<h<<','<<s<<','<<p<<','<<a<<','<<oldGas.GetPressure()<<','<<oldGas.GetEnthalpy()<<','<<oldGas.GetEntropy()<<",legacy-entropy-TOMS748\n";
      hs<<index<<",newS,"<<h<<','<<s<<','<<p<<','<<b<<','<<newGas.GetPressure()<<','<<newGas.GetEnthalpy()<<','<<newGas.GetEntropy()<<",inherited-legacy-entropy-TOMS748\n";++index;
    }
    hs.close();if(!hs)throw std::runtime_error("HS_WRITE");
    std::ofstream finish(folder+"/native_completion.json");finish<<"{\"normalReturn\":true,\"samples\":"<<samples<<",\"stateRows\":"<<2*samples
      <<",\"controls\":"<<count<<",\"controlAcceptanceMismatches\":"<<mismatches<<",\"entropyTemperatureBitMismatches\":"<<sBitMismatches
      <<",\"primitiveInputsBitwiseMatched\":true,\"oldTracedMappingBitwiseOriginal22Fields\":true,\"longDoubleDigits\":"<<std::numeric_limits<long double>::digits
      <<",\"CFDExecuted\":false,\"physicalQualification\":false,\"newQualifiedAreaM2\":0}\n";finish.close();if(!finish)throw std::runtime_error("COMPLETION_WRITE");
    return 0;
  }catch(const std::exception& e){std::ofstream f(folder+"/native_failure.txt");f<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}
}
