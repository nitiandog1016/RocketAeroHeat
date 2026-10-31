// 独立状态查询工具：与受控SU2链接同一固定CoolProp库，不修改物性或裁剪状态。
// 输入PT/DE/HS三种状态，输出SI物性及偏导；用于初始化和真实守恒场的EOS审核。
#include "AbstractState.h"
#include "CoolProp.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
int main(int argc,char** argv) {
  try {
    if(argc!=3) throw std::runtime_error("usage: RAH_CoolProp_State input.csv output.csv");
    std::ifstream input(argv[1]); if(!input) throw std::runtime_error("cannot read input");
    std::string line;std::getline(input,line);if(!line.empty()&&line.back()=='\r')line.pop_back();
    if(line!="mode,a,b")throw std::runtime_error("expected mode,a,b header");
    std::unique_ptr<CoolProp::AbstractState> fluid(CoolProp::AbstractState::factory("HEOS","Air"));
    std::ofstream output(argv[2]);if(!output)throw std::runtime_error("cannot create output");
    output<<"row,p_Pa,T_K,rho_kg_m3,e_J_kg,h_J_kg,cp_J_kgK,cv_J_kgK,a_m_s,mu_Pa_s,k_W_mK,dpdrho_e,dpde_rho,dTdrho_e,dTde_rho,s_J_kgK,phase\n";
    output<<std::setprecision(17);size_t row=0;
    while(std::getline(input,line)) {
      if(line.empty())continue;
      std::replace(line.begin(),line.end(),',',' ');std::istringstream stream(line);
      std::string mode,extra;double a,b;if(!(stream>>mode>>a>>b)||(stream>>extra)||!std::isfinite(a)||!std::isfinite(b))throw std::runtime_error("invalid state row");
      if(mode=="PT")fluid->update(CoolProp::PT_INPUTS,a,b);
      else if(mode=="DE")fluid->update(CoolProp::DmassUmass_INPUTS,a,b);
      else if(mode=="HS")fluid->update(CoolProp::HmassSmass_INPUTS,a,b);
      else throw std::runtime_error("unsupported input mode");
      // EOS求值成功不等于位于数据库声明范围；拒绝外推和液相/两相，不能裁温后继续。
      if(fluid->T()<fluid->Tmin() || fluid->T()>fluid->Tmax() || fluid->p()>fluid->pmax())
        throw std::runtime_error("state outside declared HEOS Air EOS range");
      const auto phase=fluid->phase();
      if(phase!=CoolProp::iphase_gas && phase!=CoolProp::iphase_supercritical_gas && phase!=CoolProp::iphase_supercritical)
        throw std::runtime_error("state is not single-phase gas/supercritical Air");
      std::vector<double> v={fluid->p(),fluid->T(),fluid->rhomass(),fluid->umass(),fluid->hmass(),fluid->cpmass(),fluid->cvmass(),fluid->speed_sound(),
        fluid->viscosity(),fluid->conductivity(),
        fluid->first_partial_deriv(CoolProp::iP,CoolProp::iDmass,CoolProp::iUmass),fluid->first_partial_deriv(CoolProp::iP,CoolProp::iUmass,CoolProp::iDmass),
        fluid->first_partial_deriv(CoolProp::iT,CoolProp::iDmass,CoolProp::iUmass),fluid->first_partial_deriv(CoolProp::iT,CoolProp::iUmass,CoolProp::iDmass),fluid->smass()};
      for(double x:v)if(!std::isfinite(x))throw std::runtime_error("nonfinite property");
      if(v[0]<=0||v[1]<=0||v[2]<=0||v[5]<=0||v[6]<=0||v[7]<=0)throw std::runtime_error("invalid thermodynamic state");
      output<<row++;for(double x:v)output<<','<<x;output<<','<<int(fluid->phase())<<'\n';
    }
    if(!output||row==0)throw std::runtime_error("empty/incomplete output");
    std::cout<<"CoolProp "<<CoolProp::get_global_param_string("version")<<" commit "<<CoolProp::get_global_param_string("gitrevision")<<" rows "<<row<<'\n';
    return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
