// 模块用途。使用真实SU2基类检验物性适配、冻结k守恒导数与准确原生投影Jacobian。
#include "../third_party/ideal_air/CSU2IdealAir.hpp"
#include "../third_party/ideal_air/SU2EnergyState.hpp"
#include "NativeGeneralJacobian.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>

using Vector=std::array<double,5>;

double relative(double actual,double expected) {
  return std::abs(actual-expected)/std::max(1.,std::abs(expected));
}

// 独立守恒状态恢复，作为原生矩阵的差分参照。
Vector flux(const Vector& u,double k,const std::array<double,3>& normal,double* temperature=nullptr) {
  const double velocity_squared=(u[1]*u[1]+u[2]*u[2]+u[3]*u[3])/(u[0]*u[0]);
  CSU2IdealAir fluid;
  fluid.SetTDState_rhoe(u[0],u[4]/u[0]-.5*velocity_squared-k);
  if (temperature) *temperature=fluid.GetTemperature();
  const double mass=u[1]*normal[0]+u[2]*normal[1]+u[3]*normal[2];
  return {{mass,mass*u[1]/u[0]+fluid.GetPressure()*normal[0],
    mass*u[2]/u[0]+fluid.GetPressure()*normal[1],
    mass*u[3]/u[0]+fluid.GetPressure()*normal[2],mass*(u[4]+fluid.GetPressure())/u[0]}};
}

int main(int argc,char** argv) {
  try {
    if (argc!=2) throw std::runtime_error("expected output CSV");
    std::ofstream out(argv[1]); out<<std::setprecision(17);
    out<<"T,rho,k,state_error,sound_error,jacobian_scaled_error,temperature_derivative_error,conductivity_derivative_error,heat_frozen_coefficient_error,atomic_rejections,thread_error\n";
    unsigned cases=0;
    for (double t : {80.,300.,1000.,1900.}) for (double rho : {.001,.02,1.}) for (double k : {0.,200.,10000.}) {
      CSU2IdealAir fluid;
      fluid.SetTDState_rhoT(rho,t);
      const double p=fluid.GetPressure(),h=fluid.GetEnthalpy(),s=fluid.GetEntropy(),e=fluid.GetStaticEnergy();
      double state_error=0;
      auto check=[&]() {
        state_error=std::max({state_error,relative(fluid.GetTemperature(),t),relative(fluid.GetDensity(),rho),
            relative(fluid.GetPressure(),p),relative(fluid.GetStaticEnergy(),e),relative(fluid.GetEnthalpy(),h)});
      };
      CFluidModel& base=fluid;
      base.SetTDState_PT(p,t); check(); base.SetTDState_Prho(p,rho); check();
      base.SetTDState_rhoe(rho,e); check(); base.SetTDState_hs(h,s); check();
      base.SetTDState_Ps(p,s); check(); base.SetEnergy_Prho(p,rho); check();
      base.ComputeDerivativeNRBC_Prho(p,rho); check();
      const std::array<double,3> velocity{{700.,-100.,35.}},normal{{.3,-.4,.5}};
      const double velocity_squared=700.*700.+100.*100.+35.*35.;
      const double total_energy=e+.5*velocity_squared+k;
      const double total_enthalpy=total_energy+p/rho;
      const double kappa=fluid.GetdPde_rho()/rho;
      const double chi=rah_air::conservative_chi(fluid,k);
      const double sound_error=relative(chi+kappa*(total_enthalpy-.5*velocity_squared),fluid.GetSoundSpeed2());
      double storage[5][5]{}; double* matrix[5];
      for (unsigned i=0;i<5;++i) matrix[i]=storage[i];
      NativeGeneralJacobian native;
      native.GetInviscidProjJac(velocity.data(),&total_enthalpy,&chi,&kappa,normal.data(),1.,matrix);
      Vector u{{rho,rho*700.,-rho*100.,rho*35.,rho*total_energy}};
      const double speed=std::sqrt(velocity_squared)+fluid.GetSoundSpeed();
      const Vector state_scale{{rho,rho*speed,rho*speed,rho*speed,rho*std::abs(total_enthalpy)}};
      const double mass_scale=rho*speed*std::sqrt(.5);
      const Vector flux_scale{{mass_scale,mass_scale*speed,mass_scale*speed,mass_scale*speed,mass_scale*std::abs(total_enthalpy)}};
      double jacobian_error=0,temperature_error=0;
      for (unsigned j=0;j<5;++j) {
        const double step=1e-5*state_scale[j]; Vector minus=u,plus=u;
        minus[j]-=step; plus[j]+=step;
        double tminus=0,tplus=0;
        const auto fminus=flux(minus,k,normal,&tminus),fplus=flux(plus,k,normal,&tplus);
        for (unsigned i=0;i<5;++i)
          jacobian_error=std::max(jacobian_error,std::abs((fplus[i]-fminus[i])/(2*step)-matrix[i][j])*state_scale[j]/flux_scale[i]);
        double expected=0;
        if (j==0) expected=rah_air::temperature_density_derivative(fluid,total_energy,velocity_squared);
        else if (j==4) expected=rah_air::temperature_energy_derivative(fluid);
        else expected=-fluid.GetdTde_rho()*velocity[j-1]/rho;
        temperature_error=std::max(temperature_error,relative((tplus-tminus)/(2*step),expected));
      }
      const double step=t*1e-5;
      fluid.SetTDState_rhoT(rho,t-step); const double low=fluid.GetThermalConductivity();
      fluid.SetTDState_rhoT(rho,t+step); const double high=fluid.GetThermalConductivity();
      fluid.SetTDState_rhoT(rho,t); const double lambda=fluid.GetThermalConductivity();
      const double derivative_error=relative((high-low)/(2*step),fluid.GetdktdT_rho());
      NativeGeneralHeatJacobian heat;
      const double primitive[7]={t,velocity[0],velocity[1],velocity[2],p,rho,total_enthalpy};
      const double secondary[2]={fluid.GetdTdrho_e(),fluid.GetdTde_rho()};
      heat.SetHeatFluxJacobian(primitive,secondary,1e-4,lambda,fluid.GetCp(),.01);
      double heat_error=0;
      for (unsigned j=0;j<5;++j) {
        double expected=0;
        if (j==0) expected=rah_air::temperature_density_derivative(fluid,total_energy,velocity_squared);
        else if (j==4) expected=rah_air::temperature_energy_derivative(fluid);
        else expected=-fluid.GetdTde_rho()*velocity[j-1]/rho;
        expected*=(lambda+fluid.GetCp()*1e-4/.9)/.01;
        heat_error=std::max(heat_error,relative(heat.heat_flux_jac_i[j],expected));
      }
      fluid.SetEddyViscosity(.1);
      if (fluid.GetThermalConductivity()!=lambda) throw std::runtime_error("turbulent conductivity double counted");
      unsigned rejected=0;
      try {fluid.SetTDState_PT(p,2001.);} catch(const std::domain_error&) {++rejected;}
      try {fluid.SetTDState_rhoe(-rho,e);} catch(const std::domain_error&) {++rejected;}
      check();
      double thread_error=0;
      #pragma omp parallel for num_threads(2) reduction(max:thread_error)
      for (int i=0;i<20;++i) {
        CSU2IdealAir independent; independent.SetTDState_PT(p,t+i);
        thread_error=std::max(thread_error,relative(independent.GetTemperature(),t+i));
      }
      out<<t<<','<<rho<<','<<k<<','<<state_error<<','<<sound_error<<','<<jacobian_error<<','
         <<temperature_error<<','<<derivative_error<<','<<heat_error<<','<<rejected<<','<<thread_error<<'\n';
      if (state_error>2e-10 || sound_error>2e-11 || jacobian_error>3e-6 || temperature_error>3e-6 ||
          derivative_error>3e-6 || heat_error>2e-11 || rejected!=2 || thread_error>2e-11)
        throw std::runtime_error("adapter numerical gate failed");
      ++cases;
    }
    if (!out) throw std::runtime_error("failed CSV output");
    std::cout<<"adapter cases="<<cases<<" passed\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n'; return 1;
  }
}
