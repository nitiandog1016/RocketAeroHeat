// 模块职责：复用相态标准Air数据，提供SU2的具名能量坐标与严格联合域。
// 不求流场、不裁状态；反算使用已有Boost 1.87 TOMS748。
#pragma once
#include "PhaseAirTable.hpp"
#include <boost/math/tools/toms748_solve.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace rah_phase_air {
constexpr double R=287.0491267875601, shift=300000., Tmin=60., Tmax=2000., Pref=100000.;
inline void require(bool yes,const char* why) { if(!yes) throw std::domain_error(why); }
inline double ceiling(double t) {
  require(std::isfinite(t)&&t>=Tmin&&t<=Tmax,"RAH_PHASE_T_DOMAIN");
  const double ratio=t/100.;double p=std::min(1e6,1e4*((ratio*ratio)*ratio));
  if(t<=132.6312) {
    const double th=1-t/132.6312;
    const double value=-.1567266*std::sqrt(th)-5.539635*th+.7567212*std::pow(th,2.5)-3.514322*std::pow(th,4);
    p=std::min(p,.1*(3785020.*std::exp(132.6312/t*value)));
  }
  return p;
}
inline void accept(double p,double t) {
  require(std::isfinite(p)&&p>0&&p<=ceiling(t),"RAH_PHASE_PT_DOMAIN");
}
struct Thermal { double cp,h,s,mu,k,dcp,dmu,dk; };
inline Thermal thermal(double t) {
  require(std::isfinite(t)&&t>=Tmin&&t<=Tmax,"RAH_PHASE_T_DOMAIN");
  auto next=std::upper_bound(table.begin(),table.end(),t,[](double x,const Row& y){return x<y.t;});
  size_t i=next==table.end()?table.size()-2:static_cast<size_t>(next-table.begin()-1);
  const auto &a=table[i],&b=table[i+1];const double dt=t-a.t,width=b.t-a.t;
  const double dc=(b.cp-a.cp)/width,dm=(b.mu-a.mu)/width,dk=(b.k-a.k)/width;
  return {a.cp+dc*dt,a.h+a.cp*dt+.5*dc*dt*dt,
          a.s+(a.cp-dc*a.t)*std::log(t/a.t)+dc*dt,a.mu+dm*dt,a.k+dk*dt,dc,dm,dk};
}
// 模块二：量热数学反算与最终物理接受分开；试探只在已发表温区内。
enum class Quantity { E,H,S };
inline double inverse(double target,Quantity kind,double p=Pref) {
  require(std::isfinite(target)&&std::isfinite(p)&&p>0,"RAH_PHASE_INVERSE_INPUT");
  auto value=[=](double t) {
    auto a=thermal(t);
    if(kind==Quantity::E)return a.h-R*t+shift;
    if(kind==Quantity::H)return a.h+shift;
    return a.s-R*(std::log(p)-std::log(Pref));
  };
  const double low=value(Tmin),high=value(Tmax);
  require(target>=low&&target<=high,"RAH_PHASE_ENERGY_ENTROPY_DOMAIN");
  if(target==low)return Tmin;if(target==high)return Tmax;
  boost::uintmax_t count=80;
  auto bracket=boost::math::tools::toms748_solve([&](double t){return value(t)-target;},
      Tmin,Tmax,low-target,high-target,boost::math::tools::eps_tolerance<double>(48),count);
  require(boost::math::tools::eps_tolerance<double>(48)(bracket.first,bracket.second),"RAH_PHASE_INVERSE_CONVERGENCE");
  return bracket.first+.5*(bracket.second-bracket.first);
}
struct State {
  double p,t,rho,e,h,s,cp,cv,a2,dpdrho,dpde,dtdrho,dtde,dhdrho,dhdp,dsdrho,dsdp,mu,k,dmu,dk;
};
inline State PT(double p,double t) {
  accept(p,t);const auto a=thermal(t);const double rho=p/(R*t),cv=a.cp-R;
  State q{p,t,rho,a.h-R*t+shift,a.h+shift,a.s-R*(std::log(p)-std::log(Pref)),a.cp,cv,
          a.cp/cv*R*t,R*t,rho*R/cv,0,1/cv,-a.cp*t/rho,a.cp/(rho*R),-a.cp/rho,cv/p,a.mu,a.k,a.dmu,a.dk};
  // 全部字段先构造并核验，调用方再原子提交，不能留下半更新状态。
  const double fields[]={q.p,q.t,q.rho,q.e,q.h,q.s,q.cp,q.cv,q.a2,q.dpdrho,q.dpde,q.dtdrho,q.dtde,
    q.dhdrho,q.dhdp,q.dsdrho,q.dsdp,q.mu,q.k,q.dmu,q.dk};
  for(double x:fields)require(std::isfinite(x),"RAH_PHASE_NUMERICAL_RANGE");
  require(rho>0&&cv>0&&a.mu>0&&a.k>0,"RAH_PHASE_POSITIVITY");return q;
}
inline State rhoT(double rho,double t) {
  require(std::isfinite(rho)&&rho>0,"RAH_PHASE_DENSITY");return PT(rho*R*t,t);
}
inline State state(const std::string& kind,double x,double y) {
  require(std::isfinite(x)&&std::isfinite(y),"RAH_PHASE_NONFINITE");
  if(kind=="PT")return PT(x,y);
  if(kind=="rhoT")return rhoT(x,y);
  if(kind=="Prho") {require(x>0&&y>0,"RAH_PHASE_PRHO");return PT(x,x/(y*R));}
  if(kind=="rhoe")return rhoT(x,inverse(y,Quantity::E));
  if(kind=="Ps")return PT(x,inverse(y,Quantity::S,x));
  if(kind=="hs") {
    const double t=inverse(x,Quantity::H);
    const double p=Pref*std::exp((thermal(t).s-y)/R);return PT(p,t);
  }
  throw std::invalid_argument("RAH_PHASE_STATE_PAIR");
}
}
