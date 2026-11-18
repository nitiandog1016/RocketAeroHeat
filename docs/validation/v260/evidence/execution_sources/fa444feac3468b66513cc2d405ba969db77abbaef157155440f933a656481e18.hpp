// 模块职责：同一相态Air表的稳定E/H温度反算；不改物理模型/温域/相态拒绝。
#pragma once
#include "PhaseAwareAir.hpp"

namespace rah_phase_air_stable {
using rah_phase_air::Quantity;
struct InverseInfo { const char* method="not-called";size_t interval=0; };

// 模块一：严格复用已有前向量与节点目标，避免变更原端点准入语义。
inline double caloricValue(double t,Quantity kind){
  auto a=rah_phase_air::thermal(t);
  return kind==Quantity::E?a.h-rah_phase_air::R*t+rah_phase_air::shift:a.h+rah_phase_air::shift;
}
inline double inverse(double target,Quantity kind,double p=rah_phase_air::Pref,InverseInfo* info=nullptr){
  using namespace rah_phase_air;
  require(std::isfinite(target)&&std::isfinite(p)&&p>0,"RAH_PHASE_INVERSE_INPUT");
  if(kind==Quantity::S){if(info)info->method="legacy-entropy-TOMS748";return rah_phase_air::inverse(target,kind,p);}
  const double low=caloricValue(Tmin,kind),high=caloricValue(Tmax,kind);
  require(target>=low&&target<=high,"RAH_PHASE_ENERGY_ENTROPY_DOMAIN");
  size_t left=0,right=table.size();
  while(left<right){const auto mid=left+(right-left)/2;const double value=caloricValue(table[mid].t,kind);
    if(value<target)left=mid+1;else right=mid;}
  require(left<table.size(),"RAH_PHASE_CALORIC_SEARCH");
  if(target==caloricValue(table[left].t,kind)){
    if(info){info->method="exact-forward-node";info->interval=left==table.size()-1?left-1:left;}return table[left].t;
  }
  require(left>0,"RAH_PHASE_CALORIC_SEARCH");const size_t i=left-1;const auto& a=table[i];const auto& b=table[i+1];
  // 模块二：有理化二次根。扩展精度只是中间算术，输入表/接口仍为原double。
  const long double width=static_cast<long double>(b.t)-a.t;
  const long double slope=(static_cast<long double>(b.cp)-a.cp)/width;
  const long double capacity=static_cast<long double>(a.cp)-(kind==Quantity::E?R:0.);
  const long double anchor=static_cast<long double>(a.h)-(kind==Quantity::E?static_cast<long double>(R)*a.t:0.L)+shift;
  const long double delta=static_cast<long double>(target)-anchor;
  const long double discriminant=capacity*capacity+2*slope*delta;
  require(capacity>0&&discriminant>0,"RAH_PHASE_CALORIC_MONOTONICITY");
  const long double dt=2*delta/(capacity+std::sqrt(discriminant));
  const double t=static_cast<double>(static_cast<long double>(a.t)+dt);
  if(info)info->interval=i;
  // 模块三：保存的锚点舍入可能形成极小间隙；显式原方法回退，不裁态或外推段。
  if(!(t>=a.t&&t<=b.t)||!std::isfinite(t)){
    if(info)info->method="legacy-knot-rounding-fallback";return rah_phase_air::inverse(target,kind,p);
  }
  if(info)info->method="stable-piecewise-quadratic";return t;
}
}
