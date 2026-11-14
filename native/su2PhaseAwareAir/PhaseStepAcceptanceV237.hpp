// 模块职责：在成熟局部松弛上限内作具名相态可行性回溯；不计算通量、Jacobian或裁状态。
#pragma once
#include "fluid/PhaseAwareAir.hpp"
#include <array>
#include <cstring>
#include <cstdint>
namespace rah_phase_step_v237 {
using State=std::array<double,5>;
struct Result {double alpha;unsigned reductions;State accepted;};
using Observer=void(*)(unsigned long,unsigned long,const State&,const State&,double,const Result&);
inline Observer observer=nullptr;
inline uint64_t checked=0,limited=0;
inline bool enabled=true;

inline bool domain_rejection(const std::string& s){
 return s=="RAH_PHASE_ENERGY_ENTROPY_DOMAIN"||s=="RAH_PHASE_PT_DOMAIN"||s=="RAH_PHASE_DENSITY"||s=="RAH_PHASE_T_DOMAIN";
}
inline bool admissible(const State& u,double energyOffset=300000.){
 for(double x:u)if(!std::isfinite(x))return false;
 if(!(u[0]>0))return false;
 double v2=0;for(int j=1;j<4;++j){const double v=u[j]/u[0];v2+=v*v;}
 const double e=u[4]/u[0]-.5*v2+(300000.-energyOffset);
 if(!std::isfinite(e))return false;
 try{rah_phase_air::state("rhoe",u[0],e);return true;}
 catch(const std::domain_error& error){if(domain_rejection(error.what()))return false;throw;}
}
inline State candidate(const State& old,const State& delta,double alpha){
 State out{};for(int j=0;j<5;++j)out[j]=old[j]+alpha*delta[j];return out;
}
inline Result select(const State& old,const State& delta,double upper,double energyOffset=300000.){
 if(!std::isfinite(upper)||upper<0||upper>1||!std::isfinite(energyOffset))throw std::invalid_argument("PHASE_STEP_UPPER_OR_REFERENCE");
 for(double v:delta)if(!std::isfinite(v))throw std::invalid_argument("PHASE_STEP_NONFINITE_INCREMENT");
 if(!admissible(old,energyOffset))throw std::domain_error("PHASE_STEP_ILLEGAL_OLD");
 double alpha=upper;bool nonzero=false;for(double v:delta)nonzero=nonzero||v!=0;
 if(nonzero&&alpha<1e-10)throw std::domain_error("PHASE_STEP_NO_ADMISSIBLE_PROGRESS");
 for(unsigned k=0;k<=34;++k){
   State next=candidate(old,delta,alpha);
   if(admissible(next,energyOffset)){
     if(nonzero&&std::memcmp(next.data(),old.data(),sizeof(State))==0)throw std::domain_error("PHASE_STEP_NO_ADMISSIBLE_PROGRESS");
     return {alpha,k,next};
   }
   alpha*=.5;if(alpha<1e-10)break;
 }
 throw std::domain_error("PHASE_STEP_NO_ADMISSIBLE_PROGRESS");
}
inline double native_limit(unsigned long update,unsigned long point,const double* u,const double* d,double upper){
 if(!enabled)return upper;
 State old{},delta{};for(int j=0;j<5;++j){old[j]=u[j];delta[j]=d[j];}
 const auto r=select(old,delta,upper);++checked;
 if(r.reductions){++limited;if(observer)observer(update,point,old,delta,upper,r);}
 return r.alpha;
}
}
