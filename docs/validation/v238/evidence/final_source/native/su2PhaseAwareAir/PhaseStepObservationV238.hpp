// 模块职责：为既有可行性守卫补具名失败观察；不改变选择算法或任何数值结果。
#pragma once
#include "PhaseStepAcceptanceV237.hpp"
namespace rah_phase_step_v238 {
using Failure=void(*)(unsigned long,unsigned long,const double*,const double*,double,const char*);
inline Failure failure=nullptr;
inline double native_limit(unsigned long update,unsigned long point,const double* u,const double* d,double upper){
 try{return rah_phase_step_v237::native_limit(update,point,u,d,upper);}
 catch(const std::exception& e){if(failure)failure(update,point,u,d,upper,e.what());throw;}
}
}
