// 模块职责：std::terminate退出前保留原异常及可能部分更新的原始守恒态；不恢复或继续推进。
#pragma once
#include <exception>
#include <cstdlib>
namespace observation {
CGeometry* failure_geometry=nullptr;CSolver* failure_flow=nullptr;CConfig* failure_config=nullptr;
const char* failure_phase="initialization";
[[noreturn]] inline void failed() noexcept {
  try {
    std::ofstream report(folder+"/termination.txt");
    report<<"phase="<<failure_phase<<"\nlast_observed_complete_update="<<last<<"\n";
    if(failure_config)report<<"native_inner_iteration="<<failure_config->GetInnerIter()<<"\n";
    if(auto error=std::current_exception()) {
      try{std::rethrow_exception(error);}
      catch(const std::exception& e){report<<"exception_what="<<e.what()<<"\n";}
      catch(...){report<<"exception_what=non-standard exception\n";}
    }else report<<"exception_what=no current exception available\n";
    report<<"state_role=possibly partially updated; not accepted; no continuation\n";report.flush();
    if(failure_geometry&&failure_flow){
      std::ofstream raw(folder+"/termination_raw_state.csv");raw<<std::setprecision(17)<<"PointID,Density,Momentum_x,Momentum_y,Momentum_z,Energy\n";
      for(unsigned long i=0;i<failure_geometry->GetnPoint();++i){raw<<failure_geometry->nodes->GetGlobalIndex(i);
        for(int j=0;j<5;++j)raw<<','<<failure_flow->GetNodes()->GetSolution(i,j);raw<<'\n';}
      raw.flush();
    }
  }catch(...){/* 原异常退出优先；记录失败不可触发第二次terminate或进入数值恢复。 */}
  std::_Exit(86);
}
}
