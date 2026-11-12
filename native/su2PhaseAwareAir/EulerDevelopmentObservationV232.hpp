// 模块职责：只读核对已验证轨迹，记录最后原生预处理的拒绝点数，不参与更新。
#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
namespace observation {
std::vector<std::string> baseline_physical,baseline_native;
inline std::vector<std::string> complete_lines(const std::string& name) {
  std::ifstream in(name);std::vector<std::string> rows;std::string line;
  while(std::getline(in,line)){
    if(in.eof())break; // 未落盘的末尾半行不作为可核对数据。
    if(!line.empty()&&line.back()=='\r')line.pop_back();rows.push_back(line);
  }
  return rows;
}
inline void load_baseline() {
  baseline_physical=complete_lines(folder+"/baseline_physical.csv");
  baseline_native=complete_lines(folder+"/baseline_native.csv");
  if(baseline_physical.size()!=31||baseline_native.size()!=31)throw std::runtime_error("BASELINE_30_ROWS_REQUIRED");
}
inline void verify_native_prefix(unsigned long completed) {
  if(completed>30)return;
  const auto rows=complete_lines(folder+"/history.csv");
  for(size_t i=0;i<rows.size()&&i<31;++i)
    if(rows[i]!=baseline_native[i])throw std::runtime_error("NATIVE_PREFIX_CHANGED_ROW_"+std::to_string(i));
  if(completed==30){
    if(rows.size()<31)throw std::runtime_error("NATIVE_PREFIX_30_INCOMPLETE");
    std::ofstream out(folder+"/prefix_verified.json");out<<"{\"first30NativeAndPhysicalLinesEqual\":true}\n";
    if(!out)throw std::runtime_error("PREFIX_REPORT_WRITE");
  }
}
inline void record_rejections_and_check(CConfig* cfg) {
  const auto step=cfg->GetInnerIter()+1;
  if(step<=30){
    const auto rows=complete_lines(folder+"/physical_history.csv");
    if(rows.size()!=step+1||rows[0]!=baseline_physical[0]||rows[step]!=baseline_physical[step])
      throw std::runtime_error("PHYSICAL_PREFIX_CHANGED_STEP_"+std::to_string(step));
  }
  std::ofstream out(folder+"/native_rejections.csv",step==1?std::ios::out:std::ios::app);
  if(step==1)out<<"updates,last_preprocessing_nonphysical_points\n";
  out<<step<<','<<cfg->GetNonphysical_Points()<<'\n';
  if(!out)throw std::runtime_error("NATIVE_REJECTIONS_WRITE");
}
}
