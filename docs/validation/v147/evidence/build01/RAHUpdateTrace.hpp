// RocketAeroHeat：只读保存一次隐式更新的真实输入/增量/松弛，不修改数值状态。
#pragma once
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>

template<class Nodes, class Vector, class Density>
inline void RAHWriteUpdateTrace(const char* name, bool after, unsigned long count,
                               unsigned short variables, Nodes* nodes,
                               const Vector& increments, Density density) {
  const char* enabled = std::getenv("RAH_SST_UPDATE_TRACE");
  if (!enabled || std::string(enabled) != "1") return;
  // 调用点由SU2全局同步保护，只有主线程读写；重复调用拒绝覆盖证据。
  const std::string filename = std::string(name) + (after ? "_update_after.csv" : "_update_before.csv");
  if (std::ifstream(filename).good()) SU2_MPI::Error("Update trace already exists", CURRENT_FUNCTION);
  std::ofstream out(filename);
  if (!out) SU2_MPI::Error("Cannot create update trace", CURRENT_FUNCTION);
  out << std::setprecision(17) << std::scientific << "PointID";
  if (!after) out << ",DensityUsed,Relaxation";
  for (unsigned short v=0; v<variables; ++v) out << ",Solution_" << v;
  if (!after) {
    for (unsigned short v=0; v<variables; ++v) out << ",Old_" << v;
    for (unsigned short v=0; v<variables; ++v) out << ",Increment_" << v;
  }
  out << '\n';
  for (unsigned long i=0; i<count; ++i) {
    out << i;
    if (!after) out << ',' << SU2_TYPE::GetValue(density(i)) << ',' << SU2_TYPE::GetValue(nodes->GetUnderRelaxation(i));
    for (unsigned short v=0; v<variables; ++v) out << ',' << SU2_TYPE::GetValue(nodes->GetSolution(i,v));
    if (!after) {
      for (unsigned short v=0; v<variables; ++v) out << ',' << SU2_TYPE::GetValue(nodes->GetSolution_Old(i,v));
      for (unsigned short v=0; v<variables; ++v) out << ',' << SU2_TYPE::GetValue(increments(i,v));
    }
    out << '\n';
  }
  out.close();
  if (!out) SU2_MPI::Error("Incomplete update trace", CURRENT_FUNCTION);
}
