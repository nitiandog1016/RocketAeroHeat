// Short Test 65 诊断副本专用。日志只读取主方程与线性系统，不写入求解状态。
#pragma once
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <utility>

namespace ShortProbe {
struct Choice {
  long step, point;
  double ceiling;
  Choice() {
    auto* f = std::fopen("probe.in", "r");
    if (!f || std::fscanf(f, "%ld %ld %lf", &step, &point, &ceiling) != 3) std::abort();
    std::fclose(f);
  }
};
inline const Choice& request() {
  static const Choice value;
  return value;
}
template<class... Args>
void log(const char* name, const char* format, Args... args) {
#pragma omp critical(short_probe_file)
  {
    auto* f = std::fopen(name, "ab");
    if (!f) std::abort();
    std::fprintf(f, format, args...);
    std::fclose(f);
  }
}
}
