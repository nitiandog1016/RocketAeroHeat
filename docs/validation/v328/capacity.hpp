// Short Test 65 专用。预算以原生标量通量为单位，同一内部边只乘一次系数。
#pragma once
#include <algorithm>
#include <cmath>

namespace Short65 {
struct Budget {
  double loss = 0, gain = 0, give = 1, take = 1;
  void add(double outward) {
    if (outward > 0) loss += outward;
    else gain -= outward;
  }
  void close(double old, double predictor, double mass, double upper, bool prescribed) {
    if (prescribed) return;
    const auto floorRoom = std::max(0., std::min(old, predictor));
    const auto ceilingRoom = std::max(0., upper - std::max(old, predictor));
    if (loss > 0) give = std::min(1., mass * floorRoom / loss);
    if (gain > 0) take = std::min(1., mass * ceilingRoom / gain);
  }
};
inline double factor(double flux, const Budget& i, const Budget& j) {
  return flux >= 0 ? std::min(i.give, j.take) : std::min(i.take, j.give);
}
}
