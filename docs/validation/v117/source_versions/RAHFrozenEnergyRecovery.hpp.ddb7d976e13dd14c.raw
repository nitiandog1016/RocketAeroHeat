#pragma once
// 仅为原生零松弛提供正内能恢复候选；不改变守恒增量、物性或空间离散。
// 原密度变化界限保留；正常非零松弛、能量下降或非法方向严格保持原行为。
#include <algorithm>
#include <cmath>
#include <array>
namespace rah {
struct FrozenEnergyDecision {
  double relaxation = 0.0, densityRelaxation = 0.0;
  double internalBefore = 0.0, internalTrial = 0.0;
  bool eligible = false;
};
inline double staticEnergyDensity(const std::array<double,5>& q, unsigned ndim) {
  double e=q[ndim+1];
  for (unsigned j=1;j<=ndim;++j) e-=0.5*q[j]*q[j]/q[0];
  return e;
}
inline FrozenEnergyDecision frozenEnergyRecovery(const std::array<double,5>& q,
    const std::array<double,5>& delta, unsigned ndim, double original, double allowable) {
  FrozenEnergyDecision d;d.relaxation=original;
  if (original!=0.0 || (ndim!=2 && ndim!=3) || !std::isfinite(allowable) || allowable<=0 || allowable>=1) return d;
  for (unsigned j=0;j<ndim+2;++j) if (!std::isfinite(q[j]) || !std::isfinite(delta[j])) return d;
  if (q[0]<=0) return d;
  d.internalBefore=staticEnergyDensity(q,ndim);
  if (!(d.internalBefore>0) || !std::isfinite(d.internalBefore)) return d;
  d.densityRelaxation=delta[0]==0 ? 1.0 : std::min(1.0,allowable*q[0]/std::abs(delta[0]));
  if (!std::isfinite(d.densityRelaxation) || d.densityRelaxation<1e-10) return d;
  std::array<double,5> trial=q;
  for (unsigned j=0;j<ndim+2;++j) trial[j]+=d.densityRelaxation*delta[j];
  for (unsigned j=0;j<ndim+2;++j) if (!std::isfinite(trial[j])) return d;
  if (trial[0]<=0) return d;
  d.internalTrial=staticEnergyDensity(trial,ndim);
  if (!std::isfinite(d.internalTrial) || d.internalTrial<=d.internalBefore) return d;
  d.eligible=true;d.relaxation=d.densityRelaxation;return d;
}
}
