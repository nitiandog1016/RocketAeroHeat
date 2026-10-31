#pragma once
// 有界记录真实守恒状态/线性增量；调用者以SU2临界区串行调用，不改变诊断模式的更新。
#include "RAHFrozenEnergyRecovery.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace rah {
class FrozenEnergyLog {
  std::ofstream rows, counts;
  unsigned long written=0;
 public:
  FrozenEnergyLog() : rows("rah_frozen_directions.csv"), counts("rah_frozen_counts.csv") {
    if (!rows || !counts) throw std::runtime_error("Cannot create frozen-energy diagnostics");
    rows << "iteration,localNode,ndim,nativeAlpha,allowable,eligible,appliedAlpha,internalBefore,internalTrial";
    for (int j=0;j<5;++j) rows << ",q" << j;
    for (int j=0;j<5;++j) rows << ",delta" << j;
    rows << '\n' << std::setprecision(17);
    counts << "iteration,checked,frozen,eligible,applied\n";
  }
  void record(unsigned long iter,unsigned long node,unsigned ndim,double original,double allowable,
              const std::array<double,5>& q,const std::array<double,5>& delta,
              const FrozenEnergyDecision& d,double actual) {
    if(written>=4096) return;
    rows << iter << ',' << node << ',' << ndim << ',' << original << ',' << allowable << ','
         << d.eligible << ',' << actual << ',' << d.internalBefore << ',' << d.internalTrial;
    for(double v:q) rows << ',' << v;
    for(double v:delta) rows << ',' << v;
    rows << '\n';++written;
  }
  void count(unsigned long iter,unsigned long checked,unsigned long frozen,unsigned long eligible,unsigned long applied) {
    // 每线程每迭代一行，独立审核时按迭代聚合；详细方向最多4096行。
    counts << iter << ',' << checked << ',' << frozen << ',' << eligible << ',' << applied << '\n';
    counts.flush();rows.flush();
  }
};
}
