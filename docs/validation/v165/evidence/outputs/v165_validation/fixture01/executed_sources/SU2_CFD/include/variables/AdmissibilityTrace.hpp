#pragma once
// 模块用途。对实际物性回退保存有限、线程安全的守恒态证据。
#include <cstdio>
#include <mutex>
#include <stdexcept>

namespace rah_admissibility {
inline void record(unsigned long local_point, unsigned long dimension, const double* trial,
                   const double* restored, double trial_energy, double restored_energy,
                   double kinetic, bool recovered, const char* reason) {
  struct Trace {
    std::mutex mutex;
    std::FILE* file = nullptr;
    unsigned long count = 0;
    ~Trace() { if (file) std::fclose(file); }
  };
  static Trace trace;
  const std::lock_guard<std::mutex> guard(trace.mutex);
  if (!trace.file) {
    trace.file = std::fopen("admissibility.csv", "wb");
    if (!trace.file) throw std::runtime_error("cannot create admissibility.csv");
    std::fprintf(trace.file, "event,local_point,dimension,recovered,k,trial_rho,trial_rhou,trial_rhov,trial_rhow,trial_rhoE,trial_e,old_rho,old_rhou,old_rhov,old_rhow,old_rhoE,old_e,reason\n");
  }
  if (trace.count >= 4096) throw std::runtime_error("admissibility trace limit reached after 4096 events");
  ++trace.count;
  const int written = std::fprintf(trace.file,
      "%lu,%lu,%lu,%u,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,\"%s\"\n",
      trace.count, local_point, dimension, recovered ? 1u : 0u, kinetic,
      trial[0], trial[1], trial[2], trial[3], trial[4], trial_energy,
      restored[0], restored[1], restored[2], restored[3], restored[4], restored_energy, reason);
  if (written < 0 || std::fflush(trace.file) != 0)
    throw std::runtime_error("cannot flush admissibility state record");
}
}  // namespace rah_admissibility
