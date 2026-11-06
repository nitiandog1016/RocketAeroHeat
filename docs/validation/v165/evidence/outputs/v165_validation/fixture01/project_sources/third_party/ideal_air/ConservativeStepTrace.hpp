#pragma once
// 模块用途。记录有限守恒限步及真实调用上下文，不修改流场或几何。
#include <array>
#include <cstdio>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace rah_step {
using State = std::array<double,5>;
struct Context {
  bool known=false, domain=false, output=false;
  unsigned long global=0, inner=0, time=0, mesh=0, rk=0;
  std::array<double,3> coordinates{};
  double volume=0, alpha=0;
  State linear{}, residual{};
};
// 每个OpenMP线程只读取其本次预处理的几何与求解器；退出作用域后恢复。
inline thread_local std::function<Context(unsigned long)> provider;
struct ContextGuard {
  std::function<Context(unsigned long)> previous;
  explicit ContextGuard(std::function<Context(unsigned long)> current):previous(std::move(provider)) {
    provider=std::move(current);
  }
  ~ContextGuard() { provider=std::move(previous); }
};

inline void record(unsigned long local, unsigned long dimension, const State& trial,
                   const State& old, const State& accepted, double trial_e, double old_e,
                   double accepted_e, double k, double theta, unsigned attempts,
                   bool old_valid, const char* reason) {
  const Context context=provider ? provider(local) : Context{};
  struct Trace {
    std::mutex mutex;
    std::FILE* file=nullptr;
    unsigned long count=0;
    ~Trace() { if(file) std::fclose(file); }
  };
  static Trace trace;
  const std::lock_guard<std::mutex> lock(trace.mutex);
  if(!trace.file) {
    trace.file=std::fopen("conservative_steps.csv","wb");
    if(!trace.file) throw std::runtime_error("cannot create conservative_steps.csv");
    std::fprintf(trace.file,"event,local_point,dimension,old_valid,theta,attempts,k,context_known,global_point,domain,inner_iter,time_iter,mesh,rk,output,x,y,z,volume,alpha");
    for(const char* prefix:{"trial","old","accepted","linear","stored_residual"})
      for(const char* name:{"rho","rhou","rhov","rhow","rhoE"})
        std::fprintf(trace.file,",%s_%s",prefix,name);
    std::fprintf(trace.file,",trial_e,old_e,accepted_e,reason\n");
  }
  if(trace.count>=4096) throw std::runtime_error("conservative step trace limit reached after 4096 events");
  ++trace.count;
  std::fprintf(trace.file,"%lu,%lu,%lu,%u,%.17g,%u,%.17g,%u,%lu,%u,%lu,%lu,%lu,%lu,%u,%.17g,%.17g,%.17g,%.17g,%.17g",
    trace.count,local,dimension,old_valid?1u:0u,theta,attempts,k,context.known?1u:0u,
    context.global,context.domain?1u:0u,context.inner,context.time,context.mesh,context.rk,
    context.output?1u:0u,context.coordinates[0],context.coordinates[1],context.coordinates[2],context.volume,context.alpha);
  for(const State* state:{&trial,&old,&accepted,&context.linear,&context.residual})
    for(double value:*state) std::fprintf(trace.file,",%.17g",value);
  std::fprintf(trace.file,",%.17g,%.17g,%.17g,\"%s\"\n",trial_e,old_e,accepted_e,reason);
  if(std::ferror(trace.file) || std::fflush(trace.file)!=0)
    throw std::runtime_error("cannot flush conservative step record");
}
} // namespace rah_step
