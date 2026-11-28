// 固定低速算例的只读观察。每条有效 SIMD 边使用独立槽，阶段写出前同步。
#pragma once
#include <array>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace TailBudget {
struct Edge {
  unsigned long left=0, right=0;
  double central=0, roe=0, total=0;
  bool stored=false;
};
inline bool capturing=false;
inline std::vector<Edge> edges;
inline thread_local std::array<double,16> central{};

inline bool selected(unsigned long id) {
  const auto col=id%109, row=id/109;
  return col>=11 && col<=14 && row>0 && row<88;
}

template<class C> bool watch(const C& c) {
  static const bool enabled=std::getenv("RAH_M030_TAIL_BUDGET")!=nullptr;
  const auto i=c.GetInnerIter();
  return enabled && c.GetKind_TimeIntScheme()==EULER_IMPLICIT && c.GetTimeIter()==1 &&
         (i==0 || i==1 || i==249);
}

inline FILE* output(unsigned long iter,const char* name) {
  auto path=std::to_string(iter)+"_"+name+".csv";
  auto* f=std::fopen(path.c_str(),"wb");
  if(!f) std::abort();
  return f;
}

template<class G,class S,class C>
void sample(const char* stage,G* g,S* flow,C* c) {
  if(!watch(*c)) return;
  // 所有线程在同一个原生阶段到达后，才允许单线程读取整个向量。
  #pragma omp barrier
  #pragma omp single
  {
    const std::string name(stage);
    if(g->GetnPointDomain()!=9701 || flow->GetnVar()!=4 || !c->GetAxisymmetric()) std::abort();
    if(name=="start") {
      edges.assign(g->GetnEdge(),Edge{});
      capturing=true;
    }
    auto* f=output(c->GetInnerIter(),stage);
    std::fprintf(f,"point,x,r,energyResidual,energy,previousEnergy,volume,dtND\n");
    for(unsigned long i=0;i<g->GetnPointDomain();++i) {
      const auto id=g->nodes->GetGlobalIndex(i);
      if(!selected(id)) continue;
      auto* xy=g->nodes->GetCoord(i);
      std::fprintf(f,"%lu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",id,
        double(xy[0]),double(xy[1]),double(flow->LinSysRes(i,3)),double(flow->GetNodes()->GetSolution(i,3)),
        double(flow->GetNodes()->GetSolution_time_n(i)[3]),double(g->nodes->GetVolume(i)),double(c->GetDelta_UnstTimeND()));
    }
    std::fclose(f);
    if(name=="transport") {
      auto* e=output(c->GetInnerIter(),"edges");
      std::fprintf(e,"edge,left,right,central,roe,total\n");
      for(unsigned long j=0;j<edges.size();++j) {
        const auto& v=edges[j];
        if(v.stored) std::fprintf(e,"%lu,%lu,%lu,%.17g,%.17g,%.17g\n",j,v.left,v.right,v.central,v.roe,v.total);
      }
      std::fclose(e);
      capturing=false;
    }
  }
}
} // namespace TailBudget
