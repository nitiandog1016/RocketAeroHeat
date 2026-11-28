// 在固定的 6.25 µs 终态读取 SU2 原生组装。按全局点号或边号保存各线程独占的槽。
#pragma once
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

namespace HeatLedger {
inline bool capture=false;
inline bool at_far=false;
inline std::map<const double*,unsigned long> primitive_ids;
inline std::vector<std::array<double,6>> edges;
inline std::vector<std::array<double,10>> far;
inline std::vector<unsigned char> visc_seen;

inline void far_window(bool active) {
  #pragma omp barrier
  #pragma omp single
  { at_far=active && capture; }
}

template<class T> void dump(const std::string& name,const std::vector<T>& rows) {
  auto* f=std::fopen((name+".bin").c_str(),"wb");
  if(!f || std::fwrite(rows.data(),sizeof(T),rows.size(),f)!=rows.size() || std::fclose(f)) std::abort();
}

template<class G,class S,class C> void begin(G* g,S* s,C* c) {
  if(!std::getenv("RAH_M030_STEP_LEDGER")) return;
  #pragma omp barrier
  #pragma omp single
  {
    const auto n=g->GetnPointDomain();
    if(capture || n!=9701 || n!=g->GetnPoint() || !c->GetAxisymmetric() || s->GetnVar()!=4 ||
       c->GetTimeIter()!=1 || c->GetInnerIter()!=0 || std::abs(c->GetMach()-.30)>1e-12) std::abort();
    // 求解器从步末 CSV 载入当前状态。这里仅还原后向欧拉方程的上一时间层。
    std::vector<std::array<double,4>> previous(n);
    auto* f=std::fopen("previous.bin","rb");
    if(!f || std::fread(previous.data(),sizeof(previous[0]),n,f)!=n || std::fgetc(f)!=EOF || std::fclose(f)) std::abort();
    std::vector<std::array<double,11>> state(n);
    for(unsigned long p=0;p<n;++p) {
      const auto id=g->nodes->GetGlobalIndex(p);
      if(id>=n) std::abort();
      if(!primitive_ids.emplace(s->GetNodes()->GetPrimitive(p),id).second) std::abort();
      auto& row=state[id];
      const auto* xy=g->nodes->GetCoord(p);
      row[0]=xy[0];row[1]=xy[1];row[2]=g->nodes->GetVolume(p);
      for(unsigned j=0;j<4;++j) {
        s->GetNodes()->Set_Solution_time_n(p,j,previous[id][j]);
        row[3+j]=s->GetNodes()->GetSolution(p,j);
        row[7+j]=s->GetNodes()->GetSolution_time_n(p)[j];
      }
    }
    dump("state",state);
    std::vector<std::array<double,7>> ref={{{c->GetDensity_Ref(),c->GetVelocity_Ref(),c->GetPressure_Ref(),
          c->GetEnergy_Ref(),c->GetLength_Ref(),c->GetTime_Ref(),c->GetDelta_UnstTimeND()}}};
    dump("scales",ref);
    edges.resize(g->GetnEdge(),{{-1,-1,0,0,0,0}});
    far.resize(n);
    visc_seen.resize(n);
    capture=true;
  }
}

template<class G,class S> void sample(const char* name,G* g,S* s) {
  if(!capture) return;
  #pragma omp barrier
  #pragma omp single
  {
    std::vector<std::array<double,4>> rows(g->GetnPointDomain());
    for(unsigned long p=0;p<g->GetnPointDomain();++p) {
      auto& row=rows[g->nodes->GetGlobalIndex(p)];
      for(unsigned j=0;j<4;++j) row[j]=s->LinSysRes(p,j);
    }
    dump(name,rows);
    if(std::string(name)=="rhs") {
      for(unsigned long p=0;p<far.size();++p) if(far[p][0]!=visc_seen[p]) std::abort();
      dump("edges",edges);
      dump("far",far);
      capture=false;
    }
  }
}
}
