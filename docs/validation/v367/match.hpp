// 将已保存的官方工作状态装入单次观测。原生物理上一步和网格保持既有输入身份。
#pragma once
#include <array>
#include <cstdio>
#include <cstdlib>
#include <vector>
namespace TailBudget {
template<class G,class S,class C> void bind_state(G* g,S* flow,C* c) {
  if(!std::getenv("RAH_M030_ASSEMBLY_STATE") || c->GetInnerIter()!=0) return;
  #pragma omp barrier
  #pragma omp single
  {
    const size_t n=g->GetnPointDomain();
    if(n!=9701 || c->GetTimeIter()!=1 || flow->GetnVar()!=4) std::abort();
    std::vector<std::array<double,4>> target(n);
    FILE* f=std::fopen("assembly_state.bin","rb");
    if(!f || std::fread(target.data(),sizeof(target[0]),n,f)!=n || std::fgetc(f)!=EOF) std::abort();
    std::fclose(f);
    std::vector<std::array<double,8>> verified(n);
    for(size_t p=0;p<n;++p) {
      const auto id=g->nodes->GetGlobalIndex(p);
      if(id>=n) std::abort();
      for(size_t k=0;k<4;++k) {
        flow->GetNodes()->SetSolution(p,k,target[id][k]);
        verified[id][k]=flow->GetNodes()->GetSolution(p,k);
        verified[id][k+4]=flow->GetNodes()->GetSolution_time_n(p)[k];
      }
    }
    f=std::fopen("bound_state.bin","wb");
    if(!f || std::fwrite(verified.data(),sizeof(verified[0]),n,f)!=n) std::abort();
    std::fclose(f);
  }
}
}
