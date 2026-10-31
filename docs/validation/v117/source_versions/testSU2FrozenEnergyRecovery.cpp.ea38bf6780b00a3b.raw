// 独立构造正/负恢复方向、非零原松弛、非法输入和随机线段；不依赖真实流场裁剪。
#include "../third_party/su2_thermal/RAHFrozenEnergyRecovery.hpp"
#include <cassert>
#include <iostream>
#include <random>
#include <limits>
int main() {
  std::array<double,5> q={{1,800,0,0,320000+1e-8}},up={{0,0,0,0,100}};
  const auto d=rah::frozenEnergyRecovery(q,up,3,0,.2);
  assert(d.eligible && d.relaxation==1 && d.internalTrial>d.internalBefore);
  assert(rah::frozenEnergyRecovery(q,up,3,.01,.2).relaxation==.01);
  up[4]=-100;assert(!rah::frozenEnergyRecovery(q,up,3,0,.2).eligible);
  up[4]=std::numeric_limits<double>::quiet_NaN();assert(!rah::frozenEnergyRecovery(q,up,3,0,.2).eligible);
  q[0]=-1;assert(!rah::frozenEnergyRecovery(q,up,3,0,.2).eligible);
  std::mt19937_64 rng(117);std::uniform_real_distribution<double> unit(-1,1);
  unsigned accepted=0;
  for (unsigned k=0;k<20000;++k) {
    const unsigned nd=k%2?2:3;std::array<double,5> a{},delta{};
    a[0]=std::exp(3*unit(rng));a[nd+1]=a[0]*std::exp(8+6*unit(rng));
    for(unsigned j=1;j<=nd;++j) {a[j]=a[0]*1000*unit(rng);a[nd+1]+=.5*a[j]*a[j]/a[0];}
    for(unsigned j=0;j<nd+2;++j) delta[j]=2*std::max(1.0,std::abs(a[j]))*unit(rng);
    auto r=rah::frozenEnergyRecovery(a,delta,nd,0,.2);
    if(!r.eligible) continue;++accepted;
    assert(r.relaxation>0 && r.relaxation<=1);
    for(unsigned n=0;n<=20;++n) {
      const double alpha=r.relaxation*n/20;auto b=a;
      for(unsigned j=0;j<nd+2;++j)b[j]+=alpha*delta[j];
      assert(b[0]>0 && std::abs(b[0]/a[0]-1)<=.2+1e-14);
      assert(rah::staticEnergyDensity(b,nd)>0);
    }
  }
  assert(accepted>1000);std::cout<<"eligible_random_directions="<<accepted<<"; endpoints and full segments positive\n";
}
