// 与真实驱动共用 Budget 与 factor。在有余量的有界线性场及受限转移上验收。
#include "capacity.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
  using namespace Short65;
  const double k = 1./611.0086208842873, upper = .4683950617283952;
  for (const int n : {3, 9, 27}) {
    double total = 0, smallest = 1;
    // 剪切坐标 x=ξ+0.5η，θ=0.15+0.1η，ξ 截面非零通量完全来自梯度修正。
    // 在整个规则对偶网汇集两个方向的横向修正，再由两端预算决定共享系数。
    struct Link { int i,j,x; double flux; };
    const int side=n+1;
    const auto id=[=](int x,int y,int z){return (x*side+y)*side+z;};
    const auto w=[=](int p){return (p==0 || p==n) ? .5 : 1.;};
    std::vector<Budget> budgets(side*side*side);
    std::vector<Link> links;
    for(int x=0;x<=n;++x) for(int y=0;y<=n;++y) for(int z=0;z<=n;++z) {
      const int i=id(x,y,z);
      if(x<n) links.push_back({i,id(x+1,y,z),x,.05*k*w(y)*w(z)/(n*n)});
      // y 边两点投影给出完整通量的 0.8，其余 0.2 为斜网修正。
      if(y<n) links.push_back({i,id(x,y+1,z),-1,-.02*k*w(x)*w(z)/(n*n)});
    }
    for(const auto& e:links) {budgets[e.i].add(e.flux);budgets[e.j].add(-e.flux);}
    for(int x=0;x<=n;++x) for(int y=0;y<=n;++y) for(int z=0;z<=n;++z) {
      const double theta=.15+.1*y/n, vol=w(x)*w(y)*w(z)/(n*n*n);
      budgets[id(x,y,z)].close(theta,theta,vol/(.1/n),upper,false);
    }
    for(const auto& e:links) {
      const double a=factor(e.flux,budgets[e.i],budgets[e.j]);
      smallest=std::min(a,smallest);
      if(e.x==n/2) total+=a*e.flux;
    }
    std::printf("%d,%.17g,%.17g,%.17g\n",n,total,.05*k,smallest);
    if(std::abs(total-.05*k)>1e-12*k || smallest!=1.) return 2;
  }
  Budget i,j,kNode;
  i.add(.4);i.add(.6);j.add(-.4);kNode.add(-.6);
  i.close(.1,.05,2.,1.,false);j.close(.99,.98,1.,1.,false);kNode.close(.2,.3,1.,1.,false);
  const double a=.4*factor(.4,i,j), b=.6*factor(.6,i,kNode);
  std::printf("exchange,%.17g,%.17g,%.17g\n",a,b,.05-(a+b)/2);
  return std::abs(a-.01)<1e-14 && std::abs(b-.06)<1e-14 && .05-(a+b)/2>=0 ? 0 : 3;
}
