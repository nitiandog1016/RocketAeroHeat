// 双支管银球单域使用：流动端登记本次 FDS 边界通量，标量端按同一有向顶点消费。
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace SilverFar {
struct Vertex {
  unsigned long global=0;
  std::array<double,3> normal{};
  double mass=0, used=0, adv=0, inside=0, outside=0;
  unsigned char status=0;
};
struct Exchange {
  const void* geometry=nullptr;
  unsigned long step=0;
  unsigned short marker=0;
  std::vector<Vertex> vertex;
};
inline Exchange& data() { static Exchange exchange;return exchange; }
inline void begin(const void* geo,unsigned long step,unsigned short marker,unsigned long count) {
  if(count!=1346) std::abort();
  auto& x=data();x.geometry=geo;x.step=step;x.marker=marker;x.vertex.assign(count,{});
}
inline void publish(unsigned long v,unsigned long global,const double* normal,double mass) {
  auto& f=data().vertex.at(v);
  if(f.status!=0 || !std::isfinite(mass)) std::abort();
  f.global=global;std::copy_n(normal,3,f.normal.begin());f.mass=mass;f.status=1;
}
inline double consume(const void* geo,unsigned long step,unsigned short marker,unsigned long v,
                      unsigned long global,const double* normal) {
  auto& x=data();auto& f=x.vertex.at(v);
  if(geo!=x.geometry || x.step!=step || x.marker!=marker || f.status!=1 || f.global!=global ||
     !std::equal(f.normal.begin(),f.normal.end(),normal)) std::abort();
  f.status=2;return f.mass;
}
inline bool capture() {
  static const bool enabled=[](){auto* f=std::fopen("capture.flag","rb");if(!f)return false;std::fclose(f);return true;}();
  return enabled;
}
inline void finish() {
  auto& x=data();double err=0;unsigned long inward=0,outward=0;
  FILE* detail=nullptr;
  if(capture() && x.step==0) {detail=std::fopen("flow_far.csv","wb");if(!detail)std::abort();}
  for(const auto& f:x.vertex) {
    if(f.status!=2 || f.used!=f.mass || !std::isfinite(f.adv) || f.outside!=0) std::abort();
    const double expected=f.mass*(f.mass>=0?f.inside:f.outside);
    err=std::max(err,std::abs(f.adv-expected));
    inward+=f.mass<0;outward+=f.mass>0;
    if(detail)std::fprintf(detail,"%lu,%.17g\n",f.global,f.mass);
  }
  if(detail)std::fclose(detail);
  if(err>1e-15)std::abort();
  auto* f=std::fopen("far_pairs.csv","ab");if(!f)std::abort();
  std::fprintf(f,"%lu,%zu,%lu,%lu,%.17g\n",x.step,x.vertex.size(),inward,outward,err);std::fclose(f);
}
}
