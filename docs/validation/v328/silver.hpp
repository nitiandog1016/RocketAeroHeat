// 单进程、四线程、定温双支管银球。只在专用标量对象中引用。
#pragma once
#include "capacity.hpp"
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

namespace Short65 {
struct Edge {
  unsigned long i=0,j=0;
  double mass=0, high=0, correction=0, extraDiagonal=0, alpha=1;
};
struct Node {
  unsigned long global=0;
  double old=0, low=0, final=0, volume=0, dt=0, native=0, raw=0, correction=0, limited=0;
  double farMass=0, farFlux=0, divergence=0, weight=0;
  bool wall=false;
  Budget budget;
};
static std::vector<Edge> edges;
static std::vector<Node> nodesData;
static double upper=0, maxRowError=0;
static unsigned long wallCount=0;
inline FILE* file(const char* path,const char* mode) {
  auto* f=std::fopen(path,mode);
  if(!f) std::abort();
  return f;
}
inline void begin(unsigned long count,unsigned long edgeCount) {
  nodesData.assign(count,Node{});
  edges.resize(edgeCount);
  upper=0;wallCount=0;maxRowError=0;
}
inline void reduceEdges() {
  for(const auto& e:edges) {
    auto& a=nodesData[e.i];auto& b=nodesData[e.j];
    a.raw+=e.high;b.raw-=e.high;
    a.divergence-=a.old*e.mass;b.divergence+=b.old*e.mass;
    a.correction+=e.correction;b.correction-=e.correction;
    a.budget.add(e.correction);b.budget.add(-e.correction);
  }
  for(auto& a:nodesData) {
    a.divergence-=a.old*a.farMass;
    a.raw+=a.divergence+a.farFlux;
    if(!a.wall) maxRowError=std::max(maxRowError,std::abs(a.raw-a.native));
  }
}
inline void limit() {
  for(auto& n:nodesData) n.budget.close(n.old,n.low,n.weight,upper,n.wall);
  for(auto& e:edges) {
    auto& a=nodesData[e.i];auto& b=nodesData[e.j];
    e.alpha=factor(e.correction,a.budget,b.budget);
    a.limited+=e.alpha*e.correction;b.limited-=e.alpha*e.correction;
  }
}
inline void snapshot(unsigned long worst) {
  auto* f=file("failure_nodes.csv","wb");
  std::fprintf(f,"id,old,low,final,wall,volume,dt,native,raw,correction,limited,farMass,farFlux,divergence,loss,gain,give,take\n");
  for(const auto& n:nodesData) std::fprintf(f,"%lu,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
      n.global,n.old,n.low,n.final,int(n.wall),n.volume,n.dt,n.native,n.raw,n.correction,n.limited,
      n.farMass,n.farFlux,n.divergence,n.budget.loss,n.budget.gain,n.budget.give,n.budget.take);
  std::fclose(f);
  f=file("failure_edges.csv","wb");
  std::fprintf(f,"edge,i,j,mass,high,correction,alpha,extraDiagonal\n");
  for(size_t k=0;k<edges.size();++k) {
    const auto& e=edges[k];if(e.i!=worst && e.j!=worst) continue;
    std::fprintf(f,"%zu,%lu,%lu,%.17g,%.17g,%.17g,%.17g,%.17g\n",k,nodesData[e.i].global,nodesData[e.j].global,
        e.mass,e.high,e.correction,e.alpha,e.extraDiagonal);
  }
  std::fclose(f);
}
inline void audit(unsigned long step,bool predicted) {
  double lo=1e30,hi=-1e30,maxWall=0,freeSum=0,wallSum=0,far=0,source=0,mag=0,highRms=0,schemeRms=0;
  unsigned long below=0,above=0,cut=0,minNode=0,maxNode=0;
  for(size_t i=0;i<nodesData.size();++i) {
    const auto& n=nodesData[i];
    const double val=predicted?n.low:n.final;
    if(!std::isfinite(val) || !std::isfinite(n.raw) || !(n.weight>0)) std::abort();
    if(val<lo) {lo=val;minNode=i;} if(val>hi) {hi=val;maxNode=i;}
    below+=val < -1e-12;above+=val > upper+1e-12;
    const double correction=predicted?0:n.limited;
    const double raw=n.raw-n.correction+correction;
    if(n.wall) {wallSum+=raw;maxWall=std::max(maxWall,std::abs(val-n.old));}
    else {const double assembled=n.native-n.correction+correction;freeSum+=assembled;
      highRms+=n.native*n.native;schemeRms+=assembled*assembled;}
    far+=n.farFlux;source+=n.divergence;
    mag+=std::abs(n.raw)+std::abs(n.correction)+std::abs(correction)+std::abs(n.farFlux)+std::abs(n.divergence);
  }
  const double balance=freeSum+wallSum-far-source;
  if(!predicted) for(const auto& e:edges) if(e.alpha<1 && e.correction!=0) ++cut;
  auto* f=file("stages.csv","ab");
  if(step==0 && predicted) std::fprintf(f,"iteration,stage,min,max,minID,maxID,below,above,wallDrift,limitedEdges,freeSpatial,wallSpatial,far,divergence,balance,scale,rowError,rawRMS,schemeRMS\n");
  std::fprintf(f,"%lu,%s,%.17g,%.17g,%lu,%lu,%lu,%lu,%.17g,%lu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
       step,predicted?"low":"actual",lo,hi,nodesData[minNode].global,nodesData[maxNode].global,below,above,maxWall,cut,
       freeSum,wallSum,far,source,balance,mag,maxRowError,std::sqrt(highRms/nodesData.size()),std::sqrt(schemeRms/nodesData.size()));
  std::fclose(f);
  if(below || above || maxWall>2e-14 || maxRowError>2e-13 || std::abs(balance)>2e-12*std::max(1.,mag)) {
    snapshot(above && hi-upper > -lo?maxNode:minNode);
    std::exit(42);
  }
}
}
