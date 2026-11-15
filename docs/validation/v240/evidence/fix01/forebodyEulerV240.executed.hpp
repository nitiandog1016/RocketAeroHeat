// 模块职责：真实前体运行的出口阶段检查、原面压力收益停止与完整终点记录。
#pragma once
#include "RiemannOutletObserverV240.hpp"
#include <array>
#include <vector>
#include <cstring>
namespace region240 {
struct Face {unsigned long id,a,b,c;double area;};
std::vector<Face> faces;
std::vector<std::array<double,5>> rms;
std::vector<double> previous_pressure;
double latest_change=std::numeric_limits<double>::quiet_NaN();
unsigned long checked=0,stage_update=0;
double min_mach=1e100,max_mach=-1e100;
bool benefit_stop=false;

inline void boundary(unsigned long update,unsigned long point,double un,double a,
                     const double* ui,const double* ub,const double* normal){
  if(stage_update!=update){stage_update=update;checked=0;min_mach=1e100;max_mach=-1e100;}
  const double mn=un/a;min_mach=std::min(min_mach,mn);max_mach=std::max(max_mach,mn);
  if(!(std::isfinite(mn)&&a>0&&mn>1)||std::memcmp(ui,ub,5*sizeof(double))!=0){
    std::ofstream out(observation::folder+"/outlet_failure.csv");out<<std::setprecision(17)
      <<"updates,point,un,a,mn,Nx,Ny,Nz,rho,rhou,rhov,rhow,rhoE,b0,b1,b2,b3,b4\n"
      <<update<<','<<point<<','<<un<<','<<a<<','<<mn;
    for(int d=0;d<3;++d)out<<','<<normal[d];
    for(auto p:{ui,ub})for(int j=0;j<5;++j)out<<','<<p[j];out<<'\n';out.flush();
    throw std::runtime_error("REGION_OUTLET_INCOMING_OR_STATE_CHANGED_POINT_"+std::to_string(point));
  }
  ++checked;
}

inline void initialize(CGeometry* g,CSolver* f,CConfig* cfg){
  std::ifstream input(observation::folder+"/original_faces.csv");std::string line;std::getline(input,line);
  while(std::getline(input,line)){
    for(char& c:line)if(c==',')c=' ';std::istringstream row(line);Face q;
    if(!(row>>q.id>>q.a>>q.b>>q.c>>q.area)||q.a>=g->GetnPoint()||q.b>=g->GetnPoint()||q.c>=g->GetnPoint())
      throw std::runtime_error("REGION_FACE_MAP");faces.push_back(q);
  }
  if(faces.size()!=4434||cfg->GetRestart())throw std::runtime_error("REGION_INITIAL_IDENTITY");
  rah_outlet_v240::hook=boundary;
}

inline void record(CGeometry* g,CSolver* f,CConfig* cfg){
  const auto step=cfg->GetInnerIter()+1;
  if(stage_update!=step||checked!=807)throw std::runtime_error("REGION_NATIVE_OUTLET_COUNT");
  std::ofstream log(observation::folder+"/outlet_history.csv",step==1?std::ios::out:std::ios::app);
  if(step==1)log<<"updates,nativeBoundaryNodes,minNormalMach,maxNormalMach,bitEqualComponents,lastPreprocessingRejectedPoints\n";
  log<<std::setprecision(17)<<step<<','<<checked<<','<<min_mach<<','<<max_mach<<','<<5*checked<<','<<cfg->GetNonphysical_Points()<<'\n';
  std::array<double,5> r{};
  for(unsigned long i=0;i<g->GetnPointDomain();++i){
    double value[5];for(int j=0;j<5;++j)value[j]=f->LinSysRes(i,j);value[4]-=300000*value[0];
    for(int j=0;j<5;++j)r[j]+=value[j]*value[j];
  }
  for(double& x:r)x=std::sqrt(x/g->GetnPointDomain());rms.push_back(r);
  if(step%50==0){
    std::vector<double> p(g->GetnPoint());auto* n=f->GetNodes();
    for(unsigned long i=0;i<g->GetnPoint();++i){
      const double rho=n->GetSolution(i,0);double v2=0;
      for(int d=1;d<4;++d){const double v=n->GetSolution(i,d)/rho;v2+=v*v;}
      p[i]=rah_phase_air::state("rhoe",rho,n->GetSolution(i,4)/rho-.5*v2).p;
    }
    std::vector<double> surface;double difference=0,denominator=0;
    for(size_t i=0;i<faces.size();++i){auto q=faces[i];double value=(p[q.a]+p[q.b]+p[q.c])/3;surface.push_back(value);
      if(!previous_pressure.empty()){difference+=q.area*std::abs(value-previous_pressure[i]);denominator+=q.area*std::abs(previous_pressure[i]);}}
    latest_change=previous_pressure.empty()?std::numeric_limits<double>::quiet_NaN():difference/denominator;
    previous_pressure=surface;
    std::ofstream out(observation::folder+"/pressure_trend.csv",step==50?std::ios::out:std::ios::app);
    if(step==50)out<<"updates,areaPressureChangeFromPrevious50\n";
    out<<std::setprecision(17)<<step<<','<<latest_change<<'\n';
  }
  if(step>=200&&step%100==0){
    std::array<double,5> old{},now{},ratio{};double gm=0;bool worse=false;
    for(size_t i=step-100;i<step-50;++i)for(int j=0;j<5;++j)old[j]+=rms[i][j];
    for(size_t i=step-50;i<step;++i)for(int j=0;j<5;++j)now[j]+=rms[i][j];
    for(int j=0;j<5;++j){ratio[j]=now[j]/old[j];gm+=std::log(ratio[j])/5;worse=worse||ratio[j]>1.25;}
    gm=std::exp(gm);benefit_stop=worse||(gm>=.95&&latest_change>.01);
    std::ofstream out(observation::folder+"/benefit.csv",step==200?std::ios::out:std::ios::app);
    if(step==200)out<<"updates,ratio0,ratio1,ratio2,ratio3,ratio4,geometricMean,pressureChange,stop\n";
    out<<std::setprecision(17)<<step;for(double v:ratio)out<<','<<v;
    out<<','<<gm<<','<<latest_change<<','<<benefit_stop<<'\n';
  }
}
}
