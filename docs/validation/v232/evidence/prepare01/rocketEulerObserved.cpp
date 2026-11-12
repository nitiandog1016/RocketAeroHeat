// 模块职责：复用原生成熟循环；仅观察完整状态/物理残差/边界账，不另写推进器。
#include "drivers/CSinglezoneDriver.hpp"
#include "iteration/CFluidIteration.hpp"
#include "variables/CEulerVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <limits>

namespace observation {
std::string folder;unsigned long last=0,saved=0;
void frame(CGeometry* g,CSolver* f,unsigned long step,bool computed){
  std::ostringstream name;name<<folder<<"/state_"<<std::setw(6)<<std::setfill('0')<<step<<".csv";
  std::ofstream out(name.str());out<<std::setprecision(17)<<"PointID,x,y,z,Density,Momentum_x,Momentum_y,Momentum_z,Energy,Volume,RHS0,RHS1,RHS2,RHS3,RHS4\n";
  auto* n=f->GetNodes();
  for(unsigned long i=0;i<g->GetnPoint();++i){out<<g->nodes->GetGlobalIndex(i);
    for(int d=0;d<3;++d)out<<','<<g->nodes->GetCoord(i,d);
    for(int j=0;j<5;++j)out<<','<<n->GetSolution(i,j);
    out<<','<<g->nodes->GetVolume(i);
    for(int j=0;j<5;++j)out<<','<<(computed?f->LinSysRes(i,j):std::numeric_limits<double>::quiet_NaN());out<<'\n';}
  if(!out)throw std::runtime_error("RAH_STATE_WRITE_FAILED");saved=step;
}
void observe(CGeometry* g,CSolver* f,CConfig* cfg){
  const auto step=cfg->GetInnerIter()+1;last=step;auto* n=f->GetNodes();
  double sum[5]={},sq[5]={},flux[5]={},wall[3]={},trunc=0;unsigned long zeroDt=0;
  for(unsigned long i=0;i<g->GetnPointDomain();++i){
    if(n->GetDelta_Time(i)==0)++zeroDt;
    double r[5];for(int j=0;j<5;++j){r[j]=-f->LinSysRes(i,j);trunc=std::max(trunc,std::abs(n->GetResTruncError(i)[j]));}
    r[4]-=300000*r[0];for(int j=0;j<5;++j){sum[j]+=r[j];sq[j]+=r[j]*r[j];}
  }
  for(unsigned m=0;m<cfg->GetnMarker_All();++m){const auto tag=cfg->GetMarker_All_TagBound(m);
    if(tag!="wall"&&tag!="farfield")continue;
    for(unsigned j=0;j<g->nVertex[m];++j){auto i=g->vertex[m][j]->GetNode();
      if(!g->nodes->GetDomain(i))continue;const auto* normal=g->vertex[m][j]->GetNormal();
      if(tag=="wall"){for(int d=0;d<3;++d)wall[d]-=n->GetPressure(i)*normal[d];continue;}
      const auto* b=f->GetCharacPrimVar(m,j);double un=0;for(int d=0;d<3;++d)un-=b[d+1]*normal[d];
      const double mass=b[5]*un;flux[0]+=mass;
      for(int d=0;d<3;++d)flux[d+1]+=mass*b[d+1]-b[4]*normal[d];
      flux[4]+=mass*(b[6]-300000);
    }
  }
  std::ofstream h(folder+"/physical_history.csv",step==1?std::ios::out:std::ios::app);h<<std::setprecision(17);
  if(step==1){h<<"updates,nodes,zeroDt,maxTrunc";for(auto label:{"rms","sum","outer"})for(int j=0;j<5;++j)h<<','<<label<<j;h<<",wallFx,wallFy,wallFz\n";}
  h<<step<<','<<g->GetnPointDomain()<<','<<zeroDt<<','<<trunc;
  for(int j=0;j<5;++j)h<<','<<std::sqrt(sq[j]/g->GetnPointDomain());
  for(auto value:sum)h<<','<<value;for(auto value:flux)h<<','<<value;for(auto value:wall)h<<','<<value;h<<'\n';
  if(!h)throw std::runtime_error("RAH_HISTORY_WRITE_FAILED");
  if(step==1||step%50==0)frame(g,f,step,true);
}
}

#include "rocketEulerFailure.hpp"
#include "EulerDevelopmentObservationV232.hpp"

// 模块二：只重载观测时机；父类Solve仍负责完整迭代、收敛和文件调度。
class ObservedIteration final : public CFluidIteration {
 public:
  using CFluidIteration::CFluidIteration;
  void Iterate(COutput* out,CIntegration**** integration,CGeometry**** geometry,CSolver***** solver,
    CNumerics****** numerics,CConfig** config,CSurfaceMovement** surface,CVolumetricMovement*** grid,
    CFreeFormDefBox*** box,unsigned short zone,unsigned short inst) override {
    observation::verify_native_prefix(config[zone]->GetInnerIter());
    observation::failure_phase="native-iterate";
    CFluidIteration::Iterate(out,integration,geometry,solver,numerics,config,surface,grid,box,zone,inst);
    observation::failure_phase="post-update-observation";
    observation::observe(geometry[zone][inst][0],solver[zone][inst][0][FLOW_SOL],config[zone]);
    observation::record_rejections_and_check(config[zone]);
    observation::failure_phase="native-monitor-output";
  }
};

class RocketEulerDriver final : public CSinglezoneDriver {
 public:
  using CSinglezoneDriver::CSinglezoneDriver;
  void launch(){
    auto* cfg=config_container[0];auto* g=geometry_container[0][0][0];auto* f=solver_container[0][0][0][FLOW_SOL];
    if(cfg->GetKind_FluidModel()!=RAH_PHASE_AIR_EULER||cfg->GetKind_Turb_Model()!=TURB_MODEL::NONE||
       cfg->GetnMGLevels()!=0||g->GetnPoint()!=18027||f->GetnVar()!=5||
       !dynamic_cast<CSU2PhaseAwareAir*>(f->GetFluidModel())||
       !dynamic_cast<CFluidIteration*>(iteration_container[0][0]))throw std::runtime_error("RAH_REAL_CASE_IDENTITY");
    // 单区静网格原iteration仅保存构造配置和计时；替换为只加观测的同父类。
    delete iteration_container[0][0];iteration_container[0][0]=new ObservedIteration(cfg);
    observation::failure_geometry=g;observation::failure_flow=f;observation::failure_config=cfg;
    observation::frame(g,f,0,false);
    observation::load_baseline();
    StartSolver();
    if(observation::saved!=observation::last)observation::frame(g,f,observation::last,true);
    std::ofstream result(observation::folder+"/native_completion.json");
    result<<"{\"normalReturn\":true,\"updates\":"<<observation::last<<",\"initialSaved\":true,\"physicalQualification\":false}\n";
  }
};

int main(int argc,char** argv){try{
  SU2_MPI::Init(&argc,&argv);if(argc!=3)throw std::runtime_error("CONFIG_AND_OUTPUT_REQUIRED");
  observation::folder=argv[2];std::set_terminate(observation::failed);
  RocketEulerDriver driver(argv[1],1,SU2_MPI::GetComm());driver.launch();
  observation::failure_phase="finalization";observation::failure_flow=nullptr;
  observation::failure_geometry=nullptr;observation::failure_config=nullptr;
  driver.Finalize();SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
