#include "numerics/flow/convection/hllc.hpp"
#include "PhaseStepAcceptanceV237.hpp"
#include "output/COutput.hpp"
// 模块职责：复用完整原生单区flow-SA循环；只读观察实际阶段，完整保存状态，不更改推进。
#include "drivers/CSinglezoneDriver.hpp"
#include "iteration/CFluidIteration.hpp"
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include "SingleDomainHooksV235.hpp"
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <limits>
#include <cmath>
#include <cstdint>

namespace observe {
std::string folder,phase="construction";CGeometry* geo=nullptr;CSolver** sol=nullptr;CConfig* cfg=nullptr;
unsigned long completed=0,saved=~0ul,interval=100,rejections=0;bool enabled=false;
std::ofstream ledger,events,steps,walls;std::vector<std::array<double,6>> before;
const char* columns="globalID,x,y,z,volume,rho,rhou,rhov,rhow,rhoE_internal,nuTilde,T,u,v,w,p,rho_primitive,H_internal,a,mu,mut_flow,k,Cp,mut_SA,wallDistance,flowRelaxation,flowCFL,flowRHS0,flowRHS1,flowRHS2,flowRHS3,flowRHS4_internal,saRHS";
void state_frame(unsigned long step){
  std::ostringstream name;name<<folder<<"/state_"<<std::setw(6)<<std::setfill('0')<<step<<".bin";
  const std::string tmp=name.str()+".tmp";std::ofstream f(tmp,std::ios::binary);
  const char magic[8]={'R','A','H','2','3','5','L','E'};const uint64_t N=geo->GetnPoint(),K=34;
  f.write(magic,8);f.write(reinterpret_cast<const char*>(&N),8);f.write(reinterpret_cast<const char*>(&K),8);
  auto* n=sol[FLOW_SOL]->GetNodes();auto* t=sol[TURB_SOL]->GetNodes();
  for(unsigned long i=0;i<N;++i){std::array<double,34> r{};int j=0;r[j++]=geo->nodes->GetGlobalIndex(i);
    for(int d=0;d<3;++d)r[j++]=geo->nodes->GetCoord(i,d);r[j++]=geo->nodes->GetVolume(i);
    for(int d=0;d<5;++d)r[j++]=n->GetSolution(i,d);r[j++]=t->GetSolution(i,0);
    for(int d=0;d<12;++d)r[j++]=n->GetPrimitive(i,d);
    r[j++]=t->GetmuT(i);r[j++]=geo->nodes->GetWall_Distance(i);r[j++]=n->GetUnderRelaxation(i);r[j++]=n->GetLocalCFL(i);
    for(int d=0;d<5;++d)r[j++]=step?sol[FLOW_SOL]->LinSysRes(i,d):std::numeric_limits<double>::quiet_NaN();
    r[j++]=step?sol[TURB_SOL]->LinSysRes(i,0):std::numeric_limits<double>::quiet_NaN();
    if(j!=33)throw std::runtime_error("FRAME_SCHEMA_COUNT");
    // 第34列专门标识末尾保留值，避免以后静默复用位置。
    r[j]=0;f.write(reinterpret_cast<const char*>(r.data()),sizeof(r));
  }
  f.close();if(!f||std::rename(tmp.c_str(),name.str().c_str()))throw std::runtime_error("FRAME_ATOMIC_WRITE");saved=step;
}

// 模块一：空间RHS先逐节点转换物理能量，再求和/范数；相邻阶段差属于同一次原生评估。
void stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short system){
  phase=label;const auto which=c->GetContainerPosition(system);const unsigned count=which==FLOW_SOL?5:1;
  const bool afterTime=std::string(label)=="after_time";
  const auto N=g->GetnPointDomain();std::array<long double,6> sums{},squares{};auto* n=s[which]->GetNodes();
  if(std::string(label)=="before_time"){
    before.resize(N);for(unsigned long i=0;i<N;++i){for(int j=0;j<5;++j)before[i][j]=s[FLOW_SOL]->GetNodes()->GetSolution(i,j);before[i][5]=s[TURB_SOL]->GetNodes()->GetSolution(i,0);}
  }
  double maxDelta=0,minRelax=1,maxCFL=0;
  for(unsigned long i=0;i<N;++i){
    for(unsigned j=0;j<count;++j){double r=s[which]->LinSysRes(i,j);if(which==FLOW_SOL&&j==4)r-=300000*s[which]->LinSysRes(i,0);
      const auto k=which==FLOW_SOL?j:5;sums[k]+=r;squares[k]+=(long double)r*r;}
    minRelax=std::min(minRelax,n->GetUnderRelaxation(i));maxCFL=std::max(maxCFL,n->GetLocalCFL(i));
    if(afterTime)for(unsigned j=0;j<count;++j){const auto k=which==FLOW_SOL?j:5;double now=n->GetSolution(i,j);
      if(!std::isfinite(now))throw std::runtime_error("NONFINITE_NATIVE_UPDATE");maxDelta=std::max(maxDelta,std::abs(now-before[i][k])/std::max(1.,std::abs(before[i][k])));}
  }
  ledger<<c->GetInnerIter()+1<<','<<which<<','<<label<<','<<N;
  for(auto v:sums)ledger<<','<<v;for(auto v:squares)ledger<<','<<std::sqrt(v/N);
  ledger<<','<<minRelax<<','<<maxCFL<<','<<maxDelta<<'\n';ledger.flush();if(!ledger)throw std::runtime_error("STAGE_WRITE");
  if(which==FLOW_SOL&&(std::string(label)=="before_strong"||std::string(label)=="after_strong")){
    // 同阶段原始量和梯度的压力/牛顿黏性牵引仅作观察，不回写求解器载荷缓存。
    std::array<long double,3> pressure{},shear{},row{};long double qLam=0,qEff=0;double wallSA=0;
    for(unsigned m=0;m<c->GetnMarker_All();++m)if(c->GetMarker_All_KindBC(m)==ISOTHERMAL)
      for(unsigned v=0;v<g->nVertex[m];++v){auto i=g->vertex[m][v]->GetNode(),j=g->vertex[m][v]->GetNormal_Neighbor();
        auto grad=n->GetGradient_Primitive(i);const auto* normal=g->vertex[m][v]->GetNormal();double A=0,projection=0,div=0;
        for(int d=0;d<3;++d){A+=normal[d]*normal[d];projection+=normal[d]*(g->nodes->GetCoord(j,d)-g->nodes->GetCoord(i,d));div+=grad[d+1][d];}
        A=std::sqrt(A);const double dist=std::abs(projection)/A;
        const double mu=n->GetLaminarViscosity(i)+n->GetEddyViscosity(i),k=n->GetThermalConductivity(i);
        const double eff=k+n->GetSpecificHeatCp(i)*n->GetEddyViscosity(i)/.9;
        const double dt=n->GetTemperature(j)-c->GetIsothermal_Temperature(c->GetMarker_All_TagBound(m));
        qLam+=k*dt/dist*A;qEff+=eff*dt/dist*A;wallSA=std::max(wallSA,std::abs(s[TURB_SOL]->GetNodes()->GetSolution(i,0)));
        for(int d=0;d<3;++d){pressure[d]-=n->GetPressure(i)*normal[d];row[d]+=s[FLOW_SOL]->LinSysRes(i,d+1);
          for(int e=0;e<3;++e)shear[d]-=mu*(grad[d+1][e]+grad[e+1][d]-(d==e?2./3*div:0))*normal[e];}
      }
    walls<<c->GetInnerIter()+1<<','<<label;for(auto x:pressure)walls<<','<<x;for(auto x:shear)walls<<','<<x;for(auto x:row)walls<<','<<x;
    walls<<','<<qLam<<','<<qEff<<','<<wallSA<<'\n';walls.flush();if(!walls)throw std::runtime_error("WALL_STAGE_WRITE");
  }
}

// 模块二：原CNS回退发生时记录当时的SA、旧/试探/实际返回态与原面可追溯的全局点号。
void rejection(CNSVariable* n,unsigned long i,const double* trial,const double* old,bool recovered,const char* reason){
  if(n!=sol[FLOW_SOL]->GetNodes()||i>=geo->GetnPoint())throw std::runtime_error("REJECTION_NODE_IDENTITY");
  if(++rejections>100000)throw std::runtime_error("REJECTION_OUTPUT_LIMIT");
  events<<cfg->GetInnerIter()+1<<','<<phase<<','<<i<<','<<geo->nodes->GetGlobalIndex(i)<<','<<recovered<<','<<reason;
  for(auto p:{trial,old})for(int j=0;j<5;++j)events<<','<<p[j];
  for(int j=0;j<5;++j)events<<','<<n->GetSolution(i,j);
  events<<','<<sol[TURB_SOL]->GetNodes()->GetSolution(i,0)<<','<<sol[TURB_SOL]->GetNodes()->GetSolution_Old(i,0)<<'\n';events.flush();
  if(!events)throw std::runtime_error("REJECTION_WRITE");
}

void end_iteration(){
  completed=cfg->GetInnerIter()+1;double umax=0,ntmax=0,tmin=1e100,tmax=0,mutdiff=0;
  auto* n=sol[FLOW_SOL]->GetNodes();auto* t=sol[TURB_SOL]->GetNodes();
  for(unsigned long i=0;i<geo->GetnPoint();++i){
    tmin=std::min(tmin,n->GetTemperature(i));tmax=std::max(tmax,n->GetTemperature(i));
    if(!std::isfinite(t->GetSolution(i,0)))throw std::runtime_error("NONFINITE_SA_ACCEPTED");
    mutdiff=std::max(mutdiff,std::abs(n->GetEddyViscosity(i)-t->GetmuT(i)));
  }
  for(unsigned m=0;m<cfg->GetnMarker_All();++m)if(cfg->GetMarker_All_KindBC(m)==ISOTHERMAL)
    for(unsigned v=0;v<geo->nVertex[m];++v){auto i=geo->vertex[m][v]->GetNode();
      for(int d=0;d<3;++d)umax=std::max(umax,std::abs(n->GetVelocity(i,d)));ntmax=std::max(ntmax,std::abs(t->GetSolution(i,0)));}
  steps<<completed<<','<<cfg->GetNonphysical_Points()<<','<<rejections<<','<<umax<<','<<ntmax<<','<<tmin<<','<<tmax<<','<<mutdiff<<'\n';steps.flush();
  if(!steps)throw std::runtime_error("STEP_WRITE");
  if(completed==1||completed%interval==0)state_frame(completed);
}
void failed(){std::string message="unknown";try{auto p=std::current_exception();if(p)std::rethrow_exception(p);}
  catch(const std::exception& e){message=e.what();}catch(...){message="non-standard";}
  std::ofstream f(folder+"/failure.txt");f<<"phase="<<phase<<"\ncompleted="<<completed<<"\nwhat="<<message<<'\n';f.flush();std::_Exit(86);}
}

// 模块三：父类Iterate完整包含flow、SA源项、边界、隐式更新及CFL；观测不替换任何计算环节。
class ObservedIterationV235 final:public CFluidIteration {
 public:using CFluidIteration::CFluidIteration;
 void Iterate(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,
   CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   observe::phase="native-iterate/preprocessing";
   CFluidIteration::Iterate(o,in,g,s,n,c,sm,gm,b,z,inst);
   observe::phase="accepted-flow-SA/end-iteration";observe::end_iteration();
 }
};
class HLLCDriverV238 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void launch(){using namespace observe;cfg=config_container[0];geo=geometry_container[0][0][0];sol=solver_container[0][0][0];
   if(cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||cfg->GetKind_Solver()!=MAIN_SOLVER::RANS||cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||
      cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!cfg->GetSAParsedOptions().ft2||cfg->GetnMGLevels()!=0||
      cfg->GetTime_Domain()||cfg->GetMultizone_Problem()||omp_get_max_threads()!=1||
      !dynamic_cast<CSU2PhaseAwareAir*>(sol[FLOW_SOL]->GetFluidModel())||sol[FLOW_SOL]->GetnVar()!=5||sol[TURB_SOL]->GetnVar()!=1)
     throw std::runtime_error("SINGLE_PHASE_SA_IDENTITY");
   // 检查实际工厂及同一气体/SA组合，不靠配置字符串推断派发。
   auto*** num=numerics_container[0][0][0];
   if(cfg->GetKind_Upwind_Flow()!=UPWIND::HLLC||geo->GetnPoint()!=24||cfg->GetnInner_Iter()!=5||
      !dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_TERM])||
      !dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_BOUND_TERM]))
      throw std::runtime_error("HLLC_V238_NATIVE_DISPATCH");
   std::ofstream identity(folder+"/native_identity.json");
   identity<<"{\"schemeID\":\"phase-air-general-hllc-density-jacobian-fix-v238\",\"physicalGasID\":15,\"actualInteriorBoundaryHLLC\":true,\"sameSAAirIdentity\":true,\"physicalQualification\":false}\n";identity.close();
   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");steps.open(folder+"/steps.csv");walls.open(folder+"/wall_stages.csv");
   for(auto* f:{&ledger,&events,&steps,&walls})*f<<std::setprecision(17);
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;
   ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   events<<"update,phase,localID,globalID,recovered,reason";for(auto word:{"trial","old","returned"})for(int j=0;j<5;++j)events<<','<<word<<j;events<<",SAcurrent,SAold\n";
   steps<<"update,lastFlowRejectionCount,totalObservedRejections,maxWallSpeedComponent,maxWallSA,minT,maxT,maxFlowVsUpdatedSAMuTDifference\n";
   delete iteration_container[0][0];iteration_container[0][0]=new ObservedIterationV235(cfg);
   if(enabled){rah_single_observer::stage=stage;rah_single_observer::rejection=rejection;}
   std::ofstream schema(folder+"/schema.json");schema<<"{\"format\":\"RAH235LE\",\"columns\":\""<<columns<<",reserved0\",\"initialRHSUnknown\":true,\"internalEnergyOffsetJkg\":300000,\"flowPrimitiveMutBeforeSACurrentRefresh\":true,\"physicalQualification\":false}\n";schema.close();
   state_frame(0);StartSolver();if(saved!=completed)state_frame(completed);
   // 显式原生终点保存，不设置或伪造收敛标志。
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);
   if(!forced||completed!=5||rah_phase_step_v237::checked!=5*24)throw std::runtime_error("HLLC_V238_COMPLETE_CONTROL");
   std::ofstream acceptance(folder+"/acceptance_summary.json");
   acceptance<<"{\"checked\":"<<rah_phase_step_v237::checked<<",\"limited\":"<<rah_phase_step_v237::limited<<",\"forcedFinalNativeOutput\":true}\n";acceptance.close();
   rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;
   std::ofstream f(folder+"/native_completion.json");f<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"initialSaved\":true,\"observationEnabled\":"<<(enabled?"true":"false")<<",\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{
  if(argc!=5)throw std::runtime_error("CONFIG_OUTPUT_OBSERVATION_INTERVAL_REQUIRED");
  observe::folder=argv[2];observe::enabled=std::string(argv[3])=="observe";observe::interval=std::stoul(argv[4]);
  if(!observe::interval||(std::string(argv[3])!="observe"&&std::string(argv[3])!="off"))throw std::runtime_error("INVALID_OBSERVER_MODE");
  std::set_terminate(observe::failed);SU2_MPI::Init(&argc,&argv);HLLCDriverV238 driver(argv[1],1,SU2_MPI::GetComm());driver.launch();driver.Finalize();SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
