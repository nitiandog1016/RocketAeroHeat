// 模块职责：完整原域HLLC均匀初场、固定收益停止与实际接受/失败观察；成熟推进循环保持。
#include "commonBranchV238.hpp"
#include "PhaseStepObservationV238.hpp"
#include "numerics/flow/convection/hllc.hpp"
namespace events238 {
using State=rah_phase_step_v237::State;
struct Pending {unsigned long update,point;State accepted;};
std::vector<Pending> pending;std::ofstream candidates,actual,failures;
void failure(unsigned long update,unsigned long point,const double* old,const double* delta,double upper,const char* why){
 failures<<update<<','<<point<<',';
 if(observe::geo)failures<<observe::geo->nodes->GetGlobalIndex(point);else failures<<-1;
 failures<<','<<upper<<','<<why;
 for(int d=0;d<3;++d){failures<<',';if(observe::geo)failures<<observe::geo->nodes->GetCoord(point,d);else failures<<"nan";}
 for(const auto* a:{old,delta})for(int j=0;j<5;++j)failures<<','<<a[j];failures<<'\n';failures.flush();
}
void candidate(unsigned long update,unsigned long point,const State& old,const State& delta,double upper,const rah_phase_step_v237::Result& r){
 if(rah_phase_step_v237::limited>100000)throw std::runtime_error("HLLC_EVENT_CAPACITY");
 candidates<<update<<','<<point<<','<<observe::geo->nodes->GetGlobalIndex(point)<<','<<upper<<','<<r.alpha<<','<<r.reductions;
 for(const State* a:std::array<const State*,3>{&old,&delta,&r.accepted})for(double x:*a)candidates<<','<<x;
 candidates<<'\n';candidates.flush();if(!candidates)throw std::runtime_error("CANDIDATE_WRITE");
 pending.push_back({update,point,r.accepted});
}
void stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short system){
 branch236::stage236(label,g,s,c,system);
 if(std::string(label)!="after_time"||c->GetContainerPosition(system)!=FLOW_SOL)return;
 for(const auto& p:pending){
   if(p.update!=c->GetInnerIter()+1)throw std::runtime_error("ACTUAL_ACCEPTANCE_STAGE");
   actual<<p.update<<','<<p.point<<','<<g->nodes->GetGlobalIndex(p.point);
   for(double x:p.accepted)actual<<','<<x;
   for(int j=0;j<5;++j)actual<<','<<s[FLOW_SOL]->GetNodes()->GetSolution(p.point,j);actual<<'\n';
 }
 actual.flush();if(!actual)throw std::runtime_error("ACTUAL_ACCEPTANCE_WRITE");pending.clear();
}
void open(const std::string& folder){
 candidates.open(folder+"/accepted_steps.csv");actual.open(folder+"/actual_accepted_steps.csv");failures.open(folder+"/failed_step.csv");
 for(auto* f:{&candidates,&actual,&failures})*f<<std::setprecision(17);
 candidates<<"update,localID,globalID,nativeUpper,acceptedAlpha,reductions";
 for(auto prefix:{"old","delta","accepted"})for(int j=0;j<5;++j)candidates<<','<<prefix<<j;candidates<<'\n';
 actual<<"update,localID,globalID";for(auto prefix:{"candidate","actual"})for(int j=0;j<5;++j)actual<<','<<prefix<<j;actual<<'\n';
 failures<<"update,localID,globalID,nativeUpper,reason,x,y,z";for(auto prefix:{"old","delta"})for(int j=0;j<5;++j)failures<<','<<prefix<<j;failures<<'\n';
 rah_phase_step_v237::observer=candidate;rah_phase_step_v238::failure=failure;
}
}

// 模块一：同一次真实初始化逐节点检查相容均匀六量和CFL；原面索引完整保持。
void uniform_initial(const std::string& expectedFile,const std::string& facefile){
 using namespace observe;using namespace branch236;const auto N=geo->GetnPoint();
 if(N!=335726||geo->GetnPointDomain()!=N||cfg->GetRestart())throw std::runtime_error("UNIFORM_DOMAIN_IDENTITY");
 std::array<double,6> expected{};std::ifstream ef(expectedFile);for(auto& v:expected)if(!(ef>>v))throw std::runtime_error("UNIFORM_SIX_INPUT");
 local.assign(N,~0ul);unsigned long mismatches=0;double minCFL=1e100,maxCFL=0;
 for(unsigned long i=0;i<N;++i){auto id=geo->nodes->GetGlobalIndex(i);if(id>=N||local[id]!=~0ul)throw std::runtime_error("UNIFORM_GLOBAL_ID");local[id]=i;
   for(int j=0;j<6;++j){double x=j<5?sol[FLOW_SOL]->GetNodes()->GetSolution(i,j):sol[TURB_SOL]->GetNodes()->GetSolution(i,0);mismatches+=std::memcmp(&x,&expected[j],8)!=0;}
   double cfl=sol[FLOW_SOL]->GetNodes()->GetLocalCFL(i);minCFL=std::min(minCFL,cfl);maxCFL=std::max(maxCFL,cfl);
 }
 std::ofstream report(folder+"/initial_check.json");report<<std::setprecision(17)<<"{\"points\":"<<N<<",\"uniformSixStateBitMismatches\":"<<mismatches<<",\"flowCFLMin\":"<<minCFL<<",\"flowCFLMax\":"<<maxCFL<<",\"restart\":false,\"initialRHSUnknown\":true}\n";report.close();
 if(mismatches||minCFL!=.1||maxCFL!=.1)throw std::runtime_error("UNIFORM_INITIAL_STATE_MISMATCH");
 std::ifstream ff(facefile);Face face;unsigned long id=0;double total=0;
 while(ff>>id>>face.nodes[0]>>face.nodes[1]>>face.nodes[2]>>face.area){
   if(id!=faces.size()+1||!(face.area>0))throw std::runtime_error("ORIGINAL_FACE_ID");
   for(auto& n:face.nodes){if(n>=N)throw std::runtime_error("ORIGINAL_FACE_NODE");n=local[n];}total+=face.area;faces.push_back(face);
 }
 if(faces.size()!=19752||std::abs(total-1.9646292565715553)>1e-12)throw std::runtime_error("ORIGINAL_SURFACE_COVERAGE");
}

class UniformHLLCDriverV238 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void launch(const std::string& expected,const std::string& faces){using namespace observe;
   cfg=config_container[0];geo=geometry_container[0][0][0];sol=solver_container[0][0][0];auto*** num=numerics_container[0][0][0];
   if(cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||cfg->GetKind_Upwind_Flow()!=UPWIND::HLLC||cfg->GetKind_Solver()!=MAIN_SOLVER::RANS||
      cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!cfg->GetSAParsedOptions().ft2||
      cfg->GetnMGLevels()!=0||cfg->GetTime_Domain()||cfg->GetMultizone_Problem()||omp_get_max_threads()!=1||
      !dynamic_cast<CSU2PhaseAwareAir*>(sol[FLOW_SOL]->GetFluidModel())||sol[FLOW_SOL]->GetnVar()!=5||sol[TURB_SOL]->GetnVar()!=1||cfg->GetnInner_Iter()!=400||
      !dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_TERM])||!dynamic_cast<CUpwGeneralHLLC_Flow*>(num[FLOW_SOL][CONV_BOUND_TERM]))
      throw std::runtime_error("UNIFORM_HLLC_V238_IDENTITY");
   phase="uniform-initial-check";uniform_initial(expected,faces);
   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");steps.open(folder+"/steps.csv");walls.open(folder+"/wall_stages.csv");
   branch236::trends.open(folder+"/surface_trend.csv");branch236::decisions.open(folder+"/benefit_decisions.csv");
   for(auto* f:{&ledger,&events,&steps,&walls,&branch236::trends,&branch236::decisions})*f<<std::setprecision(17);
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   events<<"update,phase,localID,globalID,recovered,reason";for(auto word:{"trial","old","returned"})for(int j=0;j<5;++j)events<<','<<word<<j;events<<",SAcurrent,SAold\n";
   steps<<"update,lastFlowRejectionCount,totalObservedRejections,maxWallSpeedComponent,maxWallSA,minT,maxT,maxFlowVsUpdatedSAMuTDifference\n";
   branch236::trends<<"update,pressureAreaL1Change,heatAreaL1Change\n";
   branch236::decisions<<"check,windowEnd,ratio0,ratio1,ratio2,ratio3,ratio4,ratio5,geometricMeanRatio,oldPressureChange,oldHeatChange,newPressureChange,newHeatChange,noSignificantDeterioration,benefit\n";
   std::ofstream schema(folder+"/schema.json");schema<<"{\"format\":\"RAH235LE\",\"columns\":\""<<columns<<",reserved0\",\"internalEnergyOffsetJkg\":300000,\"physicalQualification\":false,\"branch\":\"phase-air-general-hllc-uniform-v238\"}\n";schema.close();
   delete iteration_container[0][0];iteration_container[0][0]=new BranchIterationV236(cfg);
   rah_single_observer::stage=events238::stage;rah_single_observer::rejection=rejection;enabled=true;interval=100;
   state_frame(0);branch236::wall_frame(0);StartSolver();if(saved!=completed)state_frame(completed);if(completed%50)branch236::wall_frame(completed);
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);
   if(!forced||!events238::pending.empty())throw std::runtime_error("HLLC_FINAL_SAVE_OR_ACCEPTANCE");
   std::ofstream accept(folder+"/acceptance_summary.json");accept<<"{\"checked\":"<<rah_phase_step_v237::checked<<",\"limited\":"<<rah_phase_step_v237::limited<<",\"forcedFinalNativeOutput\":true}\n";
   rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;
   std::ofstream f(folder+"/native_completion.json");f<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"stopReason\":\""<<branch236::stopReason<<"\",\"initialSaved\":true,\"uniformFreestreamInitial\":true,\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{
 if(argc==3&&std::string(argv[1])=="--failure-control"){
   observe::folder=argv[2];events238::open(argv[2]);const auto q=rah_phase_air::PT(1,60.0001);
   double u[5]={q.rho,0,0,0,q.rho*q.e},d[5]={0,0,0,0,-1e8};bool refused=false;
   try{rah_phase_step_v238::native_limit(3,7,u,d,1);}catch(const std::domain_error& e){refused=std::string(e.what())=="PHASE_STEP_NO_ADMISSIBLE_PROGRESS";}
   if(!refused)throw std::runtime_error("FAILURE_OBSERVER_CONTROL");return 0;
 }
 if(argc!=5)throw std::runtime_error("CONFIG_OUTPUT_UNIFORM_FACES_REQUIRED");observe::folder=argv[2];std::set_terminate(observe::failed);events238::open(argv[2]);
 SU2_MPI::Init(&argc,&argv);UniformHLLCDriverV238 driver(argv[1],1,SU2_MPI::GetComm());driver.launch(argv[3],argv[4]);driver.Finalize();SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::ofstream f(observe::folder+"/caught_failure.txt");f<<"phase="<<observe::phase<<"\ncompleted="<<observe::completed<<"\nwhat="<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}}
