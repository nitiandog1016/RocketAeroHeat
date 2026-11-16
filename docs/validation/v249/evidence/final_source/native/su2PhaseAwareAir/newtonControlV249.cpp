// 模块职责：调用已有原生Newton-Krylov，在真实保存状态上完成一次完整flow-SA控制；不改数值内核。
#define main rah_unused_single_main_v249
#include "singleDomainRunV235.cpp"
#undef main
#include "integration/CNewtonIntegration.hpp"

namespace nk249 {
unsigned long flowEvaluations=0,saEvaluations=0;
void stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short sys){
  if(std::string(label)=="after_sources"){
    if(c->GetContainerPosition(sys)==FLOW_SOL)++flowEvaluations;
    if(c->GetContainerPosition(sys)==TURB_SOL)++saEvaluations;
  }
  observe::stage(label,g,s,c,sys);
}
void reject(CNSVariable* n,unsigned long i,const double* trial,const double* old,bool recovered,const char* reason){
  observe::rejection(n,i,trial,old,recovered,reason);
  // 任何试探恢复都会破坏该输入的差分函数身份；记录后拒绝本次控制，不能继续解释成有效Jv。
  throw std::runtime_error("NK_CONTROL_DOMAIN_RECOVERY_REJECTED");
}
void initial(const std::string& reference){
  using namespace observe;const auto N=geo->GetnPoint();std::vector<unsigned long> local(N,~0ul);
  for(unsigned long i=0;i<N;++i){auto id=geo->nodes->GetGlobalIndex(i);if(id>=N||local[id]!=~0ul)throw std::runtime_error("NK_GLOBAL_ID");local[id]=i;}
  std::ifstream f(reference,std::ios::binary);char magic[8];uint64_t count,k;
  f.read(magic,8);f.read(reinterpret_cast<char*>(&count),8);f.read(reinterpret_cast<char*>(&k),8);
  if(std::memcmp(magic,"RAH248LE",8)||count!=N||k!=34)throw std::runtime_error("NK_REFERENCE_FORMAT");
  unsigned long mismatch=0;double volumeError=0;
  for(unsigned long row=0;row<N;++row){std::array<double,34>a{};f.read(reinterpret_cast<char*>(a.data()),sizeof(a));const auto id=(unsigned long)a[0];if(id>=N||a[0]!=id)throw std::runtime_error("NK_REFERENCE_ID");const auto i=local[id];
    for(int j=0;j<3;++j){double x=geo->nodes->GetCoord(i,j);mismatch+=std::memcmp(&x,&a[j+1],8)!=0;}
    for(int j=0;j<6;++j){double u=j<5?sol[FLOW_SOL]->GetNodes()->GetSolution(i,j):sol[TURB_SOL]->GetNodes()->GetSolution(i,0);mismatch+=std::memcmp(&u,&a[j+5],8)!=0;}
    volumeError=std::max(volumeError,std::abs(geo->nodes->GetVolume(i)-a[4])/a[4]);
    if(sol[FLOW_SOL]->GetNodes()->GetLocalCFL(i)!=.1)throw std::runtime_error("NK_INITIAL_CFL");
  }
  if(!f||mismatch||volumeError>2e-14)throw std::runtime_error("NK_INITIAL_STATE_IDENTITY");
  std::ofstream out(folder+"/initial_check.json");out<<std::setprecision(17)<<"{\"passed\":true,\"sixStateCoordinateBitMismatches\":0,\"maximumVolumeRelativeDifference\":"<<volumeError<<",\"resetCFL\":0.1,\"sourceUpdate\":200,\"newMethodControl\":true}\n";
}
}
class NewtonControl249 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void launch(const std::string& reference){using namespace observe;cfg=config_container[0];geo=geometry_container[0][0][0];sol=solver_container[0][0][0];
   auto* engine=dynamic_cast<CNewtonIntegration*>(integration_container[0][0][FLOW_SOL]);auto ip=cfg->GetNewtonKrylovIntParam();auto dp=cfg->GetNewtonKrylovDblParam();
   if(!engine||!cfg->GetNewtonKrylov()||cfg->GetnMGLevels()!=0||cfg->GetnInner_Iter()!=1||ip[0]!=0||dp[0]!=0||dp[3]!=1e-4||dp[4]!=1||
      cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||
      !cfg->GetSAParsedOptions().ft2||!cfg->GetRestart()||geo->GetnPoint()!=335726||omp_get_max_threads()!=1)throw std::runtime_error("NK_NATIVE_IDENTITY");
   nk249::initial(reference);
   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");steps.open(folder+"/steps.csv");walls.open(folder+"/wall_stages.csv");
   for(auto* f:{&ledger,&events,&steps,&walls})*f<<std::setprecision(17);
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   events<<"update,phase,localID,globalID,recovered,reason";for(auto word:{"trial","old","returned"})for(int j=0;j<5;++j)events<<','<<word<<j;events<<",SAcurrent,SAold\n";
   steps<<"update,lastFlowRejectionCount,totalObservedRejections,maxWallSpeedComponent,maxWallSA,minT,maxT,maxFlowVsUpdatedSAMuTDifference\n";
   enabled=true;interval=1;delete iteration_container[0][0];iteration_container[0][0]=new ObservedIterationV235(cfg);
   rah_single_observer::stage=nk249::stage;rah_single_observer::rejection=nk249::reject;
   state_frame(0);StartSolver();if(saved!=completed)state_frame(completed);
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);
   if(!forced||completed!=1||nk249::flowEvaluations<2||nk249::saEvaluations!=1||rejections)throw std::runtime_error("NK_COMPLETE_CONTROL");
   rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;
   std::ofstream out(folder+"/native_completion.json");out<<"{\"normalReturn\":true,\"updates\":1,\"actualNewtonFactory\":true,\"flowSpaceEvaluations\":"<<nk249::flowEvaluations<<",\"SAEvaluations\":"<<nk249::saEvaluations<<",\"forcedNativeFinal\":true,\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{
 if(argc!=4)throw std::runtime_error("CONFIG_FOLDER_REFERENCE_REQUIRED");observe::folder=argv[2];std::set_terminate(observe::failed);SU2_MPI::Init(&argc,&argv);
 NewtonControl249 driver(argv[1],1,SU2_MPI::GetComm());driver.launch(argv[3]);driver.Finalize();SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::ofstream out(observe::folder+"/caught_failure.txt");out<<"phase="<<observe::phase<<"\ncompleted="<<observe::completed<<"\nwhat="<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}}
