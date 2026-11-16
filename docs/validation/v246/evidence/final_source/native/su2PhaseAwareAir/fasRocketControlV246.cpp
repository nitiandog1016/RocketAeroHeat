// 模块职责：复用细层完整观测与成熟FAS回调，执行真实域唯一五步控制。
#define main rah_unused_single_domain_main_v246
#include "singleDomainRunV235.cpp"
#undef main
#include "FASHooksV245.hpp"

namespace real246 {
CGeometry* coarse=nullptr;CSolver** cs=nullptr;std::ofstream levels;unsigned long transfers=0;
void coarse_frame(const char* name){
  auto* n=cs[FLOW_SOL]->GetNodes();auto* t=cs[TURB_SOL]->GetNodes();std::ofstream f(observe::folder+"/"+name,std::ios::binary);
  uint64_t N=coarse->GetnPoint(),K=33;f.write("RAH246LE",8);f.write(reinterpret_cast<char*>(&N),8);f.write(reinterpret_cast<char*>(&K),8);
  for(unsigned long i=0;i<N;++i){std::array<double,33> r{};int j=0;r[j++]=coarse->nodes->GetGlobalIndex(i);
    for(int d=0;d<3;++d)r[j++]=coarse->nodes->GetCoord(i,d);r[j++]=coarse->nodes->GetVolume(i);
    for(int d=0;d<5;++d)r[j++]=n->GetSolution(i,d);r[j++]=t->GetSolution(i,0);r[j++]=t->GetmuT(i);
    for(int d=0;d<12;++d)r[j++]=n->GetPrimitive(i,d);for(int d=0;d<8;++d)r[j++]=n->GetSecondary(i,d);
    r[j++]=coarse->nodes->GetWall_Distance(i);if(j!=K)throw std::runtime_error("COARSE_FRAME_SCHEMA");f.write(reinterpret_cast<char*>(r.data()),sizeof(r));
  }f.close();if(!f)throw std::runtime_error("COARSE_FRAME_WRITE");
}
// 各层分别检查当前守恒态，输出摘要；不把更新前primitive缓存当更新后状态。
void summarize(const char* label,CGeometry* g,CSolver** s){
  observe::phase=label;const unsigned l=g==observe::geo?0:(g==coarse?1:99);if(l==99)throw std::runtime_error("UNKNOWN_LEVEL");
  auto* n=s[FLOW_SOL]->GetNodes();std::array<long double,5> integrated{};double tmin=1e100,tmax=0;const auto N=g->GetnPoint();
  for(unsigned long i=0;i<N;++i){const auto* u=n->GetSolution(i);double v2=0;for(int d=1;d<4;++d)v2+=(u[d]/u[0])*(u[d]/u[0]);
    const auto q=rah_phase_air::state("rhoe",u[0],u[4]/u[0]-.5*v2);tmin=std::min(tmin,q.t);tmax=std::max(tmax,q.t);
    for(int d=0;d<5;++d)integrated[d]+=(u[d]-(d==4?300000*u[0]:0))*g->nodes->GetVolume(i);
    if(std::string(label)=="prolong_in")for(int d=0;d<5;++d)if(!std::isfinite(s[FLOW_SOL]->LinSysRes(i,d)))throw std::runtime_error("NONFINITE_PROLONGATION");
  }
  levels<<observe::cfg->GetInnerIter()+1<<','<<label<<','<<l<<','<<N<<','<<tmin<<','<<tmax;
  for(auto x:integrated)levels<<','<<x;levels<<'\n';levels.flush();if(!levels)throw std::runtime_error("LAYER_LEDGER_WRITE");
}
void stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short sys){
  if(g==observe::geo)observe::stage(label,g,s,c,sys);
  if(c->GetContainerPosition(sys)==FLOW_SOL&&std::string(label)=="after_time")summarize("flow_after_time",g,s);
}
void transfer(const char* label,CGeometry* f,CSolver** fs,CGeometry* c,CSolver** ss,CConfig*){++transfers;summarize(label,f,fs);summarize(label,c,ss);}
void reject(CNSVariable* v,unsigned long i,const double* trial,const double* old,bool recovered,const char* reason){
  const unsigned l=v==observe::sol[FLOW_SOL]->GetNodes()?0:(v==cs[FLOW_SOL]->GetNodes()?1:99);
  std::ofstream f(observe::folder+"/fas_rejection.txt");f<<std::setprecision(17)<<"level="<<l<<" point="<<i<<" phase="<<observe::phase<<" reason="<<reason<<" recovered="<<recovered;
  for(int d=0;d<5;++d)f<<"\ntrial"<<d<<'='<<trial[d]<<" old"<<d<<'='<<old[d];f<<'\n';f.close();throw std::runtime_error("FAS_REAL_STATE_REJECTED");
}
}
class FASRocketDriverV246 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void launch(){using namespace observe;cfg=config_container[0];geo=geometry_container[0][0][0];sol=solver_container[0][0][0];
   real246::coarse=geometry_container[0][0][1];real246::cs=solver_container[0][0][1];
   if(cfg->GetnMGLevels()!=1||cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||
      cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!cfg->GetSAParsedOptions().ft2||cfg->GetRestart()||omp_get_max_threads()!=1||
      geo->GetnPoint()!=335726||real246::coarse->GetnPoint()!=44484)throw std::runtime_error("REAL_FAS_IDENTITY");
   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");steps.open(folder+"/steps.csv");walls.open(folder+"/wall_stages.csv");real246::levels.open(folder+"/levels.csv");
   for(auto* f:{&ledger,&events,&steps,&walls,&real246::levels})*f<<std::setprecision(17);
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   steps<<"update,lastFlowRejectionCount,totalObservedRejections,maxWallSpeedComponent,maxWallSA,minT,maxT,maxFlowVsUpdatedSAMuTDifference\n";
   real246::levels<<"update,stage,level,nodes,minT,maxT,volumeRho,volumeMx,volumeMy,volumeMz,volumeEphysical\n";
   delete iteration_container[0][0];iteration_container[0][0]=new ObservedIterationV235(cfg);
   if(enabled){rah_single_observer::stage=real246::stage;rah_single_observer::rejection=real246::reject;rah_fas_observer::transfer=real246::transfer;}
   state_frame(0);real246::coarse_frame("coarse_initial.bin");StartSolver();if(saved!=completed)state_frame(completed);real246::coarse_frame("coarse_final.bin");
   rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;rah_fas_observer::transfer=nullptr;
   std::ofstream done(folder+"/native_completion.json");done<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"finePoints\":335726,\"coarsePoints\":44484,\"transfersObserved\":"<<real246::transfers<<",\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{
 if(argc!=4)throw std::runtime_error("CONFIG_FOLDER_MODE_REQUIRED");observe::folder=argv[2];const std::string mode=argv[3];
 if(mode!="observe"&&mode!="off")throw std::runtime_error("MODE");observe::enabled=mode=="observe";observe::interval=5;
 std::set_terminate(observe::failed);SU2_MPI::Init(&argc,&argv);FASRocketDriverV246 driver(argv[1],1,SU2_MPI::GetComm());driver.launch();driver.Finalize();SU2_MPI::Finalize();return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
