// 模块职责：固定原真实两层FAS的一次有界发展；仅保存、观察和停止，不改成熟求解公式。
#define main rah_unused_v235_main_248
#define state_frame legacy_state_frame
#define end_iteration legacy_end_iteration
#include "singleDomainRunV235.cpp"
#undef main
#undef state_frame
#undef end_iteration
#include "FASHooksV245.hpp"
#include "output/COutput.hpp"
#include "integration/CIntegration.hpp"
#include <algorithm>
namespace observe {
void state_frame(unsigned long step){
  std::ostringstream name;name<<folder<<"/state_"<<std::setw(6)<<std::setfill('0')<<step<<".bin";
  const std::string tmp=name.str()+".tmp";std::ofstream f(tmp,std::ios::binary);
  const char magic[8]={'R','A','H','2','4','8','L','E'};const uint64_t N=geo->GetnPoint(),K=34;
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
}
namespace dev248 {
using namespace observe;
struct Face {std::array<unsigned long,3> nodes;double area;};
std::vector<Face> faces;
std::vector<unsigned long> local;
std::vector<std::array<double,6>> logs;
std::vector<std::array<double,2>> previousFaces;
std::vector<std::array<double,2>> changes;
std::array<double,6> currentRMS{};
bool stop=false;std::string stopReason="maximum-updates";
std::ofstream trends,decisions;

// 模块二：壁面小表保存原生缓存热流与当前primitive的壁法向梯度，两者保持不同阶段身份。
void wall_frame(unsigned long step){
 std::ostringstream path;path<<folder<<"/wall_"<<std::setw(6)<<std::setfill('0')<<step<<".csv";
 std::ofstream f(path.str());f<<std::setprecision(17)<<"globalID,pressurePa,wallFluidK,qCurrentGradientWm2,qNativeCacheWm2,yPlusNative,wallDistance,thermalDistance,Tneighbor,k,mutFlow,mutSA,cfX,cfY,cfZ\n";
 auto* n=sol[FLOW_SOL]->GetNodes();std::vector<double> q(geo->GetnPoint(),std::numeric_limits<double>::quiet_NaN());
 unsigned count=0;
 for(unsigned m=0;m<cfg->GetnMarker_All();++m)if(cfg->GetMarker_All_KindBC(m)==ISOTHERMAL)
 for(unsigned v=0;v<geo->nVertex[m];++v){auto i=geo->vertex[m][v]->GetNode(),j=geo->vertex[m][v]->GetNormal_Neighbor();
   auto normal=geo->vertex[m][v]->GetNormal();double a2=0,proj=0;
   for(int d=0;d<3;++d){a2+=normal[d]*normal[d];proj+=normal[d]*(geo->nodes->GetCoord(j,d)-geo->nodes->GetCoord(i,d));}
   const double dist=std::abs(proj)/std::sqrt(a2),k=n->GetThermalConductivity(i),Tj=n->GetTemperature(j);
   q[i]=k*(Tj-300)/dist;if(!std::isfinite(q[i])||!(dist>0))throw std::runtime_error("WALL_DIAGNOSTIC_NONFINITE");
   f<<geo->nodes->GetGlobalIndex(i)<<','<<n->GetPressure(i)<<','<<n->GetTemperature(i)<<','<<q[i]<<',';
   // 初始缓存未经该分支推进，不把构造器默认值当作原生热诊断。
   if(step)f<<sol[FLOW_SOL]->GetHeatFlux(m,v);else f<<"nan";
   f<<',';if(step)f<<sol[FLOW_SOL]->GetYPlus(m,v);else f<<"nan";
   f<<','<<geo->nodes->GetWall_Distance(i)<<','<<dist<<','<<Tj<<','<<k<<','<<n->GetEddyViscosity(i)<<','<<sol[TURB_SOL]->GetNodes()->GetmuT(i);
   for(int d=0;d<3;++d){f<<',';if(step)f<<sol[FLOW_SOL]->GetCSkinFriction(m,v,d);else f<<"nan";}f<<'\n';++count;
 }
 if(count!=9878)throw std::runtime_error("ORIGINAL_WALL_COUNT");
 long double diffP=0,diffQ=0,absP=0,absQ=0;std::vector<std::array<double,2>> now;
 for(const auto& face:faces){double p=0,h=0;for(auto i:face.nodes){p+=n->GetPressure(i)/3.;h+=q[i]/3.;}if(!std::isfinite(p+h))throw std::runtime_error("FACE_NOT_ON_WALL");
   if(!previousFaces.empty()){const auto old=previousFaces[now.size()];diffP+=face.area*std::abs(p-old[0]);diffQ+=face.area*std::abs(h-old[1]);absP+=face.area*std::abs(old[0]);absQ+=face.area*std::abs(old[1]);}now.push_back({p,h});}
 std::array<double,2> change{double(diffP/std::max(1.L,absP)),double(diffQ/std::max(1.L,absQ))};changes.push_back(change);previousFaces.swap(now);
 trends<<step<<','<<change[0]<<','<<change[1]<<'\n';trends.flush();f.close();if(!f||!trends)throw std::runtime_error("WALL_TREND_WRITE");
}

void fine_stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short system){
 observe::stage(label,g,s,c,system);
 if(std::string(label)!="after_time")return;
 const auto which=c->GetContainerPosition(system);const auto N=g->GetnPointDomain();
 for(unsigned j=0;j<(which==FLOW_SOL?5:1);++j){long double sum=0;
   for(unsigned long i=0;i<N;++i){double r=s[which]->LinSysRes(i,j);if(which==FLOW_SOL&&j==4)r-=300000*s[which]->LinSysRes(i,0);if(!std::isfinite(r))throw std::runtime_error("RHS_NONFINITE");sum+=(long double)r*r;}
   currentRMS[which==FLOW_SOL?j:5]=std::sqrt(sum/N);
 }
}

// 模块三：固定两窗收益规则；仅决定是否退出原循环，不调整CFL、解或残差。
void benefit(){
 std::array<double,6> l{};for(int j=0;j<6;++j)l[j]=std::log10(std::max(1e-300,currentRMS[j]));logs.push_back(l);
 if(completed%50==0)wall_frame(completed);
 if(completed!=200&&completed!=300)return;
 unsigned noBenefit=0;bool worse=false;
 for(unsigned end:{unsigned(completed-50),unsigned(completed)}){
   std::array<double,6> ratio{},later{};double logRatio=0;bool allSmall=true,notWorse=true;
   for(unsigned j=0;j<6;++j){double a=0,b=0;for(unsigned k=end-100;k<end-50;++k)a+=logs[k][j]/50.;for(unsigned k=end-50;k<end;++k)b+=logs[k][j]/50.;ratio[j]=std::pow(10.,b-a);later[j]=std::pow(10.,b);logRatio+=b-a;notWorse=notWorse&&ratio[j]<=1.25;allSmall=allSmall&&b<=-8.;}
   const auto old=changes.at(end/50-1),now=changes.at(end/50);bool areaGood=true;
   for(unsigned j=0;j<2;++j){const double threshold=j==0?.01:.05;notWorse=notWorse&&now[j]<=1.25*std::max(old[j],threshold);areaGood=areaGood&&(now[j]<=threshold||now[j]<=.95*old[j]);}
   const double gm=std::pow(10.,logRatio/6.);const bool good=notWorse&&areaGood&&(gm<=.95||allSmall);noBenefit+=!good;worse=worse||!notWorse;
   decisions<<completed<<','<<end;for(auto x:ratio)decisions<<','<<x;decisions<<','<<gm<<','<<old[0]<<','<<old[1]<<','<<now[0]<<','<<now[1]<<','<<notWorse<<','<<good<<'\n';
 }
 stop=worse||noBenefit==2;if(stop)stopReason=worse?"significant-deterioration":"two-windows-no-benefit";decisions.flush();if(!decisions)throw std::runtime_error("DECISION_WRITE");
}
}

namespace dev248 {
extern CGeometry* coarse;extern CSolver** cs;
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
}

namespace dev248 {
CGeometry* coarse=nullptr;CSolver** cs=nullptr;unsigned long transfers=0;
bool terminal=false;double started=omp_get_wtime();
void initial(const std::string& reference,const std::string& facefile){
  const auto N=geo->GetnPoint();local.assign(N,~0ul);
  for(unsigned long i=0;i<N;++i){auto id=geo->nodes->GetGlobalIndex(i);if(id>=N||local[id]!=~0ul)throw std::runtime_error("GLOBAL_ID_LAYOUT");local[id]=i;}
  std::ifstream f(reference,std::ios::binary);char magic[8];uint64_t count,k;f.read(magic,8);f.read(reinterpret_cast<char*>(&count),8);f.read(reinterpret_cast<char*>(&k),8);
  if(std::memcmp(magic,"RAH235LE",8)||count!=N||k!=34)throw std::runtime_error("UNIFORM_REFERENCE_FORMAT");
  unsigned long mismatch=0;
  for(unsigned long p=0;p<N;++p){std::array<double,34>a{};f.read(reinterpret_cast<char*>(a.data()),sizeof(a));const auto id=(unsigned long)a[0];if(id>=N||a[0]!=id)throw std::runtime_error("INITIAL_GLOBAL_ID");auto i=local[id];
    std::array<double,11>b{};b[0]=id;for(int d=0;d<3;++d)b[d+1]=geo->nodes->GetCoord(i,d);b[4]=geo->nodes->GetVolume(i);
    for(int d=0;d<5;++d)b[d+5]=sol[FLOW_SOL]->GetNodes()->GetSolution(i,d);b[10]=sol[TURB_SOL]->GetNodes()->GetSolution(i,0);
    for(int d=0;d<11;++d)mismatch+=std::memcmp(&a[d],&b[d],8)!=0;
  }
  if(!f||mismatch)throw std::runtime_error("UNIFORM_INITIAL_BIT_MISMATCH");
  std::ifstream ff(facefile);Face face;unsigned long id=0;double area=0;
  while(ff>>id>>face.nodes[0]>>face.nodes[1]>>face.nodes[2]>>face.area){if(id!=faces.size()+1||!(face.area>0))throw std::runtime_error("FACE_ID");for(auto& n:face.nodes){if(n>=N)throw std::runtime_error("FACE_NODE");n=local[n];}area+=face.area;faces.push_back(face);}
  if(faces.size()!=19752||std::abs(area-1.9646292565715553)>1e-12)throw std::runtime_error("ORIGINAL_COVERAGE");
  std::ofstream out(folder+"/initial_check.json");out<<"{\"passed\":true,\"checkedPoints\":"<<N<<",\"initialElevenColumnBitMismatches\":0,\"restartUsed\":false}\n";
}
}
#include "FASRecoveryObserverV248.hpp"
class SpaceOnly248 final:public CIntegration {public:using CIntegration::Space_Integration;};
class Iteration248 final:public CFluidIteration {
 public:using CFluidIteration::CFluidIteration;
 void Iterate(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   observe::phase="native-iterate/preprocessing";CFluidIteration::Iterate(o,in,g,s,n,c,sm,gm,b,z,inst);
   observe::phase="accepted-flow-SA/end-iteration";observe::end_iteration();dev248::finish_cycle();dev248::benefit();if(observe::completed==1){auto prev=dev248::previousFaces;auto changes=dev248::changes;dev248::wall_frame(1);dev248::previousFaces.swap(prev);dev248::changes.swap(changes);}
 }
 bool Monitor(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   const bool native=CFluidIteration::Monitor(o,in,g,s,n,c,sm,gm,b,z,inst);
   if(native&&!dev248::stop)dev248::stopReason="native-convergence-stop-unqualified";
   if(omp_get_wtime()-dev248::started>1680&&!dev248::stop){dev248::stop=true;dev248::stopReason="wallclock-reserve-stop";}
   return native||dev248::stop;
 }
};
class Driver248 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 // 模块二：保存后仅一次完整空间残差，显式记录缓存副作用，核对当前六量/SA涡黏度逐位不变。
 void terminal_residual(){using namespace observe;const auto N=geo->GetnPoint();auto* n=sol[FLOW_SOL]->GetNodes();auto* t=sol[TURB_SOL]->GetNodes();
   std::vector<std::array<double,7>> original(N);
   for(unsigned long i=0;i<N;++i){for(int j=0;j<5;++j)original[i][j]=n->GetSolution(i,j);original[i][5]=t->GetSolution(i,0);original[i][6]=t->GetmuT(i);}
   ledger.close();walls.close();ledger.open(folder+"/terminal_stages.csv");walls.open(folder+"/terminal_wall_stages.csv");ledger<<std::setprecision(17);walls<<std::setprecision(17);
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   dev248::terminal=true;phase="terminal-flow-preprocessing";cfg->SetInnerIter(completed-1);cfg->SetGlobalParam(MAIN_SOLVER::RANS,RUNTIME_FLOW_SYS);
   auto*** num=numerics_container[0][0][0];SpaceOnly248 spatial;
   sol[FLOW_SOL]->Preprocessing(geo,sol,cfg,0,NO_RK_ITER,RUNTIME_FLOW_SYS,false);
   spatial.Space_Integration(geo,sol,num[FLOW_SOL],cfg,0,NO_RK_ITER,RUNTIME_FLOW_SYS);
   phase="terminal-SA-preprocessing";cfg->SetGlobalParam(MAIN_SOLVER::RANS,RUNTIME_TURB_SYS);
   sol[TURB_SOL]->Preprocessing(geo,sol,cfg,0,0,RUNTIME_TURB_SYS,false);
   spatial.Space_Integration(geo,sol,num[TURB_SOL],cfg,0,NO_RK_ITER,RUNTIME_TURB_SYS);
   std::ofstream binary(folder+"/terminal_spatial.bin",std::ios::binary);uint64_t count=N,k=7;binary.write("RAH248RS",8);binary.write(reinterpret_cast<char*>(&count),8);binary.write(reinterpret_cast<char*>(&k),8);
   std::ofstream post(folder+"/terminal_post_state.bin",std::ios::binary);uint64_t pk=8;post.write("RAH248PS",8);post.write(reinterpret_cast<char*>(&count),8);post.write(reinterpret_cast<char*>(&pk),8);
   unsigned long mismatch=0;std::array<long double,6> sum{},sq{},absolute{};
   for(unsigned long i=0;i<N;++i){std::array<double,7> row{};row[0]=geo->nodes->GetGlobalIndex(i);
     std::array<double,8> pr{};pr[0]=row[0];for(int j=0;j<7;++j){double v=j<5?n->GetSolution(i,j):(j==5?t->GetSolution(i,0):t->GetmuT(i));pr[j+1]=v;mismatch+=std::memcmp(&v,&original[i][j],8)!=0;}post.write(reinterpret_cast<char*>(pr.data()),sizeof(pr));
     for(int j=0;j<6;++j){double raw=j<5?sol[FLOW_SOL]->LinSysRes(i,j):sol[TURB_SOL]->LinSysRes(i,0);row[j+1]=raw;double r=raw-(j==4?300000*sol[FLOW_SOL]->LinSysRes(i,0):0);
       if(!std::isfinite(r))throw std::runtime_error("TERMINAL_NONFINITE");sum[j]+=r;sq[j]+=(long double)r*r;absolute[j]+=std::abs(r);}
     binary.write(reinterpret_cast<char*>(row.data()),sizeof(row));
   }binary.close();post.close();
   std::ofstream out(folder+"/terminal_summary.json");out<<std::setprecision(17)<<"{\"completedUpdates\":"<<completed<<",\"additionalTimeUpdates\":0,\"stateAndSAmutBitChanges\":"<<mismatch<<",\"physicalRMS\":[";
   for(int j=0;j<6;++j)out<<(j?",":"")<<std::sqrt(sq[j]/N);out<<"],\"physicalSum\":[";for(int j=0;j<6;++j)out<<(j?",":"")<<sum[j];out<<"],\"physicalSumAbs\":[";for(int j=0;j<6;++j)out<<(j?",":"")<<absolute[j];out<<"],\"fullSteadyBalanceVerified\":false,\"physicalQualification\":false}\n";out.close();
   if(!binary||!post||!out||mismatch)throw std::runtime_error("TERMINAL_STATE_OR_WRITE_FAILURE");
 }
 void launch(const std::string& reference,const std::string& facefile){using namespace observe;cfg=config_container[0];geo=geometry_container[0][0][0];sol=solver_container[0][0][0];dev248::coarse=geometry_container[0][0][1];dev248::cs=solver_container[0][0][1];
   if(cfg->GetnMGLevels()!=1||cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||cfg->GetKind_Upwind_Flow()!=UPWIND::ROE||cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!cfg->GetSAParsedOptions().ft2||cfg->GetRestart()||omp_get_max_threads()!=1||geo->GetnPoint()!=335726||dev248::coarse->GetnPoint()!=44484||cfg->GetnInner_Iter()!=200)throw std::runtime_error("FAS248_IDENTITY");
   dev248::initial(reference,facefile);
   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");steps.open(folder+"/steps.csv");walls.open(folder+"/wall_stages.csv");dev248::trends.open(folder+"/surface_trend.csv");dev248::decisions.open(folder+"/benefit_decisions.csv");
   for(auto* f:{&ledger,&events,&steps,&walls,&dev248::trends,&dev248::decisions})*f<<std::setprecision(17);
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   dev248::open_recovery();events<<"eventID,update,phase,level,localID,globalID,recovered,reason,volume,prolongationMatched,prolongationFactor";for(auto word:{"trial","old","returned","prolongBefore","prolongDirection","prolongAfter"})for(int j=0;j<5;++j)events<<','<<word<<j;events<<",SAcurrent,SAold,mutSA";for(int j=0;j<5;++j)events<<",volumeDeltaPhysical"<<j;events<<'\n';
   steps<<"update,lastFlowRejectionCount,totalObservedRejections,maxWallSpeedComponent,maxWallSA,minT,maxT,maxFlowVsUpdatedSAMuTDifference\n";
   dev248::trends<<"update,pressureAreaL1Change,heatAreaL1Change\n";
   dev248::decisions<<"check,windowEnd,ratio0,ratio1,ratio2,ratio3,ratio4,ratio5,geometricMeanRatio,oldPressureChange,oldHeatChange,newPressureChange,newHeatChange,noSignificantDeterioration,benefit\n";
   delete iteration_container[0][0];iteration_container[0][0]=new Iteration248(cfg);interval=50;enabled=true;
   rah_single_observer::stage=dev248::stage;rah_single_observer::rejection=dev248::reject;rah_fas_observer::transfer=dev248::transfer;
   state_frame(0);dev248::wall_frame(0);dev248::coarse_frame("coarse_initial.bin");StartSolver();if(saved!=completed)state_frame(completed);if(completed%50)dev248::wall_frame(completed);dev248::coarse_frame("coarse_final.bin");
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);if(!forced)throw std::runtime_error("FORCED_FINAL_OUTPUT");
   std::ofstream stopped(folder+"/flow_completion.json");stopped<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"stopReason\":\""<<dev248::stopReason<<"\",\"forcedNativeFinal\":true,\"transfersObserved\":"<<dev248::transfers<<"}\n";stopped.close();
   terminal_residual();dev248::recovery_summary();rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;rah_fas_observer::transfer=nullptr;
   std::ofstream done(folder+"/native_completion.json");done<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"stopReason\":\""<<dev248::stopReason<<"\",\"terminalSpatialEvaluated\":true,\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{if(argc!=5)throw std::runtime_error("CONFIG_OUTPUT_UNIFORM_FACES_REQUIRED");observe::folder=argv[2];std::set_terminate(observe::failed);SU2_MPI::Init(&argc,&argv);Driver248 d(argv[1],1,SU2_MPI::GetComm());d.launch(argv[3],argv[4]);d.Finalize();SU2_MPI::Finalize();return 0;}
 catch(const std::exception& e){std::ofstream f(observe::folder+"/caught_failure.txt");f<<"phase="<<observe::phase<<"\ncompleted="<<observe::completed<<"\nwhat="<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}}
