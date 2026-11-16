// 模块职责：真实NK有限发展；只记录基准/Jv/实际更新身份，复用完整原生flow-SA和既有停止规则。
#define main rah_unused_single_main_v250
#include "singleDomainRunV235.cpp"
#undef main
#include "integration/CNewtonIntegration.hpp"
#include "NewtonHooksV250.hpp"
namespace dev250 {
using namespace observe;
struct Face {std::array<unsigned long,3> nodes;double area;};
std::vector<Face> faces;std::vector<unsigned long> local;
std::vector<std::array<double,6>> logs;std::vector<std::array<double,2>> previousFaces,changes;
std::array<double,6> currentRMS{};bool stop=false,terminal=false;
std::string stopReason="maximum-updates",role="initialization",eval="none";
std::ofstream trends,decisions,contexts,phases,domain;unsigned long stageRows=0,products=0;
const rah_newton_observer::Direction* direction=nullptr;double magnitude=0,fdstep=0,started=omp_get_wtime();
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

// 模块二：只读阶段/方向观察；阶段身份独立于通用空间观察器的标签。
void event(const char* raw,CGeometry* g,CSolver** s,CConfig* c,const rah_newton_observer::Direction* d,double mag,double step){
  const std::string tag=raw;
  if(tag=="jv_begin"){role="jv";direction=d;magnitude=mag;fdstep=step;++products;}
  if(tag=="default_begin"){role="baseline";eval="DEFAULT";direction=nullptr;}
  if(tag=="explicit_begin")eval="EXPLICIT";
  if(tag=="update_begin"){role="actual_update";direction=d;magnitude=mag;fdstep=step;}
  if(tag=="accepted_preprocessing_begin")role="accepted_preprocessing";
  phase="NK/"+role+"/"+tag;
  std::array<double,5> maximum{};double physicalEnergyMax=0;
  if(d)for(unsigned long i=0;i<g->GetnPointDomain();++i){
    for(int j=0;j<5;++j)maximum[j]=std::max(maximum[j],std::abs((*d)(i,j)));
    physicalEnergyMax=std::max(physicalEnergyMax,std::abs((*d)(i,4)-300000*(*d)(i,0)));
  }
  phases<<c->GetInnerIter()+1<<','<<tag<<','<<role<<','<<eval<<','<<products<<','<<mag<<','<<step;
  for(auto x:maximum)phases<<','<<x;phases<<','<<physicalEnergyMax<<'\n';phases.flush();
  if(!phases)throw std::runtime_error("NK_PHASE_WRITE");
  if(tag=="jv_end"){role="between_products";direction=nullptr;}
  if(tag=="flow_complete"){role="SA";direction=nullptr;}
}
void stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short system){
  const auto which=c->GetContainerPosition(system);const auto name=std::string(label);
  contexts<<++stageRows<<','<<c->GetInnerIter()+1<<','<<which<<','<<label<<','<<(terminal?"terminal":role)<<','<<eval<<','<<products<<'\n';contexts.flush();
  if(!contexts)throw std::runtime_error("NK_CONTEXT_WRITE");
  observe::stage(label,g,s,c,system);phase="NK/"+(terminal?std::string("terminal"):role)+"/"+label;
  if(terminal)return;
  if((which==FLOW_SOL&&role=="baseline"&&eval=="DEFAULT"&&name=="after_strong")||(which==TURB_SOL&&name=="after_time")){
    const auto N=g->GetnPointDomain();for(unsigned j=0;j<(which==FLOW_SOL?5:1);++j){long double sum=0;
      for(unsigned long i=0;i<N;++i){double r=s[which]->LinSysRes(i,j);if(which==FLOW_SOL&&j==4)r-=300000*s[which]->LinSysRes(i,0);
        if(!std::isfinite(r))throw std::runtime_error("NK_BASE_RHS_NONFINITE");sum+=(long double)r*r;}
      currentRMS[which==FLOW_SOL?j:5]=std::sqrt(sum/N);
    }
  }
}
void reject(CNSVariable* n,unsigned long i,const double* trial,const double* old,bool recovered,const char* reason){
  phase="NK/"+role+"/"+eval;observe::rejection(n,i,trial,old,recovered,reason);
  domain<<cfg->GetInnerIter()+1<<','<<i<<','<<geo->nodes->GetGlobalIndex(i)<<','<<role<<','<<eval<<','<<products<<','<<magnitude<<','<<fdstep<<','<<n->GetUnderRelaxation(i)<<','<<(direction?1:0);
  for(int j=0;j<5;++j)domain<<','<<(direction?(*direction)(i,j):std::numeric_limits<double>::quiet_NaN());domain<<'\n';domain.flush();
  if(!domain)throw std::runtime_error("NK_DOMAIN_WRITE");
  if(role=="jv")throw std::runtime_error("NK_JV_DOMAIN_RECOVERY_REJECTED");
  if(!recovered)throw std::runtime_error("NK_ILLEGAL_OLD_STATE");
}
void read_faces(const std::string& file){
  local.assign(geo->GetnPoint(),~0ul);for(unsigned long i=0;i<geo->GetnPoint();++i)local.at(geo->nodes->GetGlobalIndex(i))=i;
  std::ifstream ff(file);Face face;unsigned long id=0;double area=0;
  while(ff>>id>>face.nodes[0]>>face.nodes[1]>>face.nodes[2]>>face.area){if(id!=faces.size()+1||!(face.area>0))throw std::runtime_error("FACE_ID");
    for(auto& n:face.nodes)n=local.at(n);area+=face.area;faces.push_back(face);}
  if(faces.size()!=19752||std::abs(area-1.9646292565715553)>1e-12)throw std::runtime_error("ORIGINAL_COVERAGE");
}
}
class Iteration250 final:public CFluidIteration {
 public:using CFluidIteration::CFluidIteration;
 void Iterate(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   observe::phase="native-iterate/preprocessing";CFluidIteration::Iterate(o,in,g,s,n,c,sm,gm,b,z,inst);
   observe::phase="accepted-flow-SA/end-iteration";observe::end_iteration();dev250::benefit();
   if(observe::completed==1){auto prev=dev250::previousFaces;auto changes=dev250::changes;dev250::wall_frame(1);dev250::previousFaces.swap(prev);dev250::changes.swap(changes);}
 }
 bool Monitor(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   const bool native=CFluidIteration::Monitor(o,in,g,s,n,c,sm,gm,b,z,inst);
   if(native&&!dev250::stop)dev250::stopReason="native-convergence-stop-unqualified";
   if(omp_get_wtime()-dev250::started>1680&&!dev250::stop){dev250::stop=true;dev250::stopReason="wallclock-reserve-stop";}
   return native||dev250::stop;
 }
};
class SpaceOnly250 final:public CIntegration {public:using CIntegration::Space_Integration;};
class Driver250 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void terminal_residual(){using namespace observe;const auto N=geo->GetnPoint();auto* n=sol[FLOW_SOL]->GetNodes();auto* t=sol[TURB_SOL]->GetNodes();
   std::vector<std::array<double,7>> original(N);
   for(unsigned long i=0;i<N;++i){for(int j=0;j<5;++j)original[i][j]=n->GetSolution(i,j);original[i][5]=t->GetSolution(i,0);original[i][6]=t->GetmuT(i);}
   ledger.close();walls.close();ledger.open(folder+"/terminal_stages.csv");walls.open(folder+"/terminal_wall_stages.csv");ledger<<std::setprecision(17);walls<<std::setprecision(17);
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   dev250::terminal=true;phase="terminal-flow-preprocessing";cfg->SetInnerIter(completed-1);cfg->SetGlobalParam(MAIN_SOLVER::RANS,RUNTIME_FLOW_SYS);
   auto*** num=numerics_container[0][0][0];SpaceOnly250 spatial;
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
 void launch(const std::string& reference,const std::string& facefile){using namespace observe;cfg=config_container[0];geo=geometry_container[0][0][0];sol=solver_container[0][0][0];
   auto* engine=dynamic_cast<CNewtonIntegration*>(integration_container[0][0][FLOW_SOL]);auto ip=cfg->GetNewtonKrylovIntParam();auto dp=cfg->GetNewtonKrylovDblParam();
   if(!engine||!cfg->GetNewtonKrylov()||cfg->GetnMGLevels()!=0||cfg->GetnInner_Iter()!=200||ip[0]!=0||dp[0]!=0||dp[3]!=1e-4||dp[4]!=1||
      cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||cfg->GetKind_Upwind_Flow()!=UPWIND::ROE||cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||
      !cfg->GetSAParsedOptions().ft2||!cfg->GetRestart()||geo->GetnPoint()!=335726||omp_get_max_threads()!=1)throw std::runtime_error("NK250_NATIVE_IDENTITY");
   dev250::initial(reference);dev250::read_faces(facefile);
   ledger.open(folder+"/stages.csv");events.open(folder+"/rejections.csv");steps.open(folder+"/steps.csv");walls.open(folder+"/wall_stages.csv");
   dev250::trends.open(folder+"/surface_trend.csv");dev250::decisions.open(folder+"/benefit_decisions.csv");
   dev250::contexts.open(folder+"/evaluation_context.csv");dev250::phases.open(folder+"/newton_phases.csv");dev250::domain.open(folder+"/domain_context.csv");
   for(auto* f:{&ledger,&events,&steps,&walls,&dev250::trends,&dev250::decisions,&dev250::contexts,&dev250::phases,&dev250::domain})*f<<std::setprecision(17);
   walls<<"update,stage,pNoutX,pNoutY,pNoutZ,tauNoutX,tauNoutY,tauNoutZ,wallRowX,wallRowY,wallRowZ,laminarHeatIntoWallW,effectiveHeatIntoWallW,maxWallSA\n";
   ledger<<"update,system,stage,nodes";for(auto word:{"sum","rms"})for(int j=0;j<6;++j)ledger<<','<<word<<j;ledger<<",minRelax,maxCFL,maxScaledUpdate\n";
   events<<"update,phase,localID,globalID,recovered,reason";for(auto word:{"trial","old","returned"})for(int j=0;j<5;++j)events<<','<<word<<j;events<<",SAcurrent,SAold\n";
   steps<<"update,lastFlowRejectionCount,totalObservedRejections,maxWallSpeedComponent,maxWallSA,minT,maxT,maxFlowVsUpdatedSAMuTDifference\n";
   dev250::trends<<"update,pressureAreaL1Change,heatAreaL1Change\n";
   dev250::decisions<<"check,windowEnd,ratio0,ratio1,ratio2,ratio3,ratio4,ratio5,geometricMeanRatio,oldPressureChange,oldHeatChange,newPressureChange,newHeatChange,noSignificantDeterioration,benefit\n";
   dev250::contexts<<"row,update,system,stage,role,evaluation,productCount\n";
   dev250::phases<<"update,event,role,evaluation,productCount,magnitude,finiteDiffStep,maxDir0,maxDir1,maxDir2,maxDir3,maxDir4,maxPhysicalEnergyDirection\n";
   dev250::domain<<"update,localID,globalID,role,evaluation,productCount,magnitude,finiteDiffStep,localRelaxation,directionAvailable,dir0,dir1,dir2,dir3,dir4\n";
   enabled=true;interval=50;delete iteration_container[0][0];iteration_container[0][0]=new Iteration250(cfg);
   rah_single_observer::stage=dev250::stage;rah_single_observer::rejection=dev250::reject;rah_newton_observer::event=dev250::event;
   state_frame(0);dev250::wall_frame(0);StartSolver();if(saved!=completed)state_frame(completed);if(completed%50)dev250::wall_frame(completed);
   const bool forced=output_container[0]->SetResultFiles(geo,cfg,sol,cfg->GetInnerIter(),true);if(!forced)throw std::runtime_error("FORCED_FINAL_OUTPUT");
   std::ofstream done(folder+"/flow_completion.json");done<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"stopReason\":\""<<dev250::stopReason<<"\",\"forcedNativeFinal\":true}\n";done.close();
   terminal_residual();rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;rah_newton_observer::event=nullptr;
   std::ofstream complete(folder+"/native_completion.json");complete<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"stopReason\":\""<<dev250::stopReason<<"\",\"actualNewtonFactory\":true,\"matrixFreeProducts\":"<<dev250::products<<",\"terminalSpatialEvaluated\":true,\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{if(argc!=5)throw std::runtime_error("CONFIG_FOLDER_REFERENCE_FACES_REQUIRED");observe::folder=argv[2];std::set_terminate(observe::failed);SU2_MPI::Init(&argc,&argv);Driver250 d(argv[1],1,SU2_MPI::GetComm());d.launch(argv[3],argv[4]);d.Finalize();SU2_MPI::Finalize();return 0;}
catch(const std::exception& e){std::ofstream f(observe::folder+"/caught_failure.txt");f<<"phase="<<observe::phase<<"\ncompleted="<<observe::completed<<"\nwhat="<<e.what()<<'\n';std::cerr<<e.what()<<'\n';return 1;}}
