// 模块职责：在已验证原生入口上添加检查点身份与有界收益停止；不替换任何流动/SA公式。
#define main rah_v235_control_main
#include "singleDomainRunV235.cpp"
#undef main
#include "output/COutput.hpp"
#include <algorithm>

namespace branch236 {
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

// 模块一：在唯一真实任务中，核对原生实际加载的六状态和全局坐标；不重新初始化或再评估残差。
void initial(const std::string& checkpoint,const std::string& facefile){
 const auto N=geo->GetnPoint();
 if(N!=335726||geo->GetnPointDomain()!=N||!cfg->GetRestart())throw std::runtime_error("BRANCH_MESH_RESTART_IDENTITY");
 local.assign(N,~0ul);
 for(unsigned long i=0;i<N;++i){auto id=geo->nodes->GetGlobalIndex(i);if(id>=N||local[id]!=~0ul)throw std::runtime_error("GLOBAL_ID_LAYOUT");local[id]=i;}
 std::ifstream f(checkpoint);std::string line;std::getline(f,line);
 const std::string prefix="\"PointID\",\"x\",\"y\",\"z\",\"Density\",\"Momentum_x\",\"Momentum_y\",\"Momentum_z\",\"Energy\",\"Nu_Tilde\",";
 if(line.rfind(prefix,0)!=0)throw std::runtime_error("CHECKPOINT_COLUMNS");
 std::vector<bool> seen(N,false);unsigned long count=0,mismatch=0,coordMismatch=0;double cflMin=1e100,cflMax=0;
 while(count<N&&std::getline(f,line)){
   std::istringstream row(line);std::string token;std::array<double,10> a{};
   for(auto& v:a){if(!std::getline(row,token,','))throw std::runtime_error("CHECKPOINT_TRUNCATED");v=std::stod(token);if(!std::isfinite(v))throw std::runtime_error("CHECKPOINT_NONFINITE");}
   const auto id=static_cast<unsigned long>(a[0]);if(id>=N||a[0]!=id||seen[id])throw std::runtime_error("CHECKPOINT_ID");seen[id]=true;auto i=local[id];
   for(int d=0;d<3;++d){double value=geo->nodes->GetCoord(i,d);coordMismatch+=std::memcmp(&value,&a[d+1],8)!=0;}
   for(int j=0;j<6;++j){double value=j<5?sol[FLOW_SOL]->GetNodes()->GetSolution(i,j):sol[TURB_SOL]->GetNodes()->GetSolution(i,0);mismatch+=std::memcmp(&value,&a[j+4],8)!=0;}
   const double cfl=sol[FLOW_SOL]->GetNodes()->GetLocalCFL(i);cflMin=std::min(cflMin,cfl);cflMax=std::max(cflMax,cfl);++count;
 }
 std::ofstream report(folder+"/initial_check.json");report<<std::setprecision(17)<<"{\"points\":"<<count<<",\"stateBitMismatches\":"<<mismatch<<",\"coordinateBitMismatches\":"<<coordMismatch
   <<",\"flowCFLMin\":"<<cflMin<<",\"flowCFLMax\":"<<cflMax<<",\"exactHistoricalContinuation\":false,\"initialRHSUnknown\":true}\n";report.close();
 if(count!=N||mismatch||coordMismatch||cflMin!=cfg->GetCFL(0)||cflMax!=cfg->GetCFL(0))throw std::runtime_error("RESTART_INITIAL_STATE_MISMATCH");
 std::ifstream ff(facefile);Face face;unsigned long id=0;double total=0;
 while(ff>>id>>face.nodes[0]>>face.nodes[1]>>face.nodes[2]>>face.area){
   if(id!=faces.size()+1||!(face.area>0))throw std::runtime_error("ORIGINAL_FACE_ID");
   for(auto& n:face.nodes){if(n>=N)throw std::runtime_error("ORIGINAL_FACE_NODE");n=local[n];}total+=face.area;faces.push_back(face);
 }
 if(faces.size()!=19752||std::abs(total-1.9646292565715553)>1e-12)throw std::runtime_error("ORIGINAL_FACE_COVERAGE");
}

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

void stage236(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short system){
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

class BranchIterationV236 final:public CFluidIteration {
 public:using CFluidIteration::CFluidIteration;
 void Iterate(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   observe::phase="native-iterate/preprocessing";CFluidIteration::Iterate(o,in,g,s,n,c,sm,gm,b,z,inst);
   observe::phase="accepted-flow-SA/end-iteration";observe::end_iteration();branch236::benefit();
 }
 bool Monitor(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst) override {
   const bool native=CFluidIteration::Monitor(o,in,g,s,n,c,sm,gm,b,z,inst);
   if(native&&!branch236::stop)branch236::stopReason="native-convergence-stop-unqualified";
   return native||branch236::stop;
 }
};

