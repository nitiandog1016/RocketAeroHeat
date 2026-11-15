// 模块职责：完整原域HLLC均匀初场、固定收益停止与实际接受/失败观察；成熟推进循环保持。
#include "commonBranchV239.hpp"
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

