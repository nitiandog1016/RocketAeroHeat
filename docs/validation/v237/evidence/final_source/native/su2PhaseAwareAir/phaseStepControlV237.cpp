// 模块职责：回放真实越域事件，调用同一生产守卫与真实CNSVariable；不构造体网格或推进流场。
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include "SingleDomainHooksV235.hpp"
#include "PhaseStepAcceptanceV237.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <limits>

int main(int argc,char** argv){try{
 if(argc!=4)throw std::runtime_error("CONFIG_INPUT_OUTPUT_REQUIRED");SU2_MPI::Init(&argc,&argv);
 CConfig cfg(argv[1],SU2_COMPONENT::SU2_CFD,false);CSU2PhaseAwareAir gas;
 const auto q=rah_phase_air::PT(1000.,300.);double velocity[3]={0,0,0};CNSVariable node(q.rho,velocity,q.e,1,3,5,&cfg);
 using namespace rah_phase_step_v237;
 std::ifstream input(argv[2]);std::string line;std::getline(input,line);
 std::ofstream out(std::string(argv[3])+"/accepted.csv"),reject(std::string(argv[3])+"/rejected.csv");
 out<<std::setprecision(17)<<"event,update,globalID,alpha,reductions,shiftedAlpha,shiftedMaxError,legalAlpha,zeroAlpha,p,T";
 for(auto prefix:{"old","delta","accepted"})for(int j=0;j<5;++j)out<<','<<prefix<<j;out<<"\n";reject<<"case,error\n";
 unsigned count=0,limitedCount=0;State firstOld{},firstTrial{};
 while(std::getline(input,line)){
   std::istringstream row(line);std::vector<std::string> cells;std::string cell;while(std::getline(row,cell,','))cells.push_back(cell);
   if(cells.size()!=23)throw std::runtime_error("REAL_EVENT_COLUMNS");State old{},trial{},delta{};
   for(int j=0;j<5;++j){trial[j]=std::stod(cells[j+6]);old[j]=std::stod(cells[j+11]);delta[j]=trial[j]-old[j];}
   if(!count){firstOld=old;firstTrial=trial;}
   if(admissible(trial)||!admissible(old))throw std::runtime_error("REAL_EVENT_DOMAIN");
   const auto accepted=select(old,delta,1.);if(accepted.reductions)++limitedCount;
   node.SetSolution_Old(0,old.data());node.SetSolution(0,accepted.accepted.data());
   if(!node.SetPrimVar(0,0.,0.,&gas))throw std::runtime_error("ACCEPTED_STILL_REVERTED");
   for(int j=0;j<5;++j)if(std::memcmp(&accepted.accepted[j],node.GetSolution(0)+j,8))throw std::runtime_error("NATIVE_ACCEPTED_STATE_CHANGED");
   State shiftedOld=old,shiftedDelta=delta;shiftedOld[4]-=300000*old[0];shiftedDelta[4]-=300000*delta[0];
   const auto shifted=select(shiftedOld,shiftedDelta,1.,0.);if(shifted.alpha!=accepted.alpha)throw std::runtime_error("GUARD_REFERENCE_ALPHA");
   double maxError=0;for(int j=0;j<5;++j){double converted=accepted.accepted[j]-(j==4?300000*accepted.accepted[0]:0.);maxError=std::max(maxError,std::abs(shifted.accepted[j]-converted)/std::max(1.,std::abs(converted)));}
   State legal{},zero{};for(int j=0;j<5;++j)legal[j]=accepted.accepted[j]-old[j];const auto legalResult=select(old,legal,1.),zeroResult=select(old,zero,1.);
   if(legalResult.alpha!=1.||zeroResult.alpha!=1.||std::memcmp(zeroResult.accepted.data(),old.data(),sizeof(State)))throw std::runtime_error("LEGAL_UPDATE_CHANGED");
   out<<++count<<','<<cells[0]<<','<<cells[3]<<','<<accepted.alpha<<','<<accepted.reductions<<','<<shifted.alpha<<','<<maxError<<','<<legalResult.alpha<<','<<zeroResult.alpha<<','<<node.GetPressure(0)<<','<<node.GetTemperature(0);
   for(const State* values:std::array<const State*,3>{&old,&delta,&accepted.accepted})for(double x:*values)out<<','<<x;out<<'\n';
 }
 if(count!=466||limitedCount!=466)throw std::runtime_error("REAL_EVENT_COUNT");
 unsigned rejected=0;auto expect=[&](const char* name,const State& old,const State& d,double alpha){
   try{select(old,d,alpha);}catch(const std::exception& e){reject<<name<<','<<e.what()<<'\n';++rejected;return;}throw std::runtime_error("INVALID_CONTROL_ACCEPTED");};
 State nonzero{};nonzero[4]=1;State nonfinite=nonzero;nonfinite[4]=std::numeric_limits<double>::quiet_NaN();
 expect("invalid_old",firstTrial,nonzero,1);expect("nonfinite_delta",firstOld,nonfinite,1);
 expect("negative_upper",firstOld,nonzero,-1);expect("above_one_upper",firstOld,nonzero,1.01);expect("nan_upper",firstOld,nonzero,std::numeric_limits<double>::quiet_NaN());
 expect("zero_progress",firstOld,nonzero,0);
 const auto cold=rah_phase_air::PT(1.,60.0001);State boundary{cold.rho,0,0,0,cold.rho*cold.e},huge{};huge[4]=-1e8;expect("no_admissible_progress",boundary,huge,1.);
 if(domain_rejection("RAH_PHASE_INVERSE_CONVERGENCE")||domain_rejection("RAH_PHASE_NUMERICAL_RANGE")||domain_rejection("unknown"))throw std::runtime_error("ERROR_CLASSIFICATION");
 std::ofstream report(std::string(argv[3])+"/report.json");report<<"{\"passed\":true,\"realEvents\":"<<count<<",\"limitedEvents\":"<<limitedCount<<",\"rejections\":"<<rejected<<",\"nativeCNSAcceptedWithoutFallback\":true,\"flowUpdates\":0,\"productionCallsiteExecuted\":false,\"physicalQualification\":false}\n";
 if(!out||!reject||!report)throw std::runtime_error("CONTROL_OUTPUT");SU2_MPI::Finalize();return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
