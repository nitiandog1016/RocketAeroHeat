// 模块职责：被动记录原生合法恢复及真实延拓方向；不修改任一求解状态或接受公式。
#pragma once
namespace dev248 {
struct RecoveryEvent {unsigned long id,update,point;unsigned layer;std::array<double,5> returned;bool checked=false;};
std::vector<RecoveryEvent> cycle_events;
std::vector<std::array<double,5>> prolong_before,prolong_direction,prolong_after;
unsigned long prolong_step=0,event_count=0,return_checks=0,cycle_checks=0;
double prolong_factor=0,max_prolong_error=0;
std::array<long double,5> disturbance{};
std::ofstream recovery_checks;

void record_recovery_state(const RecoveryEvent& e,const char* label){
  auto* g=e.layer==0?geo:coarse;auto** s=e.layer==0?sol:cs;auto* n=s[FLOW_SOL]->GetNodes();auto* t=s[TURB_SOL]->GetNodes();
  recovery_checks<<e.id<<','<<e.update<<','<<e.layer<<','<<g->nodes->GetGlobalIndex(e.point)<<','<<label;
  for(int j=0;j<5;++j)recovery_checks<<','<<n->GetSolution(e.point,j);
  recovery_checks<<','<<t->GetSolution(e.point,0)<<','<<t->GetSolution_Old(e.point,0)<<','<<t->GetmuT(e.point);
  for(int j:{0,4,6,9,10,11})recovery_checks<<','<<n->GetPrimitive(e.point,j);
  recovery_checks<<'\n';recovery_checks.flush();if(!recovery_checks)throw std::runtime_error("RECOVERY_CHECK_WRITE");
}
// 首个原生后续阶段核对returned仍在；完整循环终点另记，不能把两个阶段混成同一状态。
void verify_returns(CGeometry* g,const char* label){
  const unsigned layer=g==geo?0:1;auto** s=layer==0?sol:cs;auto* n=s[FLOW_SOL]->GetNodes();
  for(auto& e:cycle_events)if(e.layer==layer&&!e.checked){
    for(int j=0;j<5;++j){double u=n->GetSolution(e.point,j);if(std::memcmp(&u,&e.returned[j],8))throw std::runtime_error("RECOVERY_RETURN_CHANGED_BEFORE_OBSERVATION");}
    record_recovery_state(e,label);e.checked=true;++return_checks;
  }
}
void finish_cycle(){
  for(const auto& e:cycle_events){if(!e.checked||e.update!=completed)throw std::runtime_error("RECOVERY_STAGE_NOT_COMPLETED");record_recovery_state(e,"complete_flow_SA");++cycle_checks;}
  cycle_events.clear();
}
void stage(const char* label,CGeometry* g,CSolver** s,CConfig* c,unsigned short sys){
  verify_returns(g,label);observe::phase=std::string(g==geo?"fine/":"coarse/")+label;
  if(g!=geo)return;
  if(terminal)observe::stage(label,g,s,c,sys);else fine_stage(label,g,s,c,sys);
}
void transfer(const char* label,CGeometry* f,CSolver**,CGeometry* c,CSolver**,CConfig*){
  ++transfers;verify_returns(f,label);verify_returns(c,label);observe::phase=std::string("transfer/")+label;
  const auto N=geo->GetnPointDomain();auto* n=sol[FLOW_SOL]->GetNodes();
  if(std::string(label)=="prolong_in"){
    prolong_before.resize(N);prolong_direction.resize(N);prolong_after.resize(N);prolong_step=cfg->GetInnerIter()+1;prolong_factor=cfg->GetDamp_Correc_Prolong();
    if(prolong_factor!=.5)throw std::runtime_error("PROLONG_FACTOR_CHANGED");
    for(unsigned long i=0;i<N;++i)for(int j=0;j<5;++j){prolong_before[i][j]=n->GetSolution(i,j);prolong_direction[i][j]=sol[FLOW_SOL]->LinSysRes(i,j);if(!std::isfinite(prolong_before[i][j]+prolong_direction[i][j]))throw std::runtime_error("PROLONG_NONFINITE");}
  }
  if(std::string(label)=="prolong_out"){
    if(prolong_step!=cfg->GetInnerIter()+1)throw std::runtime_error("PROLONG_STAGE_ID");
    for(unsigned long i=0;i<N;++i)for(int j=0;j<5;++j){const double actual=n->GetSolution(i,j),expected=prolong_before[i][j]+prolong_factor*prolong_direction[i][j];prolong_after[i][j]=actual;
      const double error=std::abs(actual-expected)/std::max(1.,std::abs(expected));max_prolong_error=std::max(max_prolong_error,error);if(!std::isfinite(actual)||error>2e-14)throw std::runtime_error("PROLONG_UPDATE_MISMATCH");}
  }
}
void reject(CNSVariable* n,unsigned long i,const double* trial,const double* old,bool recovered,const char* reason){
  const unsigned layer=n==sol[FLOW_SOL]->GetNodes()?0:(n==cs[FLOW_SOL]->GetNodes()?1:99);if(layer>1)throw std::runtime_error("RECOVERY_LAYER_ID");
  auto* g=layer==0?geo:coarse;auto** s=layer==0?sol:cs;auto* t=s[TURB_SOL]->GetNodes();const auto step=cfg->GetInnerIter()+1;
  if(++event_count>100000)throw std::runtime_error("RECOVERY_OUTPUT_LIMIT");++rejections;
  bool matched=layer==0&&step==prolong_step&&phase=="transfer/prolong_out";
  if(matched)for(int j=0;j<5;++j)matched=matched&&!std::memcmp(&trial[j],&prolong_after[i][j],8);
  RecoveryEvent e{event_count,step,i,layer,{}};const double volume=g->nodes->GetVolume(i);
  events<<e.id<<','<<step<<','<<phase<<','<<layer<<','<<i<<','<<g->nodes->GetGlobalIndex(i)<<','<<recovered<<','<<reason<<','<<volume<<','<<matched<<','<<prolong_factor;
  for(auto p:{trial,old})for(int j=0;j<5;++j)events<<','<<p[j];
  for(int j=0;j<5;++j){e.returned[j]=n->GetSolution(i,j);events<<','<<e.returned[j];}
  for(auto* a:{&prolong_before,&prolong_direction,&prolong_after})for(int j=0;j<5;++j)events<<','<<(matched?(*a)[i][j]:std::numeric_limits<double>::quiet_NaN());
  events<<','<<t->GetSolution(i,0)<<','<<t->GetSolution_Old(i,0)<<','<<t->GetmuT(i);
  for(int j=0;j<5;++j){const double delta=((e.returned[j]-trial[j])-(j==4?300000*(e.returned[0]-trial[0]):0))*volume;events<<','<<delta;disturbance[j]+=delta;}
  events<<'\n';events.flush();if(!events)throw std::runtime_error("RECOVERY_EVENT_WRITE");
  if(!recovered)throw std::runtime_error("FAS_NATIVE_OLD_STATE_REJECTED");
  for(int j=0;j<5;++j)if(!std::isfinite(e.returned[j])||std::memcmp(&e.returned[j],&old[j],8))throw std::runtime_error("NATIVE_RECOVERY_IDENTITY");
  cycle_events.push_back(e); // 只记录，继续原生SetPrimVar的输运刷新与正常返回。
}
void open_recovery(){
  recovery_checks.open(folder+"/recovery_checks.csv");recovery_checks<<std::setprecision(17);
  recovery_checks<<"eventID,update,level,globalID,checkpoint";for(int j=0;j<5;++j)recovery_checks<<",U"<<j;
  recovery_checks<<",SAcurrent,SAold,mutSA,T,p,Hinternal,mutFlow,k,Cp\n";
}
void recovery_summary(){
  std::ofstream out(folder+"/recovery_summary.json");out<<std::setprecision(17)<<"{\"events\":"<<event_count<<",\"firstReturnChecks\":"<<return_checks<<",\"completeCycleChecks\":"<<cycle_checks<<",\"pendingEvents\":"<<cycle_events.size()<<",\"maximumProlongationScaledDifference\":"<<max_prolong_error<<",\"physicalVolumeStateDisturbance\":[";
  for(int j=0;j<5;++j)out<<(j?",":"")<<disturbance[j];out<<"],\"stateDisturbanceIsBoundaryFluxBalance\":false}\n";
}
}
