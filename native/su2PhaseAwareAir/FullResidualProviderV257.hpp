#pragma once
// 模块职责：原成熟flow-SA完整空间残差F(x)接口；不调用时间推进/旧Newton，不造人工出口。
#define main rah_unused_main_257_provider
#include "singleDomainRunV235.cpp"
#undef main
#include <algorithm>

class Space256 final : public CIntegration {public:using CIntegration::Space_Integration;};
class FullResidualProviderV257 : public CSinglezoneDriver {
protected:
  CGeometry* g=nullptr;CSolver** s=nullptr;CConfig* c=nullptr;CNumerics*** num=nullptr;
  std::vector<unsigned char> wall;std::vector<std::pair<unsigned long,unsigned>> entries;
  unsigned long calls=0,physicalCalls=0,rejected=0,probeNode=~0ul;bool nativeRecovery=false;
  std::string destination,lastReason;std::ofstream log;
  static FullResidualProviderV257* active;
  static void reject_hook(CNSVariable*,unsigned long,const double*,const double*,bool,const char* why){
    if(active){active->nativeRecovery=true;active->lastReason=why;}
  }
  void dump(const std::string& name,const std::vector<double>& values){
    std::ofstream f(destination+"/"+name,std::ios::binary);uint64_t n=values.size();f.write("RAH256V1",8);f.write(reinterpret_cast<char*>(&n),8);
    f.write(reinterpret_cast<const char*>(values.data()),values.size()*8);f.close();if(!f)throw std::runtime_error("VECTOR_WRITE");
  }
  std::vector<double> pack()const{
    std::vector<double> x;x.reserve(entries.size());
    for(auto e:entries){const auto i=e.first;const auto j=e.second;
      double value=j<5?s[FLOW_SOL]->GetNodes()->GetSolution(i,j):s[TURB_SOL]->GetNodes()->GetSolution(i,0);
      if(j==4)value-=300000*s[FLOW_SOL]->GetNodes()->GetSolution(i,0);x.push_back(value);
    }return x;
  }
  std::vector<std::array<double,6>> decode(const std::vector<double>& x)const{
    if(x.size()!=entries.size())throw std::invalid_argument("DOF_SHAPE");
    std::vector<std::array<double,6>> u(g->GetnPoint());
    for(size_t k=0;k<x.size();++k){if(!std::isfinite(x[k]))throw std::domain_error("NONFINITE_INPUT");u[entries[k].first][entries[k].second]=x[k];}
    for(auto& a:u){
      if(!(a[0]>0))throw std::domain_error("NONPOSITIVE_DENSITY");
      a[4]+=300000*a[0];const double kinetic=.5*(a[1]*a[1]+a[2]*a[2]+a[3]*a[3])/(a[0]*a[0]);
      rah_phase_air::state("rhoe",a[0],a[4]/a[0]-kinetic);
    }return u;
  }
  void commit(const std::vector<std::array<double,6>>& u){
    auto* f=s[FLOW_SOL]->GetNodes();auto* t=s[TURB_SOL]->GetNodes();
    for(unsigned long i=0;i<u.size();++i){for(unsigned j=0;j<5;++j){f->SetSolution(i,j,u[i][j]);f->SetSolution_Old(i,j,u[i][j]);}
      t->SetSolution(i,0,u[i][5]);t->SetSolution_Old(i,0,u[i][5]);}
  }
  void assert_unmodified(const std::vector<std::array<double,6>>& u)const{
    for(unsigned long i=0;i<u.size();++i)for(unsigned j=0;j<6;++j){const double actual=j<5?s[FLOW_SOL]->GetNodes()->GetSolution(i,j):s[TURB_SOL]->GetNodes()->GetSolution(i,0);
      if(std::memcmp(&actual,&u[i][j],8))throw std::runtime_error("NATIVE_CHANGED_INPUT_STATE");}
  }
  void state(const std::string& name){
    std::ofstream f(destination+"/"+name,std::ios::binary);uint64_t n=g->GetnPoint(),cols=19;f.write("RAH256S1",8);f.write(reinterpret_cast<char*>(&n),8);f.write(reinterpret_cast<char*>(&cols),8);
    auto* flow=s[FLOW_SOL]->GetNodes();auto* sa=s[TURB_SOL]->GetNodes();
    for(unsigned long i=0;i<n;++i){std::array<double,19>a{};a[0]=g->nodes->GetGlobalIndex(i);for(unsigned j=0;j<3;++j)a[1+j]=g->nodes->GetCoord(i,j);
      a[4]=g->nodes->GetVolume(i);for(unsigned j=0;j<5;++j)a[5+j]=flow->GetSolution(i,j);a[9]-=300000*a[5];
      a[10]=sa->GetSolution(i,0);a[11]=flow->GetTemperature(i);a[12]=flow->GetPressure(i);a[13]=flow->GetThermalConductivity(i);a[14]=flow->GetLaminarViscosity(i);
      a[15]=sa->GetmuT(i);a[16]=flow->GetEddyViscosity(i);a[17]=wall[i];a[18]=g->nodes->GetWall_Distance(i);f.write(reinterpret_cast<const char*>(a.data()),sizeof(a));
    }f.close();if(!f)throw std::runtime_error("STATE_WRITE");
  }
  void full_residual(const std::string& label){
    std::ofstream f(destination+"/"+label+"_full_residual.bin",std::ios::binary);uint64_t n=g->GetnPoint(),cols=7;f.write("RAH256R1",8);f.write(reinterpret_cast<char*>(&n),8);f.write(reinterpret_cast<char*>(&cols),8);
    for(unsigned long i=0;i<n;++i){std::array<double,7>a{};a[0]=g->nodes->GetGlobalIndex(i);for(unsigned j=0;j<6;++j)a[j+1]=j<5?s[FLOW_SOL]->LinSysRes(i,j):s[TURB_SOL]->LinSysRes(i,0);
      a[5]-=300000*a[1];f.write(reinterpret_cast<const char*>(a.data()),sizeof(a));}
    f.close();if(!f)throw std::runtime_error("RESIDUAL_WRITE");
  }
public:
  using CSinglezoneDriver::CSinglezoneDriver;
  // 完整输入先核验，正返回值为可恢复拒绝且不发布残差；负返回值为内部错误。
  int evaluate(const std::vector<double>& x,std::vector<double>& result,const std::string& label){
    ++calls;std::vector<std::array<double,6>> decoded;
    try{decoded=decode(x);}catch(const std::domain_error& e){++rejected;lastReason=e.what();log<<calls<<','<<label<<",1,"<<lastReason<<",0\n";log.flush();return 1;}
    commit(decoded);nativeRecovery=false;active=this;rah_single_observer::rejection=reject_hook;
    c->SetInnerIter(0);c->SetGlobalParam(MAIN_SOLVER::RANS,RUNTIME_FLOW_SYS);
    // 当前flow密度/分子物性→当前SA涡黏度→flow刷新。无冻结上次mu_t。
    s[FLOW_SOL]->Preprocessing(g,s,c,0,NO_RK_ITER,RUNTIME_FLOW_SYS,false);
    if(nativeRecovery)throw std::runtime_error("NATIVE_RECOVERY_DURING_INPUT_REFRESH");
    s[TURB_SOL]->Postprocessing(g,s,c,0);
    s[FLOW_SOL]->Preprocessing(g,s,c,0,NO_RK_ITER,RUNTIME_FLOW_SYS,false);
    if(nativeRecovery)throw std::runtime_error("NATIVE_RECOVERY_DURING_MUT_REFRESH");
    Space256 spatial;spatial.Space_Integration(g,s,num[FLOW_SOL],c,0,NO_RK_ITER,RUNTIME_FLOW_SYS);
    c->SetGlobalParam(MAIN_SOLVER::RANS,RUNTIME_TURB_SYS);
    s[TURB_SOL]->Preprocessing(g,s,c,0,NO_RK_ITER,RUNTIME_TURB_SYS,false);
    spatial.Space_Integration(g,s,num[TURB_SOL],c,0,NO_RK_ITER,RUNTIME_TURB_SYS);
    if(nativeRecovery)throw std::runtime_error("NATIVE_RECOVERY_DURING_RESIDUAL");assert_unmodified(decoded);
    std::vector<double> f;f.reserve(entries.size());
    for(auto e:entries){double value=e.second<5?s[FLOW_SOL]->LinSysRes(e.first,e.second):s[TURB_SOL]->LinSysRes(e.first,0);
      if(e.second==4)value-=300000*s[FLOW_SOL]->LinSysRes(e.first,0);if(!std::isfinite(value))throw std::runtime_error("NONFINITE_SPACE_RESIDUAL");f.push_back(value);}
    result.swap(f);++physicalCalls;log<<calls<<','<<label<<",0,accepted_operator_evaluation,0\n";log.flush();return 0;
  }
  void initialize(const std::string& folder){
    destination=folder;g=geometry_container[0][0][0];s=solver_container[0][0][0];c=config_container[0];num=numerics_container[0][0][0];
    if(g->GetnPoint()!=335726||g->GetnPointDomain()!=g->GetnPoint()||c->GetnMGLevels()!=0||c->GetNewtonKrylov()||omp_get_max_threads()!=1||
       c->GetKind_FluidModel()!=RAH_PHASE_AIR||c->GetKind_Upwind_Flow()!=UPWIND::ROE||c->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!c->GetSAParsedOptions().ft2)
      throw std::runtime_error("OPERATOR_IDENTITY");
    wall.assign(g->GetnPoint(),0);unsigned long nwall=0;
    for(unsigned m=0;m<c->GetnMarker_All();++m)if(c->GetMarker_All_KindBC(m)==ISOTHERMAL){
      if(c->GetIsothermal_Temperature(c->GetMarker_All_TagBound(m))!=300)throw std::runtime_error("GIVEN_WALL_T");
      for(unsigned v=0;v<g->nVertex[m];++v)wall[g->vertex[m][v]->GetNode()]=1;}
    for(unsigned long i=0;i<g->GetnPoint();++i){nwall+=wall[i];
      for(unsigned j=0;j<6;++j){if(wall[i]&&(j==1||j==2||j==3||j==5)){
          const double value=j<5?s[FLOW_SOL]->GetNodes()->GetSolution(i,j):s[TURB_SOL]->GetNodes()->GetSolution(i,0);if(value!=0)throw std::runtime_error("STRONG_WALL_INPUT_NOT_ZERO");
        }else entries.push_back({i,j});}
    }if(nwall!=9878||entries.size()!=6*g->GetnPoint()-4*nwall)throw std::runtime_error("ORIGINAL_WALL_DOF_COUNT");
    log.open(destination+"/calls.csv");log<<"call,label,status,reason,time_updates\n";
  }
};
inline FullResidualProviderV257* FullResidualProviderV257::active=nullptr;
