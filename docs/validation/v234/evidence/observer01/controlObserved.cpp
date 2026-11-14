// 模块职责：同一制造场的真实单区/双区接口、梯度、热壁与接受控制；不推进流动。
#include "drivers/CSinglezoneDriver.hpp"
#include "drivers/CMultizoneDriver.hpp"
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cstring>


// 仅记录未捕获异常和确切调用阶段，禁止将异常退出当作成功。
namespace failure_observer {
const char* phase="constructor";int scenario=-1,zone=-1;
void failed(){std::string message="unknown";try{auto p=std::current_exception();if(p)std::rethrow_exception(p);}
catch(const std::exception& e){message=e.what();}catch(...){message="non-standard exception";}
std::ofstream f("termination.txt");f<<"phase="<<phase<<"\nscenario="<<scenario<<"\nzone="<<zone<<"\nwhat="<<message<<"\n";f.flush();std::_Exit(86);}
}

struct Zone {
  CConfig* cfg; CGeometry* geo; CSolver** sol; CNumerics*** num;
  CSolver *flow,*sa; CNSVariable* nodes;
};

// 模块一：仅暴露成熟容器与传递入口，类布局、边界和推进实现不变。
template<class Base> class ControlAccess : public Base {
 public:
  using Base::Base;
  Zone zone(unsigned z){auto** s=this->solver_container[z][0][0];return {
    this->config_container[z],this->geometry_container[z][0][0],s,this->numerics_container[z][0][0],
    s[FLOW_SOL],s[TURB_SOL],dynamic_cast<CNSVariable*>(s[FLOW_SOL]->GetNodes())};}
};
class Multi final:public ControlAccess<CMultizoneDriver>{
 public:using ControlAccess::ControlAccess;
 void exchange(){TransferData(0,1);TransferData(1,0);}
};
class Single final:public ControlAccess<CSinglezoneDriver>{
 public:using ControlAccess::ControlAccess;void exchange(){}
};

void install(Zone z,int scenario){
  if(!z.nodes||!z.sa||z.sa->GetnVar()!=1||z.flow->GetnVar()!=5||
     !dynamic_cast<CSU2PhaseAwareAir*>(z.flow->GetFluidModel())||
     z.cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!z.cfg->GetSAParsedOptions().ft2)
    throw std::runtime_error("MULTIZONE_PHASE_SA_IDENTITY");
  const auto ref=rah_phase_air::PT(1200.,350.123);double nu=3*ref.mu/ref.rho,sign=scenario==2?-1.:1.;
  for(unsigned long i=0;i<z.geo->GetnPoint();++i){
    const auto* x=z.geo->nodes->GetCoord(i);bool affine=scenario!=0;
    const auto q=rah_phase_air::PT(1200.,350.123+(affine? .25*x[0]+1.5*x[1]-.75*x[2]:0.));
    double v[3]={40.+(affine?4*x[1]:0),30.+(affine?3*x[2]:0),-20.+(affine?2*x[0]:0)};
    double U[5]={q.rho,q.rho*v[0],q.rho*v[1],q.rho*v[2],q.rho*(q.e+.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))};
    z.nodes->SetSolution(i,U);z.nodes->SetSolution_Old(i,U);
    if(!z.nodes->SetPrimVar(i,0,0,z.flow->GetFluidModel()))throw std::runtime_error("MANUFACTURED_STATE_REJECTED");
    z.sa->GetNodes()->SetSolution(i,0,sign*nu*(1+(affine?.05*x[0]:0)));
    z.sa->GetNodes()->SetSolution_Old(i,0,z.sa->GetNodes()->GetSolution(i,0));
  }
  z.sa->Postprocessing(z.geo,z.sol,z.cfg,0);
  for(unsigned long i=0;i<z.geo->GetnPoint();++i){
    z.nodes->SetPrimVar(i,z.sa->GetNodes()->GetmuT(i),0,z.flow->GetFluidModel());
    z.nodes->SetSecondaryVar(i,z.flow->GetFluidModel());
  }
}

void gradients(Zone z,int scenario){
  z.cfg->SetGlobalParam(MAIN_SOLVER::RANS,RUNTIME_FLOW_SYS);
  z.flow->Preprocessing(z.geo,z.sol,z.cfg,0,NO_RK_ITER,RUNTIME_FLOW_SYS,false);
  z.sa->Preprocessing(z.geo,z.sol,z.cfg,0,NO_RK_ITER,RUNTIME_TURB_SYS,false);
  if(scenario==3)return; // 第四场保留真实原生梯度，供接口离散敏感性对照。
  const auto ref=rah_phase_air::PT(1200.,350.123);double nu=3*ref.mu/ref.rho,sign=scenario==2?-1.:1.;
  const double gT[3]={.25,1.5,-.75},A[3][3]={{0,4,0},{0,0,3},{2,0,0}};
  for(unsigned long i=0;i<z.geo->GetnPoint();++i){
    auto g=z.nodes->GetGradient_Primitive(i);auto sg=z.sa->GetNodes()->GetGradient(i);
    for(unsigned j=0;j<z.nodes->GetnPrimVarGrad();++j)for(int d=0;d<3;++d)g[j][d]=0;
    for(int d=0;d<3;++d)sg[0][d]=0;
    if(scenario==0)continue;
    for(int d=0;d<3;++d){
      g[0][d]=gT[d];for(int j=0;j<3;++j)g[j+1][d]=A[j][d];
      g[5][d]=-z.nodes->GetDensity(i)/z.nodes->GetTemperature(i)*gT[d];
      if(z.nodes->GetnPrimVarGrad()>6){g[6][d]=z.nodes->GetSpecificHeatCp(i)*gT[d];
        for(int j=0;j<3;++j)g[6][d]+=z.nodes->GetVelocity(i,j)*A[j][d];}
      sg[0][d]=d==0?sign*nu*.05:0;
    }
  }
}

void clear(Zone z){z.flow->LinSysRes.SetValZero();z.flow->Jacobian.SetValZero();z.sa->LinSysRes.SetValZero();z.sa->Jacobian.SetValZero();}

template<class Driver> void control(Driver& driver,unsigned nz,const std::string& folder){
  std::ofstream nodes(folder+"/nodes.csv"),interfaces(folder+"/interfaces.csv"),stages(folder+"/stages.csv"),walls(folder+"/walls.csv"),recovery(folder+"/recovery.csv");
  for(auto* out:{&nodes,&interfaces,&stages,&walls,&recovery})*out<<std::setprecision(17);
  nodes<<"scenario,zone,localID,globalID,x,y,z,volume,wallDistance";
  for(int j=0;j<5;++j)nodes<<",U"<<j;nodes<<",SA";for(int j=0;j<12;++j)nodes<<",P"<<j;
  for(int j=0;j<7;++j)for(int d=0;d<3;++d)nodes<<",G"<<j<<d;nodes<<'\n';
  interfaces<<"scenario,zone,node,x,y,z,nx,ny,nz,donors,flowWeight,saWeight,flowDonorMaxDiff,saDonorMaxDiff,R0,R1,R2,R3,R4,Rsa\n";
  stages<<"scenario,zone,stage,node,x,y,z,R0,R1,R2,R3,R4,Rsa\n";
  walls<<"scenario,zone,node,x,y,z,nx,ny,nz,area,distance,Tnear,k,TwallFluid,REwall,q,preRx,preRy,preRz,postRx,postRy,postRz,oldMx,oldMy,oldMz,oldSA\n";
  recovery<<"zone,trial,returned,threw,flowRestoredBits,SAUnchangedBits\n";
  auto record=[&](Zone z,unsigned iz,int scenario,const char* stage){
    for(unsigned long i=0;i<z.geo->GetnPointDomain();++i){stages<<scenario<<','<<iz<<','<<stage<<','<<i;
      for(int d=0;d<3;++d)stages<<','<<z.geo->nodes->GetCoord(i,d);
      for(int j=0;j<5;++j)stages<<','<<z.flow->LinSysRes(i,j);stages<<','<<z.sa->LinSysRes(i,0)<<'\n';}
  };
  unsigned long stateChanges=0;
  for(int scenario=0;scenario<4;++scenario){failure_observer::scenario=scenario;failure_observer::phase="install";
    for(unsigned z=0;z<nz;++z)install(driver.zone(z),scenario);
    failure_observer::phase="exchange-and-gradient";driver.exchange();for(unsigned z=0;z<nz;++z)gradients(driver.zone(z),scenario);driver.exchange();
    for(unsigned iz=0;iz<nz;++iz){failure_observer::zone=iz;failure_observer::phase="dump-state";auto z=driver.zone(iz);auto* fn=z.flow->GetNodes();auto* sn=z.sa->GetNodes();
      std::vector<std::array<double,6>> original(z.geo->GetnPoint());
      for(unsigned long i=0;i<z.geo->GetnPoint();++i){nodes<<scenario<<','<<iz<<','<<i<<','<<z.geo->nodes->GetGlobalIndex(i);
        for(int d=0;d<3;++d)nodes<<','<<z.geo->nodes->GetCoord(i,d);nodes<<','<<z.geo->nodes->GetVolume(i)<<','<<z.geo->nodes->GetWall_Distance(i);
        for(int j=0;j<5;++j){original[i][j]=fn->GetSolution(i,j);nodes<<','<<original[i][j];}
        original[i][5]=sn->GetSolution(i,0);nodes<<','<<original[i][5];
        for(int j=0;j<12;++j)nodes<<','<<fn->GetPrimitive(i,j);
        auto g=fn->GetGradient_Primitive(i);for(int j=0;j<7;++j)for(int d=0;d<3;++d)nodes<<','<<g[j][d];nodes<<'\n';}
      failure_observer::phase="interface";clear(z);auto** nf=z.num[FLOW_SOL];auto** ns=z.num[TURB_SOL];
      z.flow->BC_Fluid_Interface(z.geo,z.sol,nf[CONV_BOUND_TERM],nf[VISC_BOUND_TERM],z.cfg);
      z.sa->BC_Fluid_Interface(z.geo,z.sol,ns[CONV_BOUND_TERM],ns[VISC_BOUND_TERM],z.cfg);
      for(unsigned m=0;m<z.cfg->GetnMarker_All();++m)if(z.cfg->GetMarker_All_KindBC(m)==FLUID_INTERFACE){
        for(unsigned v=0;v<z.geo->nVertex[m];++v){auto i=z.geo->vertex[m][v]->GetNode();interfaces<<scenario<<','<<iz<<','<<i;
          for(int d=0;d<3;++d)interfaces<<','<<z.geo->nodes->GetCoord(i,d);
          for(int d=0;d<3;++d)interfaces<<','<<-z.geo->vertex[m][v]->GetNormal(d);
          unsigned donors=z.flow->GetnSlidingStates(m,v);double wf=0,ws=0,df=0,ds=0;
          for(unsigned k=0;k<donors;++k){wf+=z.flow->GetSlidingState(m,v,12,k);
            for(int j=0;j<12;++j)df=std::max(df,std::abs(z.flow->GetSlidingState(m,v,j,k)-fn->GetPrimitive(i,j)));}
          for(unsigned k=0;k<z.sa->GetnSlidingStates(m,v);++k){ws+=z.sa->GetSlidingState(m,v,1,k);ds=std::max(ds,std::abs(z.sa->GetSlidingState(m,v,0,k)-sn->GetSolution(i,0)));}
          interfaces<<','<<donors<<','<<wf<<','<<ws<<','<<df<<','<<ds;
          for(int j=0;j<5;++j)interfaces<<','<<z.flow->LinSysRes(i,j);interfaces<<','<<z.sa->LinSysRes(i,0)<<'\n';
        }
      }
      // 模块三：原生空间/特征/热壁分阶段；不调用时间积分或线性求解。
      failure_observer::phase="interior-flow";clear(z);z.flow->Upwind_Residual(z.geo,z.sol,nf,z.cfg,0);z.flow->Viscous_Residual(z.geo,z.sol,nf,z.cfg,0,NO_RK_ITER);
      failure_observer::phase="interior-SA";z.sa->Upwind_Residual(z.geo,z.sol,ns,z.cfg,0);z.sa->Viscous_Residual(z.geo,z.sol,ns,z.cfg,0,NO_RK_ITER);
      z.flow->BC_Fluid_Interface(z.geo,z.sol,nf[CONV_BOUND_TERM],nf[VISC_BOUND_TERM],z.cfg);
      z.sa->BC_Fluid_Interface(z.geo,z.sol,ns[CONV_BOUND_TERM],ns[VISC_BOUND_TERM],z.cfg);
      record(z,iz,scenario,"interior_plus_interface");
      for(unsigned m=0;m<z.cfg->GetnMarker_All();++m)if(z.cfg->GetMarker_All_KindBC(m)==RIEMANN_BOUNDARY){
        failure_observer::phase="farfield-flow";z.flow->BC_Riemann(z.geo,z.sol,nf[CONV_BOUND_TERM],nf[VISC_BOUND_TERM],z.cfg,m);
        failure_observer::phase="farfield-SA";z.sa->BC_Riemann(z.geo,z.sol,ns[CONV_BOUND_TERM],ns[VISC_BOUND_TERM],z.cfg,m);}
      record(z,iz,scenario,"before_wall");
      for(unsigned m=0;m<z.cfg->GetnMarker_All();++m)if(z.cfg->GetMarker_All_KindBC(m)==ISOTHERMAL){
        std::vector<std::array<double,5>> before(z.geo->nVertex[m]);
        for(unsigned v=0;v<z.geo->nVertex[m];++v){auto i=z.geo->vertex[m][v]->GetNode();for(int j=0;j<5;++j)before[v][j]=z.flow->LinSysRes(i,j);}
        failure_observer::phase="thermal-wall-flow";z.flow->BC_Isothermal_Wall(z.geo,z.sol,nf[CONV_BOUND_TERM],nf[VISC_BOUND_TERM],z.cfg,m);
        failure_observer::phase="thermal-wall-SA";z.sa->BC_Isothermal_Wall(z.geo,z.sol,ns[CONV_BOUND_TERM],ns[VISC_BOUND_TERM],z.cfg,m);
        failure_observer::phase="surface-loads";z.flow->Pressure_Forces(z.geo,z.cfg);z.flow->Friction_Forces(z.geo,z.cfg);
        for(unsigned v=0;v<z.geo->nVertex[m];++v){auto* vertex=z.geo->vertex[m][v];auto i=vertex->GetNode(),j=vertex->GetNormal_Neighbor();
          double area=0,distance=0;for(int d=0;d<3;++d){double n=vertex->GetNormal(d);area+=n*n;distance+=n*(z.geo->nodes->GetCoord(j,d)-z.geo->nodes->GetCoord(i,d));}
          area=std::sqrt(area);walls<<scenario<<','<<iz<<','<<i;for(int d=0;d<3;++d)walls<<','<<z.geo->nodes->GetCoord(i,d);
          for(int d=0;d<3;++d)walls<<','<<-vertex->GetNormal(d);
          walls<<','<<area<<','<<std::abs(distance)/area<<','<<z.nodes->GetTemperature(j)<<','<<z.nodes->GetThermalConductivity(i)<<','<<z.nodes->GetTemperature(i)
            <<','<<z.flow->LinSysRes(i,4)-before[v][4]<<','<<z.flow->GetHeatFlux(m,v);
          for(int d=1;d<4;++d)walls<<','<<before[v][d];for(int d=1;d<4;++d)walls<<','<<z.flow->LinSysRes(i,d);
          for(int d=1;d<4;++d)walls<<','<<fn->GetSolution_Old(i,d);walls<<','<<sn->GetSolution_Old(i,0)<<'\n';}
      }
      record(z,iz,scenario,"after_wall");
      for(unsigned long i=0;i<z.geo->GetnPoint();++i){for(int j=0;j<5;++j){double u=fn->GetSolution(i,j);stateChanges+=std::memcmp(&u,&original[i][j],8)!=0;}
        double u=sn->GetSolution(i,0);stateChanges+=std::memcmp(&u,&original[i][5],8)!=0;}
    }
  }
  // 模块四：已存在CNS拒绝路径及SA独立身份，只检查具名试探，不重写恢复机制。
  for(unsigned iz=0;iz<nz;++iz){auto z=driver.zone(iz);install(z,2);const unsigned i=0;std::array<double,5> good{},bad{};
    for(int j=0;j<5;++j)good[j]=bad[j]=z.nodes->GetSolution(i,j);double sa=z.sa->GetNodes()->GetSolution(i,0);bad[4]=-1e8;
    for(int invalidOld=0;invalidOld<2;++invalidOld){z.nodes->SetSolution(i,bad.data());z.nodes->SetSolution_Old(i,invalidOld?bad.data():good.data());
      int returned=-1,threw=0;try{returned=z.nodes->SetPrimVar(i,0,0,z.flow->GetFluidModel());}catch(const std::domain_error&){threw=1;}
      int restored=0;for(int j=0;j<5;++j){double value=z.nodes->GetSolution(i,j);restored+=std::memcmp(&value,&good[j],8)==0;}
      double currentSA=z.sa->GetNodes()->GetSolution(i,0);recovery<<iz<<','<<invalidOld<<','<<returned<<','<<threw<<','<<restored<<','<<(std::memcmp(&sa,&currentSA,8)==0)<<'\n';}
    install(z,0);
  }
  if(!nodes||!interfaces||!stages||!walls||!recovery)throw std::runtime_error("MULTIZONE_CONTROL_WRITE_FAILED");
  std::ofstream result(folder+"/native.json");result<<"{\"zones\":"<<nz<<",\"scenarios\":4,\"flowUpdates\":0,\"currentStateBitChanges\":"<<stateChanges<<",\"actualTransferDataCalled\":"<<(nz==2?"true":"false")<<"}\n";
  if(stateChanges)throw std::runtime_error("MANUFACTURED_CURRENT_STATE_CHANGED");
}

int main(int argc,char** argv){try{
  std::set_terminate(failure_observer::failed);SU2_MPI::Init(&argc,&argv);if(argc!=4)throw std::runtime_error("CONFIG_OUTPUT_MODE_REQUIRED");
  if(std::string(argv[3])=="multi"){Multi d(argv[1],2,SU2_MPI::GetComm());control(d,2,argv[2]);d.Finalize();}
  else if(std::string(argv[3])=="single"){Single d(argv[1],1,SU2_MPI::GetComm());control(d,1,argv[2]);d.Finalize();}
  else throw std::runtime_error("MODE_NOT_ADMITTED");
  SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
