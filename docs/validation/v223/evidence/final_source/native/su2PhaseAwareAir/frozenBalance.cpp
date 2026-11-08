// 模块职责：在原完整末场上调用成熟空间算子，记录完整阶段账；绝不推进时间或解线性系统。
#include "drivers/CSinglezoneDriver.hpp"
#include "integration/CIntegration.hpp"
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cstring>

// 仅公开现有受保护入口，不覆盖空间离散实现或更改原生类布局。
class SpatialOnly final : public CIntegration {public: using CIntegration::Space_Integration;};

class FrozenBalanceDriver final : public CSinglezoneDriver {
 public:
  using CSinglezoneDriver::CSinglezoneDriver;
  void inspect(const std::string& folder) {
    auto* cfg=config_container[0];auto* geo=geometry_container[0][0][0];
    auto** sol=solver_container[0][0][0];auto** num=numerics_container[0][0][0][FLOW_SOL];
    auto* flow=sol[FLOW_SOL];auto* sa=sol[TURB_SOL];auto* node=dynamic_cast<CNSVariable*>(flow->GetNodes());
    const auto N=geo->GetnPoint(),D=geo->GetnPointDomain();
    if(!node||!dynamic_cast<CSU2PhaseAwareAir*>(flow->GetFluidModel())||N!=335726||D!=N||sa->GetnVar()!=1)
      throw std::runtime_error("FROZEN_REAL_IDENTITY");
    unsigned wall=999,far=999;
    if(cfg->GetnMarker_All()!=2)throw std::runtime_error("FROZEN_MARKER_COUNT");
    for(unsigned m=0;m<cfg->GetnMarker_All();++m) {
      const auto tag=cfg->GetMarker_All_TagBound(m);
      if(tag=="wall"&&cfg->GetMarker_All_KindBC(m)==ISOTHERMAL)wall=m;
      else if(tag=="farfield"&&cfg->GetMarker_All_KindBC(m)==RIEMANN_BOUNDARY)far=m;
      else throw std::runtime_error("FROZEN_BOUNDARY_IDENTITY");
    }
    if(wall==999||far==999||geo->nVertex[wall]!=9878)throw std::runtime_error("FROZEN_ORIGINAL_WALL");
    cfg->SetInnerIter(199);cfg->SetGlobalParam(MAIN_SOLVER::RANS,RUNTIME_FLOW_SYS);
    std::vector<std::array<double,6>> original(N);
    for(unsigned long i=0;i<N;++i) {
      for(unsigned j=0;j<5;++j)original[i][j]=node->GetSolution(i,j);
      original[i][5]=sa->GetNodes()->GetSolution(i,0);
    }
    auto prepare=[&]() {
      flow->Preprocessing(geo,sol,cfg,MESH_0,NO_RK_ITER,RUNTIME_FLOW_SYS,false);
    };
    std::ofstream stages(folder+"/stages.csv"),fields(folder+"/field.csv"),walls(folder+"/wall.csv");
    stages<<std::setprecision(17)<<"stage,Rmass,Rmx,Rmy,Rmz,REstar,sumAbsMass,sumAbsMx,sumAbsMy,sumAbsMz,sumAbsEstar\n";
    auto observe=[&](const char* name) {
      std::array<long double,5> sum{},absolute{};
      for(unsigned long i=0;i<N;++i)for(unsigned j=0;j<5;++j) {
        const auto r=flow->LinSysRes(i,j);sum[j]+=r;absolute[j]+=std::abs(r);
      }
      stages<<name;for(auto x:sum)stages<<','<<x;for(auto x:absolute)stages<<','<<x;stages<<'\n';stages.flush();
    };
    // 仅展开当前两个物理marker的原生调度；所有通量、梯度、源项和边界公式仍来自原生对象。
    prepare();observe("prepared");
    flow->Upwind_Residual(geo,sol,num,cfg,MESH_0);observe("convective");
    flow->Viscous_Residual(geo,sol,num,cfg,MESH_0,NO_RK_ITER);observe("viscous");
    flow->Source_Residual(geo,sol,num,cfg,MESH_0);observe("source");
    auto* cb=num[CONV_BOUND_TERM];auto* vb=num[VISC_BOUND_TERM];
    flow->BC_Fluid_Interface(geo,sol,cb,vb,cfg);observe("fluid_interface");
    flow->BC_Riemann(geo,sol,cb,vb,cfg,far);observe("farfield");
    flow->Impose_Fixed_Values(geo,cfg);observe("fixed_values");
    std::vector<std::array<double,3>> beforeWall(geo->nVertex[wall]);
    for(unsigned v=0;v<geo->nVertex[wall];++v)for(unsigned j=0;j<3;++j)
      beforeWall[v][j]=flow->LinSysRes(geo->vertex[wall][v]->GetNode(),j+1);
    flow->BC_Isothermal_Wall(geo,sol,cb,vb,cfg,wall);observe("isothermal");
    std::vector<std::array<double,5>> residual(N);
    for(unsigned long i=0;i<N;++i)for(unsigned j=0;j<5;++j)residual[i][j]=flow->LinSysRes(i,j);
    // 使用原生完整调度作同场对照；不求伪时间步，不解Jacobian，不调用Run。
    prepare();SpatialOnly reference;
    reference.Space_Integration(geo,sol,num,cfg,MESH_0,NO_RK_ITER,RUNTIME_FLOW_SYS);
    observe("native_dispatch");
    flow->Pressure_Forces(geo,cfg);flow->Friction_Forces(geo,cfg);
    unsigned long mismatch=0,stateChanged=0;
    for(unsigned long i=0;i<N;++i) {
      for(unsigned j=0;j<5;++j) {
        const double r=flow->LinSysRes(i,j),u=node->GetSolution(i,j);
        mismatch+=std::memcmp(&r,&residual[i][j],sizeof(double))!=0;
        stateChanged+=std::memcmp(&u,&original[i][j],sizeof(double))!=0;
      }
      const double n=sa->GetNodes()->GetSolution(i,0);stateChanged+=std::memcmp(&n,&original[i][5],sizeof(double))!=0;
    }
    fields<<std::setprecision(17)<<"id,volume,rho,rhou,rhov,rhow,rhoEstar,nutilde,Rmass,Rmx,Rmy,Rmz,REstar\n";
    for(unsigned long i=0;i<N;++i) {
      fields<<geo->nodes->GetGlobalIndex(i)<<','<<geo->nodes->GetVolume(i);
      for(auto v:original[i])fields<<','<<v;for(auto v:residual[i])fields<<','<<v;fields<<'\n';
    }
    walls<<std::setprecision(17)<<"id,area,nx,ny,nz,normalID,distance,T,p,rho,k,Tnear,preWallRx,preWallRy,preWallRz,REstar,q,cfX,cfY,cfZ\n";
    for(unsigned v=0;v<geo->nVertex[wall];++v) {
      const auto i=geo->vertex[wall][v]->GetNode(),j=geo->vertex[wall][v]->GetNormal_Neighbor();
      const auto* normal=geo->vertex[wall][v]->GetNormal();double area=0,projection=0;
      for(unsigned d=0;d<3;++d){area+=normal[d]*normal[d];projection+=normal[d]*(geo->nodes->GetCoord(j,d)-geo->nodes->GetCoord(i,d));}
      area=std::sqrt(area);walls<<geo->nodes->GetGlobalIndex(i)<<','<<area;
      for(unsigned d=0;d<3;++d)walls<<','<<normal[d];
      walls<<','<<geo->nodes->GetGlobalIndex(j)<<','<<std::abs(projection)/area<<','<<node->GetTemperature(i)<<','<<node->GetPressure(i)
           <<','<<node->GetDensity(i)<<','<<node->GetThermalConductivity(i)<<','<<node->GetTemperature(j);
      for(auto r:beforeWall[v])walls<<','<<r;walls<<','<<residual[i][4]<<','<<flow->GetHeatFlux(wall,v);
      for(unsigned d=0;d<3;++d)walls<<','<<flow->GetCSkinFriction(wall,v,d);walls<<'\n';
    }
    std::ofstream result(folder+"/native.json");result<<"{\"points\":"<<N<<",\"stagedNativeRHSMismatches\":"<<mismatch
      <<",\"currentStateBitChanges\":"<<stateChanged<<",\"spaceEvaluations\":2,\"flowUpdates\":0,\"forceRefN\":"
      <<std::setprecision(17)<<flow->GetReferenceDynamicPressure()*cfg->GetRefArea()<<",\"bodyForceCoefficients\":["
      <<flow->GetTotal_CFx()<<','<<flow->GetTotal_CFy()<<','<<flow->GetTotal_CFz()<<"]}\n";
    std::cout<<"Frozen original rocket: mismatch="<<mismatch<<" stateChanged="<<stateChanged<<std::endl;
    if(mismatch||stateChanged)throw std::runtime_error("FROZEN_COMPARISON_FAILED");
  }
};

int main(int argc,char** argv) {
  if(argc!=3)return 2;
  try {SU2_MPI::Init(&argc,&argv);FrozenBalanceDriver driver(argv[1],1,SU2_MPI::GetComm());driver.inspect(argv[2]);driver.Finalize();SU2_MPI::Finalize();}
  catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
  return 0;
}
