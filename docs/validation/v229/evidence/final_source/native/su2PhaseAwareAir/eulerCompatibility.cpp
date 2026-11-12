// 模块一：真实Euler分派与固定三维控制；不调用时间推进或稳态迭代。
#include "drivers/CSinglezoneDriver.hpp"
#include "variables/CEulerVariable.hpp"
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

class EulerCompatibility final : public CSinglezoneDriver {
 public:
  using CSinglezoneDriver::CSinglezoneDriver;
  void check(const std::string& folder, bool restart) {
    auto* cfg=config_container[0];auto* geo=geometry_container[0][0][0];
    auto** sol=solver_container[0][0][0];auto*** num=numerics_container[0][0][0];
    auto* flow=sol[FLOW_SOL];auto* node=dynamic_cast<CEulerVariable*>(flow->GetNodes());
    if(!node || dynamic_cast<CNSVariable*>(node) || flow->GetnVar()!=5 ||
       cfg->GetKind_FluidModel()!=RAH_PHASE_AIR_EULER || cfg->GetKind_Turb_Model()!=TURB_MODEL::NONE ||
       !dynamic_cast<CSU2PhaseAwareAir*>(flow->GetFluidModel()))throw std::runtime_error("EULER_DISPATCH_REJECTED");
    const auto N=geo->GetnPoint();if(N!=24)throw std::runtime_error("CONTROL_TOPOLOGY");
    auto dump=[&](const std::string& name) {
      std::ofstream f(folder+"/"+name);f<<std::setprecision(17)<<"point,x,y,z,rho,rhou,rhov,rhow,rhoE,p,T\n";
      for(unsigned long i=0;i<N;++i){f<<i;for(int d=0;d<3;++d)f<<','<<geo->nodes->GetCoord(i,d);
        for(int j=0;j<5;++j)f<<','<<node->GetSolution(i,j);
        f<<','<<node->GetPressure(i)<<','<<node->GetTemperature(i)<<'\n';}
      if(!f)throw std::runtime_error("DUMP_WRITE");
    };
    dump(restart?"restarted.csv":"initial.csv");if(restart)return;
    unsigned wall=999,far=999;
    for(unsigned m=0;m<cfg->GetnMarker_All();++m){
      if(cfg->GetMarker_All_TagBound(m)=="wall")wall=m;
      if(cfg->GetMarker_All_TagBound(m)=="farfield")far=m;}
    if(wall==999||far==999)throw std::runtime_error("CONTROL_MARKER");
    const auto q=rah_phase_air::PT(1200.,350.123);const double a=std::sqrt(q.a2);
    double v[3]={2*a,.4*a,.3*a};double U[5]={q.rho,q.rho*v[0],q.rho*v[1],q.rho*v[2],
      q.rho*(q.e+.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))};
    auto install=[&](){
      for(unsigned long i=0;i<N;++i){node->SetSolution(i,U);node->SetSolution_Old(i,U);
        if(!node->SetPrimVar(i,flow->GetFluidModel()))throw std::runtime_error("VALID_PRIMITIVE_REJECTED");
        node->SetSecondaryVar(i,flow->GetFluidModel());}
      flow->LinSysRes.SetValZero();flow->Jacobian.SetValZero();};
    install();
    // 模块二：同一制造场的四类一般EOS特征边界，不复写边界计算。
    flow->BC_Riemann(geo,sol,num[FLOW_SOL][CONV_BOUND_TERM],num[FLOW_SOL][VISC_BOUND_TERM],cfg,far);
    std::ofstream b(folder+"/boundary.csv");b<<std::setprecision(17)<<"node,nx,ny,nz,T,u,v,w,p,rho,H,a,R0,R1,R2,R3,R4\n";
    for(unsigned j=0;j<geo->nVertex[far];++j){auto i=geo->vertex[far][j]->GetNode();
      auto* normal=geo->vertex[far][j]->GetNormal();auto* p=flow->GetCharacPrimVar(far,j);b<<i;
      for(int d=0;d<3;++d)b<<','<<-normal[d];for(int k=0;k<8;++k)b<<','<<p[k];
      for(int k=0;k<5;++k)b<<','<<flow->LinSysRes(i,k);b<<'\n';}
    // 模块三：实际静止滑移投影。固定非零残差/单位对角，比较投影而非伪造流动解。
    install();const double residual[5]={.17,.31,-.27,.19,17.};
    for(unsigned long i=0;i<N;++i){for(int k=0;k<5;++k)flow->LinSysRes(i,k)=residual[k];
      auto* block=flow->Jacobian.GetBlock(i,i);for(int k=0;k<5;++k)block[k*5+k]=1.;}
    flow->BC_Euler_Wall(geo,sol,num[FLOW_SOL][CONV_BOUND_TERM],num[FLOW_SOL][VISC_BOUND_TERM],cfg,wall);
    std::ofstream w(folder+"/wall.csv");w<<std::setprecision(17)<<"node,nx,ny,nz";
    for(int k=0;k<5;++k)w<<",U"<<k;for(int k=0;k<5;++k)w<<",R"<<k;
    for(int k=0;k<25;++k)w<<",J"<<k;w<<",p,T\n";
    for(unsigned j=0;j<geo->nVertex[wall];++j){auto i=geo->vertex[wall][j]->GetNode();
      const auto* normal=geo->vertex[wall][j]->GetNormal();double n[3];double area=0.;
      for(int d=0;d<3;++d)area+=normal[d]*normal[d];area=std::sqrt(area);
      auto it=geo->symmetryNormals[wall].find(j);
      for(int d=0;d<3;++d)n[d]=(it!=geo->symmetryNormals[wall].end())?it->second[d]:normal[d]/area;
      w<<i;for(int d=0;d<3;++d)w<<','<<n[d];
      node->SetSolution(i,node->GetSolution_Old(i));
      if(!node->SetPrimVar(i,flow->GetFluidModel()))throw std::runtime_error("PROJECTED_STATE_REJECTED");
      for(int k=0;k<5;++k)w<<','<<node->GetSolution(i,k);
      for(int k=0;k<5;++k)w<<','<<flow->LinSysRes(i,k);
      const auto* block=flow->Jacobian.GetBlock(i,i);for(int k=0;k<25;++k)w<<','<<block[k];
      w<<','<<node->GetPressure(i)<<','<<node->GetTemperature(i)<<'\n';}
    // 模块四：非法守恒态直接经过生产SetPrimVar；记录拒绝，不剪裁后继续。
    std::ofstream r(folder+"/refusals.csv");r<<"case,rejected\n";
    double valid[5];for(int k=0;k<5;++k)valid[k]=node->GetSolution(0,k);
    for(int kind=0;kind<2;++kind){double bad[5];for(int k=0;k<5;++k)bad[k]=valid[k];
      if(kind==0)bad[0]=-1.;else bad[4]=valid[0]*1e8;
      node->SetSolution(0,bad);bool rejected=false;
      try{rejected=!node->SetPrimVar(0,flow->GetFluidModel());}catch(const std::domain_error&){rejected=true;}
      r<<(kind==0?"negative_density":"out_of_domain_energy")<<','<<rejected<<'\n';
      if(!rejected)throw std::runtime_error("INVALID_STATE_ACCEPTED");
      node->SetSolution(0,valid);node->SetSolution_Old(0,valid);
      if(!node->SetPrimVar(0,flow->GetFluidModel()))throw std::runtime_error("VALID_STATE_RESTORE_FAILED");}
    dump("saved.csv");
    if(!output_container[0]->SetResultFiles(geo,cfg,sol,0,true))throw std::runtime_error("NATIVE_OUTPUT_MISSING");
    if(!b||!w||!r)throw std::runtime_error("CONTROL_WRITE");
    std::ofstream status(folder+"/dispatch.json");
    status<<"{\"passed\":true,\"fluidEnum\":16,\"conservativeVariables\":5,\"nodes\":24,\"CFDIterations\":0,\"heatFluxSolved\":false}\n";
  }
};
int main(int argc,char** argv){try{SU2_MPI::Init(&argc,&argv);
  if(argc!=4)throw std::runtime_error("CONFIG_OUTPUT_MODE_REQUIRED");
  EulerCompatibility driver(argv[1],1,SU2_MPI::GetComm());driver.check(argv[2],std::string(argv[3])=="restart");
  driver.Finalize();SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
