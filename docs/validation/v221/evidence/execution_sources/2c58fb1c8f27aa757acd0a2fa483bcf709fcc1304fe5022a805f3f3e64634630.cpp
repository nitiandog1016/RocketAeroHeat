// 模块职责：经真实driver构造一次完整小域，调用生产边界/SA/输出接口。
// 只检查固定制造状态，不调用Run/StartSolver，不推进任何物理流场。
#include "drivers/CSinglezoneDriver.hpp"
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

class CompatibilityDriver final : public CSinglezoneDriver {
 public:
  using CSinglezoneDriver::CSinglezoneDriver;
  void check(const std::string& folder,bool restart) {
    auto* cfg=config_container[0];auto* geo=geometry_container[0][0][0];
    auto** sol=solver_container[0][0][0];auto*** num=numerics_container[0][0][0];
    auto* flow=sol[FLOW_SOL];auto* sa=sol[TURB_SOL];
    auto* node=dynamic_cast<CNSVariable*>(flow->GetNodes());auto* tn=sa->GetNodes();
    if(!node||!dynamic_cast<CSU2PhaseAwareAir*>(flow->GetFluidModel())||tn->GetnVar()!=1)
      throw std::runtime_error("wrong actual driver dispatch");
    if(cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!cfg->GetSAParsedOptions().ft2)
      throw std::runtime_error("wrong actual SA dispatch");
    const auto N=geo->GetnPoint();
    auto dump=[&](const std::string& name) {
      std::ofstream f(folder+"/"+name);f<<std::setprecision(17)<<"point,x,y,z,rho,rhou,rhov,rhow,rhoE,nutilde,p,T\n";
      for(unsigned long i=0;i<N;++i) {
        f<<i;for(int j=0;j<3;++j)f<<','<<geo->nodes->GetCoord(i,j);
        for(int j=0;j<5;++j)f<<','<<node->GetSolution(i,j);
        f<<','<<tn->GetSolution(i,0)<<','<<node->GetPressure(i)<<','<<node->GetTemperature(i)<<'\n';
      }
    };
    dump(restart?"restarted.csv":"initial.csv");
    if(restart)return;
    unsigned wall=999,far=999;
    for(unsigned m=0;m<cfg->GetnMarker_All();++m) {
      if(cfg->GetMarker_All_TagBound(m)=="wall")wall=m;
      if(cfg->GetMarker_All_TagBound(m)=="farfield")far=m;
    }
    if(wall==999||far==999||N!=24)throw std::runtime_error("control topology differs");
    auto install=[&](double sign,bool final) {
      for(unsigned long i=0;i<N;++i) {
        const auto q=rah_phase_air::PT(1200.,350.123);double a=std::sqrt(q.a2);
        double v[3]={2*a,.4*a,.3*a};if(final)for(auto& x:v)x=0.;
        double U[5]={q.rho,q.rho*v[0],q.rho*v[1],q.rho*v[2],q.rho*(q.e+.5*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))};
        node->SetSolution(i,U);node->SetSolution_Old(i,U);
        if(!node->SetPrimVar(i,0,0,flow->GetFluidModel()))throw std::runtime_error("manufactured state rejected");
        tn->SetSolution(i,0,sign*3*q.mu/q.rho);tn->SetSolution_Old(i,0,tn->GetSolution(i,0));
      }
      sa->Postprocessing(geo,sol,cfg,0);
      for(unsigned long i=0;i<N;++i) {
        node->SetPrimVar(i,tn->GetmuT(i),0,flow->GetFluidModel());node->SetSecondaryVar(i,flow->GetFluidModel());
        auto grad=node->GetGradient_Primitive(i);
        for(unsigned j=0;j<node->GetnPrimVarGrad();++j)for(int d=0;d<3;++d)grad[j][d]=0.;
      }
      flow->LinSysRes.SetValZero();flow->Jacobian.SetValZero();sa->LinSysRes.SetValZero();sa->Jacobian.SetValZero();
    };
    std::ofstream boundary(folder+"/boundary.csv"),walls(folder+"/wall.csv"),source(folder+"/sa_source.csv");
    boundary<<std::setprecision(17)<<"sign,node,nx,ny,nz,T,u,v,w,p,rho,H,mu,mut,k,cp,saResidual,flowR0,flowR1,flowR2,flowR3,flowR4\n";
    walls<<std::setprecision(17)<<"sign,node,area,distance,Tnear,rho,E,u2,cp,cv,k,mut,RE,J40,J44,oldMx,oldMy,oldMz,oldSA,saResidual,saDiag\n";
    source<<std::setprecision(17)<<"sign,rho,mu,nutilde,mut,distance,volume,omega,residual,jacobian\n";
    for(double sign:{1.,-1.}) {
      install(sign,false);
      // 实际源项对象来自driver原生SA工厂；只设置公开输入，不复写源项公式。
      auto* sn=num[TURB_SOL][SOURCE_FIRST_TERM];double vort[3]={0,0,120};
      sn->SetPrimitive(node->GetPrimitive(0),node->GetPrimitive(0));
      sn->SetScalarVar(tn->GetSolution(0),tn->GetSolution(0));
      auto sg=tn->GetGradient(0);for(int d=0;d<3;++d)sg[0][d]=0;
      sn->SetScalarVarGradient(sg,sg);sn->SetPrimVarGradient(node->GetGradient_Primitive(0),node->GetGradient_Primitive(0));
      sn->SetVorticity(vort,vort);sn->SetStrainMag(120,120);sn->SetDistance(.01,.01);sn->SetRoughness(0,0);sn->SetVolume(.002);
      const auto sr=sn->ComputeResidual(cfg);
      source<<sign<<','<<node->GetDensity(0)<<','<<node->GetLaminarViscosity(0)<<','<<tn->GetSolution(0,0)<<','<<tn->GetmuT(0)
        <<",0.01,0.002,120,"<<sr[0]<<','<<sr.jacobian_i[0][0]<<'\n';
      flow->BC_Riemann(geo,sol,num[FLOW_SOL][CONV_BOUND_TERM],num[FLOW_SOL][VISC_BOUND_TERM],cfg,far);
      sa->BC_Riemann(geo,sol,num[TURB_SOL][CONV_BOUND_TERM],num[TURB_SOL][VISC_BOUND_TERM],cfg,far);
      for(unsigned v=0;v<geo->nVertex[far];++v) {
        auto i=geo->vertex[far][v]->GetNode();auto* normal=geo->vertex[far][v]->GetNormal();
        auto* b=flow->GetCharacPrimVar(far,v);boundary<<sign<<','<<i;
        for(int d=0;d<3;++d)boundary<<','<<-normal[d];
        for(int j=0;j<7;++j)boundary<<','<<b[j];
        for(int j=8;j<12;++j)boundary<<','<<b[j];
        boundary<<','<<sa->LinSysRes(i,0);for(int j=0;j<5;++j)boundary<<','<<flow->LinSysRes(i,j);boundary<<'\n';
      }
      flow->LinSysRes.SetValZero();flow->Jacobian.SetValZero();sa->LinSysRes.SetValZero();sa->Jacobian.SetValZero();
      flow->BC_Isothermal_Wall(geo,sol,num[FLOW_SOL][CONV_BOUND_TERM],num[FLOW_SOL][VISC_BOUND_TERM],cfg,wall);
      sa->BC_Isothermal_Wall(geo,sol,num[TURB_SOL][CONV_BOUND_TERM],num[TURB_SOL][VISC_BOUND_TERM],cfg,wall);
      for(unsigned v=0;v<geo->nVertex[wall];++v) {
        auto* vertex=geo->vertex[wall][v];auto i=vertex->GetNode(),j=vertex->GetNormal_Neighbor();
        auto* normal=vertex->GetNormal();double area=std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
        double distance=0,u2=0;for(int d=0;d<3;++d){distance+=normal[d]*(geo->nodes->GetCoord(j,d)-geo->nodes->GetCoord(i,d))/area;u2+=node->GetVelocity(i,d)*node->GetVelocity(i,d);}
        const auto* jac=flow->Jacobian.GetBlock(i,i);const auto q=rah_phase_air::PT(node->GetPressure(i),node->GetTemperature(i));
        walls<<sign<<','<<i<<','<<area<<','<<std::abs(distance)<<','<<node->GetTemperature(j)<<','<<node->GetDensity(i)<<','<<node->GetEnergy(i)
          <<','<<u2<<','<<node->GetSpecificHeatCp(i)<<','<<q.cv<<','<<node->GetThermalConductivity(i)<<','<<node->GetEddyViscosity(i)
          <<','<<flow->LinSysRes(i,4)<<','<<jac[20]<<','<<jac[24];
        for(int d=1;d<4;++d)walls<<','<<node->GetSolution_Old(i,d);
        walls<<','<<tn->GetSolution_Old(i,0)<<','<<sa->LinSysRes(i,0)<<','<<sa->Jacobian.GetBlock(i,i)[0]<<'\n';
      }
    }
    // 模块三：由实际输出器保存同一制造状态，供第二进程真实重启逐位核对。
    install(1.,true);dump("saved.csv");
    if(!output_container[0]->SetResultFiles(geo,cfg,sol,0,true))throw std::runtime_error("native output not written");
    if(!boundary||!walls||!source)throw std::runtime_error("control write failed");
  }
};
int main(int argc,char** argv) {
  try {
    SU2_MPI::Init(&argc,&argv);if(argc!=4)throw std::runtime_error("expected config output mode");
    CompatibilityDriver driver(argv[1],1,SU2_MPI::GetComm());driver.check(argv[2],std::string(argv[3])=="restart");
    driver.Finalize();SU2_MPI::Finalize();std::cout<<"RAH_COMPLETE_CONTROL_EXECUTED\n";return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
