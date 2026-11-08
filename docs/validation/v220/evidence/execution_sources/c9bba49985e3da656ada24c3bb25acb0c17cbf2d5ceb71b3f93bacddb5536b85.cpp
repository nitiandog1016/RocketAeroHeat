// 模块职责：通过真实CFluidModel/CNSVariable/general Roe检查新气体接入。
// 本控制不推进流场；实际热壁与混合特征边界仍须另作完整小域执行。
#include "CSU2PhaseAwareAir.hpp"
#include "numerics/flow/convection/roe.hpp"
#include "variables/CNSVariable.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

double error(double a,double b){return std::abs(a-b)/std::max(1.,std::abs(b));}
int main(int argc,char** argv) {
 try {
  if(argc!=2)throw std::runtime_error("expected output directory");
  std::istringstream text("SOLVER= RANS\nKIND_TURB_MODEL= SA\nSA_OPTIONS= (NEGATIVE, WITHFT2)\n"
    "SYSTEM_MEASUREMENTS= SI\nREF_DIMENSIONALIZATION= DIMENSIONAL\nFLUID_MODEL= RAH_PHASE_AIR\n"
    "VISCOSITY_MODEL= RAH_PHASE_LINEAR\nCONDUCTIVITY_MODEL= RAH_PHASE_LINEAR\n"
    "CONV_NUM_METHOD_FLOW= ROE\nTIME_DISCRE_FLOW= EULER_IMPLICIT\nTIME_DISCRE_TURB= EULER_IMPLICIT\n"
    "INIT_OPTION= TD_CONDITIONS\nMACH_NUMBER= 5\nFREESTREAM_PRESSURE= 1197.0279007383137\nFREESTREAM_TEMPERATURE= 226.50908361133003\n"
    "OUTPUT_FILES= (RESTART_ASCII)\n");
  CConfig config(text,SU2_COMPONENT::SU2_CFD,false);
  if(config.GetKind_FluidModel()!=RAH_PHASE_AIR || config.GetSAParsedOptions().version!=SA_OPTIONS::NEG ||
     !config.GetSAParsedOptions().ft2)throw std::runtime_error("wrong configuration dispatch");
  CSU2PhaseAwareAir air;CFluidModel& base=air;
  base.SetLaminarViscosityModel(&config);base.SetThermalConductivityModel(&config);
  std::ofstream out(std::string(argv[1])+"/states.csv");out<<std::setprecision(17);
  out<<"case,mode,p,T,rho,e,h,s,cp,cv,a,mu,k,dpdrho,dpde,dtde,dhdrho,dhdp,dsdrho,dsdp,dmu,dk\n";
  unsigned count=0,rejections=0,restorations=0;double recovery=0,primitive_error=0;
  for(double t:{60.,75.123,100.123,132.6312,180.123,226.50908361133003,298.15,300.123,1000.123,1600.123,2000.}) {
   const double p=.05*rah_phase_air::ceiling(t);base.SetTDState_PT(p,t);
   const double rho=base.GetDensity(),e=base.GetStaticEnergy(),h=base.GetEnthalpy(),s=base.GetEntropy();
   for(int mode=0;mode<8;++mode) {
    if(mode==0)base.SetTDState_PT(p,t);
    if(mode==1)base.SetTDState_rhoT(rho,t);
    if(mode==2)base.SetTDState_Prho(p,rho);
    if(mode==3)base.SetTDState_rhoe(rho,e);
    if(mode==4)base.SetTDState_hs(h,s);
    if(mode==5)base.SetTDState_Ps(p,s);
    if(mode==6)base.SetEnergy_Prho(p,rho);
    if(mode==7)base.ComputeDerivativeNRBC_Prho(p,rho);
    recovery=std::max(recovery,std::abs(base.GetTemperature()-t));
    out<<count<<','<<mode<<','<<base.GetPressure()<<','<<base.GetTemperature()<<','<<base.GetDensity()<<','
       <<base.GetStaticEnergy()<<','<<base.GetEnthalpy()<<','<<base.GetEntropy()<<','<<base.GetCp()<<','<<base.GetCv()<<','
       <<base.GetSoundSpeed()<<','<<base.GetLaminarViscosity()<<','<<base.GetThermalConductivity()<<','
       <<base.GetdPdrho_e()<<','<<base.GetdPde_rho()<<','<<base.GetdTde_rho()<<','<<base.Getdhdrho_P()<<','
       <<base.GetdhdP_rho()<<','<<base.Getdsdrho_P()<<','<<base.GetdsdP_rho()<<','<<base.GetdmudT_rho()<<','<<base.GetdktdT_rho()<<'\n';
   }
   // 在内部温度验证真实节点原始量恢复，避免把端点舍入当成新的物理余量。
   if(t>61&&t<1999) {
    double velocity[3]={-900.,10.,300.};double E=e+.5*(900.*900.+100.+300.*300.);
    CNSVariable node(rho,velocity,E,1,3,5,&config);
    if(!node.SetPrimVar(0,0.,0.,&air))throw std::runtime_error("native valid state rejected");
    primitive_error=std::max({primitive_error,error(node.GetTemperature(0),t),error(node.GetPressure(0),p)});
    std::array<double,5> old{},bad{};
    for(int j=0;j<5;++j)old[j]=node.GetSolution(0,j);
    bad=old;bad[4]=-1e9;node.SetSolution_Old(0,old.data());node.SetSolution(0,bad.data());
    if(node.SetPrimVar(0,0.,0.,&air))throw std::runtime_error("invalid candidate was accepted");
    for(int j=0;j<5;++j)if(node.GetSolution(0,j)!=old[j])throw std::runtime_error("old state was not restored");
    primitive_error=std::max({primitive_error,error(node.GetTemperature(0),t),error(node.GetPressure(0),p)});
    ++restorations;
   }
   ++count;
  }
  // 模块二：在固定内点调用实际setter作中心差分；不跨节点，不调扰动寻找通过。
  std::ofstream derivatives(std::string(argv[1])+"/derivatives.csv");derivatives<<std::setprecision(17);
  derivatives<<"T,name,analytic,finiteDifference\n";
  double derivative_error=0;
  for(double t:{75.123,180.123,300.123,1000.123,1600.123}) {
    const double p=.05*rah_phase_air::ceiling(t);base.SetTDState_PT(p,t);
    const double r=base.GetDensity(),e=base.GetStaticEnergy();
    const double dr=r*1e-5,de=std::max(1.,std::abs(e))*1e-6,dp=p*1e-5,dt=t*1e-6;
    auto check=[&](const char* name,double analytic,auto setter,auto getter,double x,double step) {
      setter(x+step);double plus=getter();setter(x-step);double minus=getter();
      const double fd=(plus-minus)/(2*step);derivative_error=std::max(derivative_error,error(fd,analytic));
      derivatives<<t<<','<<name<<','<<analytic<<','<<fd<<'\n';base.SetTDState_PT(p,t);
    };
    check("dpdrho_e",base.GetdPdrho_e(),[&](double x){base.SetTDState_rhoe(x,e);},[&]{return base.GetPressure();},r,dr);
    check("dpde_rho",base.GetdPde_rho(),[&](double x){base.SetTDState_rhoe(r,x);},[&]{return base.GetPressure();},e,de);
    check("dTdrho_e",base.GetdTdrho_e(),[&](double x){base.SetTDState_rhoe(x,e);},[&]{return base.GetTemperature();},r,dr);
    check("dTde_rho",base.GetdTde_rho(),[&](double x){base.SetTDState_rhoe(r,x);},[&]{return base.GetTemperature();},e,de);
    check("dhdrho_p",base.Getdhdrho_P(),[&](double x){base.SetTDState_Prho(p,x);},[&]{return base.GetEnthalpy();},r,dr);
    check("dhdp_rho",base.GetdhdP_rho(),[&](double x){base.SetTDState_Prho(x,r);},[&]{return base.GetEnthalpy();},p,dp);
    check("dsdrho_p",base.Getdsdrho_P(),[&](double x){base.SetTDState_Prho(p,x);},[&]{return base.GetEntropy();},r,dr);
    check("dsdp_rho",base.GetdsdP_rho(),[&](double x){base.SetTDState_Prho(x,r);},[&]{return base.GetEntropy();},p,dp);
    check("dmudT_rho",base.GetdmudT_rho(),[&](double x){base.SetTDState_rhoT(r,x);},[&]{return base.GetLaminarViscosity();},t,dt);
    check("dkdT_rho",base.GetdktdT_rho(),[&](double x){base.SetTDState_rhoT(r,x);},[&]{return base.GetThermalConductivity();},t,dt);
    check("dmudrho_T",base.Getdmudrho_T(),[&](double x){base.SetTDState_rhoT(x,t);},[&]{return base.GetLaminarViscosity();},r,dr);
    check("dkdrho_T",base.Getdktdrho_T(),[&](double x){base.SetTDState_rhoT(x,t);},[&]{return base.GetThermalConductivity();},r,dr);
  }
  base.SetTDState_PT(1000,300);
  const double before_e=base.GetStaticEnergy(),before_mu=base.GetLaminarViscosity();
  auto reject=[&](auto fn){try{fn();}catch(const std::domain_error&){++rejections;return;}throw std::runtime_error("invalid state accepted");};
  reject([&]{base.SetTDState_PT(1000,59.9);});reject([&]{base.SetTDState_PT(1,2000.1);});
  reject([&]{base.SetTDState_PT(rah_phase_air::ceiling(60)*1.01,60);});
  reject([&]{base.SetTDState_rhoe(-1,before_e);});reject([&]{base.SetTDState_Prho(1000,-1);});
  reject([&]{base.SetTDState_rhoe(.01,-1e9);});reject([&]{base.SetTDState_hs(1e10,0);});
  reject([&]{base.SetTDState_Ps(1000,1e9);});reject([&]{base.SetTDState_PT(std::numeric_limits<double>::quiet_NaN(),300);});
  if(base.GetTemperature()!=300||base.GetStaticEnergy()!=before_e||base.GetLaminarViscosity()!=before_mu)
    throw std::runtime_error("rejection was not atomic");
  base.SetEddyViscosity(.05);
  if(base.GetThermalConductivity()!=rah_phase_air::PT(1000,300).k)throw std::runtime_error("double turbulent conductivity");

  // 模块三：实际general Roe在两种显式能量坐标下的通量协变性，不复写Roe。
  CUpwGeneralRoe_Flow roe(3,5,&config);double normal[3]={.3,-.4,.5};roe.SetNormal(normal);
  double vi[12]{},vj[12]{},si[8]{},sj[8]{};
  auto fill=[&](double* v,double* secondary,double pp,double tt,double u,double w) {
   auto q=rah_phase_air::PT(pp,tt);v[0]=tt;v[1]=u;v[2]=12;v[3]=w;v[4]=pp;v[5]=q.rho;
   v[6]=q.h+.5*(u*u+144+w*w);v[7]=std::sqrt(q.a2);v[8]=q.mu;v[10]=q.k;v[11]=q.cp;
   secondary[0]=q.dpdrho;secondary[1]=q.dpde;
  };
  fill(vi,si,1197,226.509,-1418,516);fill(vj,sj,1600,300.123,-1250,450);
  roe.SetPrimitive(vi,vj);roe.SetSecondary(si,sj);
  const auto first=roe.ComputeResidual(&config);std::array<double,5> shifted{};
  for(int j=0;j<5;++j)shifted[j]=first[j];
  vi[6]-=rah_phase_air::shift;vj[6]-=rah_phase_air::shift;
  const auto second=roe.ComputeResidual(&config);double flux_error=0;
  for(int j=0;j<5;++j)flux_error=std::max(flux_error,error(shifted[j]-(j==4?rah_phase_air::shift*shifted[0]:0.),second[j]));
  // 独立制造非法Roe平均，必须明确报错，不能静默返回零通量。
  si[0]=-1e12;sj[0]=-1e12;
  reject([&]{roe.ComputeResidual(&config);});
  if(recovery>1e-7||primitive_error>2e-10||flux_error>2e-10||derivative_error>3e-6||rejections!=10||restorations!=9)
    throw std::runtime_error("native control gate failed");
  std::ofstream summary(std::string(argv[1])+"/native_report.json");summary<<std::setprecision(17)
    <<"{\"passed\":true,\"stateRows\":"<<count*8<<",\"rejections\":"<<rejections
    <<",\"temperatureRecoveryK\":"<<recovery<<",\"nativePrimitiveError\":"<<primitive_error
    <<",\"nativeRoeCovarianceError\":"<<flux_error
    <<",\"derivativeDifferenceError\":"<<derivative_error<<",\"restoredStates\":"<<restorations
    <<",\"actualBoundaryExecutionVerified\":false,\"physicalQualification\":false}\n";
  if(!out||!derivatives||!summary)throw std::runtime_error("output write failed");
  std::cout<<"RAH_PHASE_NATIVE_CONTROL_PASSED\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
