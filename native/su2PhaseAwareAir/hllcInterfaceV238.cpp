// 模块职责：调用成熟一般EOS HLLC对象，保存通量/Jacobian给独立解析审核，不自写通量公式。
#include "numerics/flow/convection/hllc.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
using Vector=std::array<double,5>;
struct Primitive {double p,t;std::array<double,3> v;};
struct Pair {std::string name;Primitive left,right;bool exactLeft;};
struct Result {Vector f;std::array<double,25> jl,jr;};

Result flux(CUpwGeneralHLLC_Flow& h,CConfig& cfg,const Primitive& l,const Primitive& r,double offset,double sign){
  double a[12]{},b[12]{},sa[8]{},sb[8]{},normal[3]={sign,0,0};
  auto fill=[&](const Primitive& p,double* v,double* s){
    const auto q=rah_phase_air::PT(p.p,p.t);double v2=0;v[0]=q.t;
    for(int d=0;d<3;++d){v[d+1]=p.v[d];v2+=p.v[d]*p.v[d];}
    v[4]=q.p;v[5]=q.rho;v[6]=q.h-300000+offset+.5*v2;v[7]=std::sqrt(q.a2);
    v[8]=q.mu;v[10]=q.k;v[11]=q.cp;s[0]=q.dpdrho;s[1]=q.dpde;
  };
  fill(l,a,sa);fill(r,b,sb);h.SetNormal(normal);h.SetPrimitive(a,b);h.SetSecondary(sa,sb);
  const auto result=h.ComputeResidual(&cfg);Result out{};
  for(int i=0;i<5;++i){out.f[i]=result[i];for(int j=0;j<5;++j){out.jl[i*5+j]=result.jacobian_i[i][j];out.jr[i*5+j]=result.jacobian_j[i][j];}}
  for(double x:out.f)if(!std::isfinite(x))throw std::runtime_error("HLLC_NONFINITE_FLUX");
  for(const auto* m:{&out.jl,&out.jr})for(double x:*m)if(!std::isfinite(x))throw std::runtime_error("HLLC_NONFINITE_JACOBIAN");
  return out;
}

int main(int argc,char** argv){try{
 if(argc!=3)throw std::runtime_error("CONFIG_OUTPUT_REQUIRED");SU2_MPI::Init(&argc,&argv);
 CConfig cfg(argv[1],SU2_COMPONENT::SU2_CFD,false);CUpwGeneralHLLC_Flow hllc(3,5,&cfg);
 if(cfg.GetKind_FluidModel()!=RAH_PHASE_AIR||cfg.GetKind_TimeIntScheme_Flow()!=EULER_IMPLICIT)throw std::runtime_error("HLLC_CONTROL_IDENTITY");
 Primitive fast{1200,350.123,{1200,40,-20}},plus{1200,350.123,{120,40,-20}},minus{1200,350.123,{-120,40,-20}};
 std::vector<Pair> cases={{"uniform_supersonic",fast,fast,true},{"uniform_subsonic_plus",plus,plus,true},{"uniform_subsonic_minus",minus,minus,true},
   {"stationary_contact",{1200,300,{0,20,-10}},{1200,600,{0,20,-10}},true},
   {"moving_contact",{1200,300,{120,20,-10}},{1200,600,{120,20,-10}},true},
   {"expansion",{1200,350,{-50,20,-10}},{800,300,{50,20,-10}},false}};
 std::ofstream f(std::string(argv[2])+"/flux.csv"),reject(std::string(argv[2])+"/rejected.csv");
 f<<std::setprecision(17)<<"case,offset,reverse,exactLeft,pL,TL,uL,vL,wL,pR,TR,uR,vR,wR";
 for(int j=0;j<5;++j)f<<",F"<<j;for(auto name:{"JL","JR"})for(int i=0;i<5;++i)for(int j=0;j<5;++j)f<<','<<name<<i<<j;f<<'\n';
 unsigned calls=0;
 for(const auto& c:cases)for(double offset:{0.,300000.})for(bool reverse:{false,true}){
   const auto r=flux(hllc,cfg,reverse?c.right:c.left,reverse?c.left:c.right,offset,reverse?-1.:1.);
   f<<c.name<<','<<offset<<','<<reverse<<','<<c.exactLeft;
   for(const Primitive* p:std::array<const Primitive*,2>{&c.left,&c.right}){f<<','<<p->p<<','<<p->t;for(double v:p->v)f<<','<<v;}
   for(double x:r.f)f<<','<<x;for(const auto* m:{&r.jl,&r.jr})for(double x:*m)f<<','<<x;f<<'\n';++calls;
 }
 reject<<"case,error\n";unsigned rejects=0;
 for(const Primitive p:std::vector<Primitive>{{1200,59,{0,0,0}},{-1,300,{0,0,0}},{1e7,300,{0,0,0}}}){
   try{flux(hllc,cfg,p,plus,300000,1);}catch(const std::domain_error& e){reject<<++rejects<<','<<e.what()<<'\n';continue;}
   throw std::runtime_error("HLLC_INVALID_INPUT_ACCEPTED");
 }
 std::ofstream report(std::string(argv[2])+"/report.json");report<<"{\"nativeCalls\":"<<calls<<",\"inputRejections\":"<<rejects<<",\"primitiveStatesFromSamePhaseAir\":true,\"actualFactoryDispatchTested\":false,\"scientificAcceptancePendingIndependentAudit\":true,\"flowUpdates\":0,\"physicalQualification\":false}\n";
 if(!f||!reject||!report)throw std::runtime_error("HLLC_CONTROL_OUTPUT");SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
