// 模块职责：完整原生两层FAS/细层SA控制；按实际层身份保存状态和层间算子。
#include "drivers/CSinglezoneDriver.hpp"
#include "iteration/CFluidIteration.hpp"
#include "variables/CNSVariable.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include "SingleDomainHooksV235.hpp"
#include "FASHooksV245.hpp"
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdint>

namespace obs245 {
std::string folder,phase="construction"; bool enabled=false;
std::array<CGeometry*,2> geo{}; std::array<CSolver**,2> sol{}; CConfig* cfg=nullptr;
std::ofstream index; unsigned long seq=0,completed=0;
unsigned level(CGeometry* g){for(unsigned l=0;l<2;++l)if(geo[l]==g)return l;throw std::runtime_error("UNKNOWN_LEVEL");}
// 列定义保持明确：当前守恒量与缓存原始量分列，不能把更新后量和更新前缓存混比。
const char* columns="id,volume,rho,mx,my,mz,Ein,SA,muTsa,T,u,v,w,p,rhoP,Hin,a,mu,muTflow,k,Cp,secondary0,secondary1,secondary2,secondary3,secondary4,secondary5,secondary6,secondary7,old0,old1,old2,old3,old4,rhs0,rhs1,rhs2,rhs3,rhs4,tau0,tau1,tau2,tau3,tau4,wallDistance";
void snap(const std::string& label,unsigned l){
  phase=label;auto* g=geo[l];auto** s=sol[l];auto* n=s[FLOW_SOL]->GetNodes();auto* t=s[TURB_SOL]->GetNodes();
  std::ostringstream name;name<<"frame_"<<std::setw(5)<<std::setfill('0')<<seq<<".bin";
  std::ofstream f(folder+"/"+name.str(),std::ios::binary);uint64_t N=g->GetnPoint(),K=45;
  f.write("RAH245LE",8);f.write(reinterpret_cast<char*>(&N),8);f.write(reinterpret_cast<char*>(&K),8);
  for(unsigned long i=0;i<N;++i){std::array<double,45> r{};unsigned j=0;r[j++]=i;r[j++]=g->nodes->GetVolume(i);
    for(int d=0;d<5;++d)r[j++]=n->GetSolution(i,d);
    r[j++]=t->GetSolution(i,0);r[j++]=t->GetmuT(i);
    for(int d=0;d<12;++d)r[j++]=n->GetPrimitive(i,d);
    for(int d=0;d<8;++d)r[j++]=n->GetSecondary(i,d);
    for(int d=0;d<5;++d)r[j++]=n->GetSolution_Old(i,d);
    for(int d=0;d<5;++d)r[j++]=s[FLOW_SOL]->LinSysRes(i,d);
    for(int d=0;d<5;++d)r[j++]=n->GetResTruncError(i)[d];
    r[j++]=g->nodes->GetWall_Distance(i);if(j!=K)throw std::runtime_error("SCHEMA_COUNT");
    // 仅检查当前守恒状态；缓存的primitive可能属于更新前，后验按阶段区别。
    double v2=0;for(int d=3;d<6;++d)v2+=(r[d]/r[2])*(r[d]/r[2]);
    rah_phase_air::state("rhoe",r[2],r[6]/r[2]-.5*v2);
    for(int d=2;d<9;++d)if(!std::isfinite(r[d]))throw std::runtime_error("NONFINITE_CURRENT_STATE");
    if(label=="prolong_in")for(int d=34;d<39;++d)if(!std::isfinite(r[d]))throw std::runtime_error("NONFINITE_CORRECTION_REJECTED");
    f.write(reinterpret_cast<char*>(r.data()),sizeof(r));
  }
  f.close();if(!f)throw std::runtime_error("SNAP_WRITE");
  index<<seq++<<','<<completed<<','<<cfg->GetInnerIter()<<','<<label<<','<<l<<','<<N<<','<<name.str()<<'\n';index.flush();
}
void stage(const char* label,CGeometry* g,CSolver**,CConfig*,unsigned short sys){
  snap(std::string(cfg->GetContainerPosition(sys)==FLOW_SOL?"flow_":"SA_")+label,level(g));
}
void transfer(const char* label,CGeometry* f,CSolver**,CGeometry* c,CSolver**,CConfig*){
  snap(label,level(f));snap(label,level(c));
}
void reject(CNSVariable* v,unsigned long i,const double*,const double*,bool,const char* reason){
  for(unsigned l=0;l<2;++l)if(sol[l][FLOW_SOL]->GetNodes()==v){
    std::ofstream f(folder+"/rejection.txt");f<<"level="<<l<<" localID="<<i<<" phase="<<phase<<" reason="<<reason<<'\n';f.close();
    throw std::runtime_error("FAS_STATE_REJECTION");
  }
  throw std::runtime_error("FAS_UNKNOWN_REJECTION_LEVEL");
}
void end(){++completed;for(unsigned l=0;l<2;++l)snap("end",l);}
void failure(){std::string message="unknown";try{auto e=std::current_exception();if(e)std::rethrow_exception(e);}
  catch(const std::exception& e){message=e.what();}catch(...){}
  std::ofstream f(folder+"/failure.txt");f<<"phase="<<phase<<"\ncompleted="<<completed<<"\nwhat="<<message<<'\n';f.close();std::_Exit(86);
}
}
// 模块二：父类完整执行flow与SA；不替代原生迭代或重新计算缓存。
class FASIterationV245 final:public CFluidIteration {
 public:using CFluidIteration::CFluidIteration;
 void Iterate(COutput* o,CIntegration**** in,CGeometry**** g,CSolver***** s,CNumerics****** n,CConfig** c,
 CSurfaceMovement** sm,CVolumetricMovement*** gm,CFreeFormDefBox*** b,unsigned short z,unsigned short inst)override{
   CFluidIteration::Iterate(o,in,g,s,n,c,sm,gm,b,z,inst);obs245::end();
 }
};
class FASDriverV245 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void launch(){using namespace obs245;cfg=config_container[0];
   if(cfg->GetnMGLevels()!=1||cfg->GetFinestMesh()!=0||cfg->GetKind_FluidModel()!=RAH_PHASE_AIR||
      cfg->GetKind_Solver()!=MAIN_SOLVER::RANS||cfg->GetKind_Turb_Model()!=TURB_MODEL::SA||
      cfg->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!cfg->GetSAParsedOptions().ft2||omp_get_max_threads()!=1)
     throw std::runtime_error("FAS_IDENTITY_OR_LEVEL_COUNT");
   index.open(folder+"/frames.csv");index<<"seq,completed,inner,label,level,points,file\n";
   std::ofstream links(folder+"/children.csv"),walls(folder+"/walls.csv");links<<"coarse,fine,volumeFine,volumeCoarse\n";
   walls<<"level,point,neighbor,nx,ny,nz,x,y,z,neighborX,neighborY,neighborZ\n";links<<std::setprecision(17);walls<<std::setprecision(17);
   for(unsigned l=0;l<2;++l){geo[l]=geometry_container[0][0][l];sol[l]=solver_container[0][0][l];
     if(!sol[l][TURB_SOL]||sol[l][FLOW_SOL]->GetnVar()!=5||!dynamic_cast<CSU2PhaseAwareAir*>(sol[l][FLOW_SOL]->GetFluidModel()))throw std::runtime_error("LAYER_MODEL_IDENTITY");
     for(unsigned m=0;m<cfg->GetnMarker_All();++m)if(cfg->GetViscous_Wall(m))for(unsigned v=0;v<geo[l]->nVertex[m];++v){
       auto i=geo[l]->vertex[m][v]->GetNode(),k=geo[l]->vertex[m][v]->GetNormal_Neighbor();walls<<l<<','<<i<<','<<k;
       for(int d=0;d<3;++d)walls<<','<<geo[l]->vertex[m][v]->GetNormal()[d];
       for(int d=0;d<3;++d)walls<<','<<geo[l]->nodes->GetCoord(i,d);
       for(int d=0;d<3;++d)walls<<','<<geo[l]->nodes->GetCoord(k,d);walls<<'\n';
     }
   }
   if(geo[1]->GetnPoint()>=geo[0]->GetnPoint())throw std::runtime_error("NO_ACTUAL_COARSENING");
   for(unsigned long c=0;c<geo[1]->GetnPoint();++c)for(unsigned j=0;j<geo[1]->nodes->GetnChildren_CV(c);++j){
     auto f=geo[1]->nodes->GetChildren_CV(c,j);links<<c<<','<<f<<','<<geo[0]->nodes->GetVolume(f)<<','<<geo[1]->nodes->GetVolume(c)<<'\n';
   }
   links.close();walls.close();
   std::ofstream schema(folder+"/schema.json");schema<<"{\"columns\":\""<<columns<<"\",\"offsetJkg\":300000,\"physicalQualification\":false}\n";schema.close();
   delete iteration_container[0][0];iteration_container[0][0]=new FASIterationV245(cfg);
   for(unsigned l=0;l<2;++l)snap("initial",l);
   if(enabled){rah_single_observer::stage=stage;rah_single_observer::rejection=reject;rah_fas_observer::transfer=transfer;}
   StartSolver();rah_single_observer::stage=nullptr;rah_single_observer::rejection=nullptr;rah_fas_observer::transfer=nullptr;
   std::ofstream done(folder+"/native_completion.json");done<<"{\"normalReturn\":true,\"updates\":"<<completed<<",\"finePoints\":"<<geo[0]->GetnPoint()<<",\"coarsePoints\":"<<geo[1]->GetnPoint()<<",\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{
 if(argc!=4)throw std::runtime_error("CONFIG_FOLDER_MODE_REQUIRED");obs245::folder=argv[2];std::string mode=argv[3];
 if(mode!="off"&&mode!="observe")throw std::runtime_error("MODE");obs245::enabled=mode=="observe";
 std::set_terminate(obs245::failure);SU2_MPI::Init(&argc,&argv);FASDriverV245 driver(argv[1],1,SU2_MPI::GetComm());driver.launch();driver.Finalize();SU2_MPI::Finalize();return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
