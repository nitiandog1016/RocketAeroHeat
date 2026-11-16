// 模块职责：真实完整域的原生聚合与初始化预检；本入口没有流动/SA更新循环。
#include "drivers/CSinglezoneDriver.hpp"
#include "integration/CMultiGridIntegration.hpp"
#include "fluid/CSU2PhaseAwareAir.hpp"
#include "SingleDomainHooksV235.hpp"
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <limits>

namespace pf246 {
std::string folder,phase="construction";unsigned long rejections=0;
void reject(CNSVariable*,unsigned long i,const double*,const double*,bool,const char* reason){
  ++rejections;std::ofstream f(folder+"/rejection.txt");f<<phase<<" point="<<i<<" reason="<<reason<<'\n';f.close();
  throw std::runtime_error("INITIAL_STATE_REJECTED");
}
void failed(){std::string message="unknown";try{auto e=std::current_exception();if(e)std::rethrow_exception(e);}
  catch(const std::exception& e){message=e.what();}catch(...){}
  std::ofstream f(folder+"/failure.txt");f<<"phase="<<phase<<"\nwhat="<<message<<'\n';f.close();std::_Exit(86);
}
void frame(CGeometry* g,CSolver** s,const std::string& name){
  phase=name;auto* n=s[FLOW_SOL]->GetNodes();auto* t=s[TURB_SOL]->GetNodes();
  std::ofstream f(folder+"/"+name+".bin",std::ios::binary);uint64_t N=g->GetnPoint(),K=33;
  f.write("RAH246LE",8);f.write(reinterpret_cast<char*>(&N),8);f.write(reinterpret_cast<char*>(&K),8);
  for(unsigned long i=0;i<N;++i){std::array<double,33> r{};int j=0;r[j++]=g->nodes->GetGlobalIndex(i);
    for(int d=0;d<3;++d)r[j++]=g->nodes->GetCoord(i,d);r[j++]=g->nodes->GetVolume(i);
    for(int d=0;d<5;++d)r[j++]=n->GetSolution(i,d);r[j++]=t->GetSolution(i,0);r[j++]=t->GetmuT(i);
    for(int d=0;d<12;++d)r[j++]=n->GetPrimitive(i,d);for(int d=0;d<8;++d)r[j++]=n->GetSecondary(i,d);
    r[j++]=g->nodes->GetWall_Distance(i);if(j!=K)throw std::runtime_error("FRAME_SCHEMA");
    f.write(reinterpret_cast<char*>(r.data()),sizeof(r));
  }
  f.close();if(!f)throw std::runtime_error("FRAME_WRITE");
}
}

class FASPreflightV246 final:public CSinglezoneDriver {
 public:using CSinglezoneDriver::CSinglezoneDriver;
 void inspect(){using namespace pf246;auto* c=config_container[0];auto** g=geometry_container[0][0];auto*** s=solver_container[0][0];
   if(c->GetnMGLevels()!=1||c->GetFinestMesh()!=0||c->GetKind_FluidModel()!=RAH_PHASE_AIR||
      c->GetKind_Solver()!=MAIN_SOLVER::RANS||c->GetKind_Turb_Model()!=TURB_MODEL::SA||c->GetRestart()||
      c->GetSAParsedOptions().version!=SA_OPTIONS::NEG||!c->GetSAParsedOptions().ft2||omp_get_max_threads()!=1||
      g[0]->GetnPoint()!=335726||g[0]->GetnElem()!=705490||g[1]->GetnPoint()>=g[0]->GetnPoint())
     throw std::runtime_error("ORIGINAL_DOMAIN_OR_FAS_IDENTITY");
   for(int l=0;l<2;++l)if(!s[l][TURB_SOL]||!dynamic_cast<CSU2PhaseAwareAir*>(s[l][FLOW_SOL]->GetFluidModel()))throw std::runtime_error("LAYER_IDENTITY");
   std::ofstream schema(folder+"/schema.json");schema<<"{\"columns\":\"globalID,x,y,z,volume,rho,mx,my,mz,Ein,SA,muTsa,T,u,v,w,p,rhoP,Hin,a,mu,muTflow,k,Cp,secondary0,secondary1,secondary2,secondary3,secondary4,secondary5,secondary6,secondary7,wallDistance\",\"offsetJkg\":300000,\"flowUpdates\":0,\"SAUpdates\":0,\"physicalQualification\":false}\n";schema.close();
   phase="parent-and-boundary-observation";
   // 模块一：保留真实局部/全球号和全部边界，不因分区或粗化丢弃任何原壁。
   std::ofstream links(folder+"/children.bin",std::ios::binary);uint64_t count=g[0]->GetnPoint(),cols=2;
   links.write("RAH246ID",8);links.write(reinterpret_cast<char*>(&count),8);links.write(reinterpret_cast<char*>(&cols),8);
   unsigned long actual=0;for(unsigned long p=0;p<g[1]->GetnPoint();++p)for(unsigned j=0;j<g[1]->nodes->GetnChildren_CV(p);++j){
     std::array<uint64_t,2> row{p,g[1]->nodes->GetChildren_CV(p,j)};links.write(reinterpret_cast<char*>(row.data()),sizeof(row));++actual;}
   links.close();if(!links||actual!=count)throw std::runtime_error("CHILDREN_COUNT");
   std::ofstream boundary(folder+"/boundary.csv");boundary<<std::setprecision(17);
   boundary<<"level,marker,kind,point,globalID,neighbor,nx,ny,nz,x,y,z,neighborX,neighborY,neighborZ,isViscous\n";
   for(unsigned l=0;l<2;++l)for(unsigned m=0;m<c->GetnMarker_All();++m)for(unsigned v=0;v<g[l]->nVertex[m];++v){
     auto i=g[l]->vertex[m][v]->GetNode(),j=g[l]->vertex[m][v]->GetNormal_Neighbor();
     boundary<<l<<','<<m<<','<<c->GetMarker_All_KindBC(m)<<','<<i<<','<<g[l]->nodes->GetGlobalIndex(i)<<','<<j;
     for(int d=0;d<3;++d)boundary<<','<<g[l]->vertex[m][v]->GetNormal()[d];
     for(int d=0;d<3;++d)boundary<<','<<g[l]->nodes->GetCoord(i,d);
     for(int d=0;d<3;++d)boundary<<','<<g[l]->nodes->GetCoord(j,d);boundary<<','<<c->GetViscous_Wall(m)<<'\n';
   }
   boundary.close();if(!boundary)throw std::runtime_error("BOUNDARY_WRITE");
   std::ofstream walls(folder+"/fine_walls.csv");walls<<"marker,face,nodes,g0,g1,g2,g3\n";unsigned long facecount=0;
   for(unsigned m=0;m<c->GetnMarker_All();++m)if(c->GetViscous_Wall(m))for(unsigned long e=0;e<g[0]->GetnElem_Bound(m);++e){
     auto* face=g[0]->bound[m][e];const auto n=face->GetnNodes();walls<<m<<','<<e<<','<<n;
     for(unsigned j=0;j<4;++j)walls<<','<<(j<n?static_cast<long long>(g[0]->nodes->GetGlobalIndex(face->GetNode(j))):-1ll);walls<<'\n';++facecount;}
   walls.close();if(facecount!=19752||!walls)throw std::runtime_error("ORIGINAL_WALL_COUNT");
   frame(g[0],s[0],"fine_constructed");frame(g[1],s[1],"coarse_constructed");
   // 模块二：只执行原生初始化和flow层间限制；不调用私有SA限制替身，不推进任何方程。
   rah_single_observer::rejection=reject;phase="native-driver-preprocess";Preprocess(0);
   phase="fine-primitive-refresh";s[0][FLOW_SOL]->Preprocessing(g[0],s[0],c,0,NO_RK_ITER,RUNTIME_FLOW_SYS,true);
   frame(g[0],s[0],"fine_initialized");
   auto* mg=dynamic_cast<CMultiGridIntegration*>(integration_container[0][0][FLOW_SOL]);if(!mg)throw std::runtime_error("NATIVE_MG_TYPE");
   phase="native-flow-restriction";mg->SetRestricted_Solution(RUNTIME_FLOW_SYS,s[0][FLOW_SOL],s[1][FLOW_SOL],g[0],g[1],c);
   phase="coarse-primitive-refresh";s[1][FLOW_SOL]->Preprocessing(g[1],s[1],c,1,NO_RK_ITER,RUNTIME_FLOW_SYS,true);
   frame(g[1],s[1],"coarse_restricted");rah_single_observer::rejection=nullptr;
   std::ofstream f(folder+"/native_completion.json");f<<"{\"normalReturn\":true,\"finePoints\":"<<g[0]->GetnPoint()<<",\"coarsePoints\":"<<g[1]->GetnPoint()<<",\"fineCells\":"<<g[0]->GetnElem()<<",\"wallFaces\":"<<facecount<<",\"flowUpdates\":0,\"SAUpdates\":0,\"rejections\":"<<rejections<<",\"physicalQualification\":false}\n";
 }
};
int main(int argc,char** argv){try{
  if(argc!=3)throw std::runtime_error("CONFIG_FOLDER_REQUIRED");pf246::folder=argv[2];std::set_terminate(pf246::failed);
  SU2_MPI::Init(&argc,&argv);FASPreflightV246 driver(argv[1],1,SU2_MPI::GetComm());driver.inspect();driver.Finalize();SU2_MPI::Finalize();return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
