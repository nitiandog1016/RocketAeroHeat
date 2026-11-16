// 单壳零孔PLC薄入口。复用TetGen1.5.1原核心；AGPL-3.0-or-later。
// 从本项目fixed_surface_plc.cpp派生，仅增加具名零孔/多marker合同，不修改数值算法。
#include "tetgen.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>

// 模块一：保存真实传入核心的数组，用于独立核对原壁、marker和零孔身份。
static void save_input(const std::string& name,const tetgenio& m) {
  std::ofstream f(name.c_str());f<<std::setprecision(17);
  f<<m.numberofpoints<<' '<<m.numberoffacets<<' '<<m.numberofholes<<'\n';
  for(int i=0;i<m.numberofpoints;++i)f<<m.pointlist[3*i]<<' '<<m.pointlist[3*i+1]<<' '<<m.pointlist[3*i+2]<<' '<<m.pointmarkerlist[i]<<'\n';
  for(int i=0;i<m.numberoffacets;++i){const int* p=m.facetlist[i].polygonlist[0].vertexlist;f<<p[0]<<' '<<p[1]<<' '<<p[2]<<' '<<m.facetmarkerlist[i]<<'\n';}
  if(!f)throw std::runtime_error("PLC_INPUT_SNAPSHOT_FAILED");
}
static void save_output(const std::string& name,const tetgenio& m) {
  std::ofstream f(name.c_str());f<<std::setprecision(17);
  if(!m.pointmarkerlist||!m.trifacemarkerlist||m.numberofcorners!=4)throw std::runtime_error("PLC_OUTPUT_FIELDS");
  f<<m.numberofpoints<<' '<<m.numberoftetrahedra<<' '<<m.numberoftrifaces<<'\n';
  for(int i=0;i<m.numberofpoints;++i)f<<m.pointlist[3*i]<<' '<<m.pointlist[3*i+1]<<' '<<m.pointlist[3*i+2]<<' '<<m.pointmarkerlist[i]<<'\n';
  for(int i=0;i<m.numberoftetrahedra;++i)for(int j=0;j<4;++j)f<<m.tetrahedronlist[4*i+j]-m.firstnumber<<(j==3?'\n':' ');
  for(int i=0;i<m.numberoftrifaces;++i)f<<m.trifacelist[3*i]-m.firstnumber<<' '<<m.trifacelist[3*i+1]-m.firstnumber<<' '<<m.trifacelist[3*i+2]-m.firstnumber<<' '<<m.trifacemarkerlist[i]<<'\n';
  if(!f)throw std::runtime_error("PLC_OUTPUT_FAILED");
}
int main(int argc,char** argv) {
  try {
    if(argc!=3)throw std::runtime_error("usage: single_shell_plc input.dat output_prefix");
    tetgenio in,out;in.firstnumber=0;std::ifstream f(argv[1]);f>>in.numberofpoints>>in.numberoffacets>>in.numberofholes;
    if(!f||in.numberofpoints<4||in.numberofpoints>400000||in.numberoffacets<4||in.numberoffacets>900000||in.numberofholes!=0)
      throw std::runtime_error("SINGLE_SHELL_ZERO_HOLE_COUNTS");
    in.pointlist=new REAL[3*in.numberofpoints];in.pointmarkerlist=new int[in.numberofpoints];
    in.facetlist=new tetgenio::facet[in.numberoffacets];in.facetmarkerlist=new int[in.numberoffacets];
    for(int i=0;i<in.numberoffacets;++i)tetgenio::init(&in.facetlist[i]);
    for(int i=0;i<in.numberofpoints;++i){
      f>>in.pointlist[3*i]>>in.pointlist[3*i+1]>>in.pointlist[3*i+2]>>in.pointmarkerlist[i];
      if(in.pointmarkerlist[i]!=i+1)throw std::runtime_error("PLC_POINT_IDENTITY");
      for(int j=0;j<3;++j)if(!std::isfinite(in.pointlist[3*i+j]))throw std::runtime_error("PLC_POINT_NONFINITE");
    }
    // 模块二：每个三角形独立facet；1原壁、2出口、3远场，无隐式表面合并。
    for(int i=0;i<in.numberoffacets;++i){
      tetgenio::facet& a=in.facetlist[i];a.numberofpolygons=1;a.polygonlist=new tetgenio::polygon[1];
      tetgenio::init(a.polygonlist);a.polygonlist[0].numberofvertices=3;a.polygonlist[0].vertexlist=new int[3];
      for(int j=0;j<3;++j){int& v=a.polygonlist[0].vertexlist[j];f>>v;if(v<0||v>=in.numberofpoints)throw std::runtime_error("PLC_FACE_INDEX");}
      f>>in.facetmarkerlist[i];if(in.facetmarkerlist[i]<1||in.facetmarkerlist[i]>3)throw std::runtime_error("PLC_FACE_MARKER");
    }
    if(!f)throw std::runtime_error("PLC_INCOMPLETE_INPUT");std::string extra;if(f>>extra)throw std::runtime_error("PLC_TRAILING_INPUT");
    const std::string prefix(argv[2]);save_input(prefix+"_import.dat",in);
    // 模块三：固定成熟开关，一次实际调用；不接受CLI参数搜索。
    char sw[]="pY1q2/5O2/7S30000zVC";tetgenbehavior b;
    if(!b.parse_commandline(sw)||!b.plc||b.refine||!b.nobisect||b.nobisect_nomerge!=1||b.convex||b.noexact||
       b.steinerleft!=30000||b.optlevel!=2||b.optscheme!=7||b.minratio!=2||b.mindihedral!=5)
      throw std::runtime_error("PLC_FIXED_SWITCHES");
    std::ofstream meta((prefix+"_options.json").c_str());
    meta<<"{\"switches\":\""<<sw<<"\",\"steinerLimit\":30000,\"inputTetrahedra\":0,\"explicitHoles\":0,\"markers\":[1,2,3]}\n";
    meta.close();if(!meta)throw std::runtime_error("PLC_OPTIONS_WRITE");
    tetrahedralize(&b,&in,&out);save_output(prefix+"_result.dat",out);return 0;
  }catch(int code){std::cerr<<"TetGen exception "<<code<<std::endl;return 10+code;}
   catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 2;}
}
