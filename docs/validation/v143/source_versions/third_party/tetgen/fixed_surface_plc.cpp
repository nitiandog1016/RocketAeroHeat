// 固定三角PLC/显式空腔入口，原TetGen核心不改。AGPL-3.0-or-later。
#include "tetgen.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// 模块1：保留实际表面输入，验证时不依赖输入文件回读冒充原生导入。
static void save_input(const std::string& name,const tetgenio& m) {
  std::ofstream f(name.c_str()); f << std::setprecision(17);
  f << m.numberofpoints << ' ' << m.numberoffacets << ' ' << m.numberofholes << '\n';
  for(int i=0;i<m.numberofpoints;++i)
    f << m.pointlist[3*i] << ' ' << m.pointlist[3*i+1] << ' ' << m.pointlist[3*i+2] << ' ' << m.pointmarkerlist[i] << '\n';
  for(int i=0;i<m.numberoffacets;++i) {
    const int* p=m.facetlist[i].polygonlist[0].vertexlist;
    f << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << m.facetmarkerlist[i] << '\n';
  }
  for(int i=0;i<m.numberofholes;++i) f << m.holelist[3*i] << ' ' << m.holelist[3*i+1] << ' ' << m.holelist[3*i+2] << '\n';
  if(!f) throw std::runtime_error("input snapshot write failure");
}
static void save_output(const std::string& name,const tetgenio& m) {
  std::ofstream f(name.c_str()); f << std::setprecision(17);
  if(!m.pointmarkerlist || !m.trifacemarkerlist || m.numberofcorners!=4) throw std::runtime_error("missing output fields");
  f << m.numberofpoints << ' ' << m.numberoftetrahedra << ' ' << m.numberoftrifaces << '\n';
  for(int i=0;i<m.numberofpoints;++i)
    f << m.pointlist[3*i] << ' ' << m.pointlist[3*i+1] << ' ' << m.pointlist[3*i+2] << ' ' << m.pointmarkerlist[i] << '\n';
  for(int i=0;i<m.numberoftetrahedra;++i)
    for(int j=0;j<4;++j) f << m.tetrahedronlist[4*i+j]-m.firstnumber << (j==3?'\n':' ');
  for(int i=0;i<m.numberoftrifaces;++i)
    f << m.trifacelist[3*i]-m.firstnumber << ' ' << m.trifacelist[3*i+1]-m.firstnumber << ' ' << m.trifacelist[3*i+2]-m.firstnumber << ' ' << m.trifacemarkerlist[i] << '\n';
  if(!f) throw std::runtime_error("result write failure");
}
int main(int argc,char** argv) {
  try {
    if(argc!=3) throw std::runtime_error("usage: fixed_surface_plc input.dat output_prefix");
    tetgenio in,out; in.firstnumber=0;
    std::ifstream f(argv[1]); f >> in.numberofpoints >> in.numberoffacets >> in.numberofholes;
    if(!f || in.numberofpoints<4 || in.numberofpoints>400000 || in.numberoffacets<4 || in.numberoffacets>900000 || in.numberofholes!=1)
      throw std::runtime_error("invalid PLC counts");
    in.pointlist=new REAL[3*in.numberofpoints]; in.pointmarkerlist=new int[in.numberofpoints];
    in.facetlist=new tetgenio::facet[in.numberoffacets]; in.facetmarkerlist=new int[in.numberoffacets];
    for(int i=0;i<in.numberoffacets;++i) tetgenio::init(&in.facetlist[i]);
    in.holelist=new REAL[3*in.numberofholes];
    for(int i=0;i<in.numberofpoints;++i) {
      f >> in.pointlist[3*i] >> in.pointlist[3*i+1] >> in.pointlist[3*i+2] >> in.pointmarkerlist[i];
      if(in.pointmarkerlist[i]!=i+1) throw std::runtime_error("boundary identity invalid");
    }
    // 模块2：每个输入三角面单独一个facet，禁止隐式改变三角剖分。
    for(int i=0;i<in.numberoffacets;++i) {
      tetgenio::facet& a=in.facetlist[i]; a.numberofpolygons=1; a.polygonlist=new tetgenio::polygon[1];
      tetgenio::init(a.polygonlist); a.polygonlist[0].numberofvertices=3; a.polygonlist[0].vertexlist=new int[3];
      for(int j=0;j<3;++j) {
        int& v=a.polygonlist[0].vertexlist[j]; f >> v;
        if(v<0 || v>=in.numberofpoints) throw std::runtime_error("facet index invalid");
      }
      f >> in.facetmarkerlist[i]; if(in.facetmarkerlist[i]!=1) throw std::runtime_error("facet marker invalid");
    }
    for(int i=0;i<3*in.numberofholes;++i) f >> in.holelist[i];
    if(!f) throw std::runtime_error("incomplete PLC");
    std::string extra; if(f >> extra) throw std::runtime_error("trailing PLC input");
    std::string prefix(argv[2]); save_input(prefix+"_import.dat",in);
    // 模块3：固定参数和单次调用；不通过命令行扫描网格设置。
    char sw[]="pY1q2/5O2/7S30000zVC"; tetgenbehavior b;
    if(!b.parse_commandline(sw) || !b.plc || b.refine || !b.nobisect || b.nobisect_nomerge!=1 ||
       b.convex || b.noexact || b.steinerleft!=30000 || b.optlevel!=2 || b.optscheme!=7 || b.minratio!=2 || b.mindihedral!=5)
      throw std::runtime_error("unexpected parsed options");
    std::ofstream meta((prefix+"_options.json").c_str());
    meta << "{\"switches\":\"" << sw << "\",\"steinerLimit\":30000,\"inputTetrahedra\":0,\"explicitHoles\":1}\n";
    meta.close(); if(!meta) throw std::runtime_error("options write failure");
    std::cout << "FIXED_SWITCHES " << sw << std::endl;
    tetrahedralize(&b,&in,&out); save_output(prefix+"_result.dat",out); return 0;
  } catch(int code) { std::cerr << "TetGen exception " << code << std::endl; return 10+code;
  } catch(const std::exception& e) { std::cerr << e.what() << std::endl; return 2; }
}
