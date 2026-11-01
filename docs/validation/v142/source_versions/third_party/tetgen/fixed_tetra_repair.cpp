// TetGen 1.5.1 的独立研究入口：导入已有带孔网格并保留原始输入/输出。
// 上游数值源码不修改；本封装与链接引擎按 AGPL-3.0-or-later 使用。
#include "tetgen.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

// 模块1：输出完整数组；零编号，17位有效数字保留double往返。
static void save(const std::string& path, const tetgenio& m) {
  std::ofstream f(path.c_str());
  if (!f || !m.pointmarkerlist || !m.trifacemarkerlist || m.numberofcorners != 4)
    throw std::runtime_error("missing output fields");
  f << m.numberofpoints << ' ' << m.numberoftetrahedra << ' ' << m.numberoftrifaces << '\n';
  f << std::setprecision(17);
  for (int i=0;i<m.numberofpoints;++i)
    f << m.pointlist[3*i] << ' ' << m.pointlist[3*i+1] << ' ' << m.pointlist[3*i+2]
      << ' ' << m.pointmarkerlist[i] << '\n';
  for (int i=0;i<m.numberoftetrahedra;++i) {
    for(int j=0;j<4;++j) f << m.tetrahedronlist[4*i+j]-m.firstnumber << (j==3?'\n':' ');
  }
  for (int i=0;i<m.numberoftrifaces;++i)
    f << m.trifacelist[3*i]-m.firstnumber << ' ' << m.trifacelist[3*i+1]-m.firstnumber
      << ' ' << m.trifacelist[3*i+2]-m.firstnumber << ' ' << m.trifacemarkerlist[i] << '\n';
  if (!f) throw std::runtime_error("output write failure");
}

int main(int argc,char** argv) {
  try {
    if (argc!=3) throw std::runtime_error("usage: fixed_tetra_repair input.dat output_prefix");
    tetgenio in,out; in.firstnumber=0; in.numberofcorners=4;
    std::ifstream f(argv[1]);
    f >> in.numberofpoints >> in.numberoftetrahedra >> in.numberoftrifaces;
    if (!f || in.numberofpoints<4 || in.numberofpoints>400000 ||
        in.numberoftetrahedra<1 || in.numberoftetrahedra>900000 ||
        in.numberoftrifaces<4 || in.numberoftrifaces>900000)
      throw std::runtime_error("invalid mesh size");
    in.pointlist=new REAL[3*in.numberofpoints]; in.pointmarkerlist=new int[in.numberofpoints];
    in.tetrahedronlist=new int[4*in.numberoftetrahedra];
    in.trifacelist=new int[3*in.numberoftrifaces]; in.trifacemarkerlist=new int[in.numberoftrifaces];
    for(int i=0;i<in.numberofpoints;++i)
      f >> in.pointlist[3*i] >> in.pointlist[3*i+1] >> in.pointlist[3*i+2] >> in.pointmarkerlist[i];
    for(int i=0;i<4*in.numberoftetrahedra;++i) {
      f >> in.tetrahedronlist[i];
      if(in.tetrahedronlist[i]<0 || in.tetrahedronlist[i]>=in.numberofpoints)
        throw std::runtime_error("tetrahedron index out of range");
    }
    for(int i=0;i<in.numberoftrifaces;++i) {
      for(int j=0;j<3;++j) {
        f >> in.trifacelist[3*i+j];
        if(in.trifacelist[3*i+j]<0 || in.trifacelist[3*i+j]>=in.numberofpoints)
          throw std::runtime_error("face index out of range");
      }
      f >> in.trifacemarkerlist[i];
    }
    if(!f) throw std::runtime_error("incomplete mesh input");
    std::string extra; if(f >> extra) throw std::runtime_error("trailing input");
    std::string prefix(argv[2]); save(prefix+"_import.dat",in);

    // 模块2：固定参数；先保存解析结果，禁止调用方静默扫描参数。
    char switches[]="rY1q2/5O2/7S10000zVC"; tetgenbehavior b;
    if(!b.parse_commandline(switches) || !b.refine || !b.nobisect || b.convex ||
       b.nobisect_nomerge!=1 || b.optlevel!=2 || b.optscheme!=7 ||
       b.steinerleft!=10000 || b.minratio!=2 || b.mindihedral!=5 || b.noexact)
      throw std::runtime_error("unexpected parsed behavior");
    std::ofstream meta((prefix+"_options.json").c_str());
    meta << "{\"switches\":\"" << switches << "\",\"optLevel\":" << b.optlevel
      << ",\"optScheme\":" << b.optscheme << ",\"maxDihedral\":" << b.optmaxdihedral
      << ",\"steinerLimit\":" << b.steinerleft << ",\"minRatio\":" << b.minratio
      << ",\"minDihedral\":" << b.mindihedral << "}\n";
    meta.close(); if(!meta) throw std::runtime_error("options write failure");
    std::cout << "FIXED_SWITCHES " << switches << std::endl;
    // 模块3：仅一次原生调用；异常/失败由外部监督器保留。
    tetrahedralize(&b,&in,&out);
    save(prefix+"_result.dat",out);
    return 0;
  } catch(int code) {
    std::cerr << "TetGen exception " << code << std::endl; return 10+code;
  } catch(const std::exception& e) {
    std::cerr << e.what() << std::endl; return 2;
  }
}
