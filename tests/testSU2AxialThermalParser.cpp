// 独立输入测试：正常插值、常量极限与缺失/非法/超覆盖拒绝，不是物理精度测试。
#include "../third_party/su2_thermal/RAHAxialWallTemperature.hpp"
#include <cassert>
#include <iostream>
static void Write(const std::string& name,const std::string& body) {
  std::ofstream out(name);out<<body;out.close();assert(out);
}
int main() {
  const std::string head="RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K\n";
  Write("linear.csv",head+"0,200\n0.2,280\n0.5,310\n");
  RAHAxialWallTemperature p("linear.csv");
  assert(p.At(0)==200 && p.At(.5)==310 && std::abs(p.At(.1)-240)<1e-12);
  assert(std::abs(p.At(.35)-295)<1e-12);
  Write("constant.csv",head+"0,300\n1,300\n");
  assert(RAHAxialWallTemperature("constant.csv").At(.378)==300);
  int rejected=0;
  for(const std::string& text:{head+"0,200\n0,210\n",head+"1,200\n0,210\n",
      head+"0,-1\n1,210\n",head+"0,nan\n1,210\n",head+"0,200junk\n1,210\n",
      head+"0,200,1\n1,210\n",head+"0,200\n",std::string("bad-units\n0,200\n1,210\n")}) {
    Write("invalid.csv",text);
    try { RAHAxialWallTemperature bad("invalid.csv"); }
    catch(const std::runtime_error&) { ++rejected; }
  }
  for(double x:{-.001,.501}) {
    try { p.At(x); } catch(const std::runtime_error&) { ++rejected; }
  }
  try { RAHAxialWallTemperature missing("does_not_exist.csv"); }
  catch(const std::runtime_error&) { ++rejected; }
  assert(rejected==11);
  std::cout<<"Thermal parser interpolation, constant limit and 11 rejection cases passed.\n";
}
