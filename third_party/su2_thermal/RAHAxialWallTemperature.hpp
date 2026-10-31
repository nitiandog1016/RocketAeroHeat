// RocketAeroHeat 原生热边界适配器；LGPL-2.1-or-later，与SU2同许可证。
// 仅负责读取/校验SI温度表和记录实际输入；不修改任何流动或热流离散公式。
#ifndef RAH_AXIAL_WALL_TEMPERATURE_HPP
#define RAH_AXIAL_WALL_TEMPERATURE_HPP
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class RAHAxialWallTemperature {
 public:
  // 模块一：没有显式文件时保持官方行为；有文件时任何错误都不能回退到标量。
  RAHAxialWallTemperature() {
    const char* file = std::getenv("RAH_WALL_TEMPERATURE_FILE");
    if (file && *file) Load(file);
  }
  explicit RAHAxialWallTemperature(const std::string& file) { Load(file); }
  bool Active() const { return !x_.empty(); }

  // 模块二：线性插值、全节点覆盖和端点严格检查，不外推、不裁剪温度。
  double At(double x) const {
    if (!Active() || !std::isfinite(x) || x < x_.front() || x > x_.back())
      throw std::runtime_error("RAH thermal table does not cover wall node X");
    const auto it = std::lower_bound(x_.begin(), x_.end(), x);
    if (it == x_.begin()) return t_.front();
    if (it == x_.end()) return t_.back();
    const auto j = static_cast<std::size_t>(it - x_.begin());
    const double w = (x - x_[j-1]) / (x_[j] - x_[j-1]);
    return (1-w)*t_[j-1] + w*t_[j];
  }

 private:
  std::vector<double> x_, t_;
  static double Number(const std::string& token) {
    std::istringstream in(token);in.imbue(std::locale::classic());double value;
    if (!(in >> value) || !std::isfinite(value))
      throw std::runtime_error("RAH thermal table requires finite numeric data");
    in >> std::ws;
    if (!in.eof()) throw std::runtime_error("RAH thermal table has trailing data");
    return value;
  }
  void Load(const std::string& file) {
    std::ifstream in(file);std::string line;
    if (!in || !std::getline(in,line)) throw std::runtime_error("RAH thermal table cannot be read");
    if (!line.empty() && line.back()=='\r') line.pop_back();
    if (line != "RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K")
      throw std::runtime_error("RAH thermal table schema or SI units mismatch");
    while (std::getline(in,line)) {
      if (!line.empty() && line.back()=='\r') line.pop_back();
      const auto comma=line.find(',');
      if (line.size()>4096 || comma==std::string::npos || line.find(',',comma+1)!=std::string::npos)
        throw std::runtime_error("RAH thermal table requires exactly two columns");
      const double x=Number(line.substr(0,comma)),t=Number(line.substr(comma+1));
      if (t<=0 || (!x_.empty() && x<=x_.back()) || x_.size()>=10000)
        throw std::runtime_error("RAH thermal table must have increasing X and positive T");
      x_.push_back(x);t_.push_back(t);
    }
    if (!in.eof() || x_.size()<2) throw std::runtime_error("RAH thermal table is incomplete");
  }
};
#endif
