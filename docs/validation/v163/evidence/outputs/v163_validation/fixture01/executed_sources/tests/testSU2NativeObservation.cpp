// 模块用途。独立验证普通/OpenMP异常记录和实际C++双精度文本逐位往返。
#include "../third_party/su2_thermal/NativeTermination.hpp"
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  rah_termination::install();
  if (argc != 2) return 2;
  const std::string mode = argv[1];
  if (mode == "std") throw std::domain_error("rah-fixture-domain-error");
  if (mode == "nonstd") throw 17;
  if (mode == "omp") {
    #pragma omp parallel for num_threads(2)
    for (int i = 0; i < 2; ++i) {
      if (i == 1) throw std::domain_error("rah-fixture-openmp-error");
    }
    return 3;
  }
  if (mode != "precision") return 2;
  int old_exact = 0, new_exact = 0;
  size_t maximum_bytes = 0;
  const double center = 2.999999692909394e+02;
  for (int offset = -4; offset <= 4; ++offset) {
    double value = center;
    for (int j = 0; j < std::abs(offset); ++j)
      value = std::nextafter(value, offset < 0 ? -INFINITY : INFINITY);
    for (int precision : {15, std::numeric_limits<double>::max_digits10 - 1}) {
      std::ostringstream data;
      data.precision(precision);
      data << std::scientific << value;
      const double recovered = std::stod(data.str());
      const bool exact = std::memcmp(&value, &recovered, sizeof(double)) == 0;
      if (precision == 15) old_exact += exact;
      else new_exact += exact;
    }
  }
  // 同时核对极端有限双精度的符号、指数和分隔符容量。
  for (double value : {std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(),
                       std::numeric_limits<double>::denorm_min(), -0.0, center}) {
    std::ostringstream data;
    data.precision(std::numeric_limits<double>::max_digits10 - 1);
    data << ", " << std::scientific << value;
    maximum_bytes = (std::max)(maximum_bytes, data.str().size());
  }
  std::cout << "{\"oldExact\":" << old_exact << ",\"newExact\":" << new_exact
            << ",\"samples\":9,\"maximumScalarBytes\":" << maximum_bytes << "}\n";
  return new_exact == 9 && old_exact < 9 && maximum_bytes <= 32 ? 0 : 1;
}
