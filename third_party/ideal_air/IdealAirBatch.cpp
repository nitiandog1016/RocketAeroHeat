// 模块用途。将已验证的热力学核提供给Python批量初场构造，全部量采用SI。
#include "IdealAir.hpp"
#include <cstdio>
#include <cstddef>

#ifdef _WIN32
#define RAH_EXPORT __declspec(dllexport)
#else
#define RAH_EXPORT __attribute__((visibility("default")))
#endif

// 每行依次输出T、rho、p、e、h、cp和声速。异常只返回失败位置，调用方丢弃整批。
extern "C" RAH_EXPORT std::size_t rah_ideal_air_batch(
    int mode, std::size_t count, const double* first, const double* second,
    double* output, char* error, std::size_t capacity) {
  for (std::size_t row = 0; row < count; ++row) {
    try {
      rah_air::State state;
      if (mode == 0) {
        state = rah_air::from_pressure_temperature(first[row], second[row]);
      } else if (mode == 1) {
        const double temperature = rah_air::temperature_from(
            second[row], rah_air::InverseQuantity::Enthalpy);
        state = rah_air::from_pressure_temperature(first[row], temperature);
      } else if (mode == 2) {
        state = rah_air::state("rhoe", first[row], second[row]);
      } else {
        throw std::invalid_argument("unsupported batch state pair");
      }
      const double values[] = {state.temperature, state.rho, state.pressure,
          state.energy, state.enthalpy, state.cp, state.sound_speed};
      std::copy(values, values + 7, output + 7 * row);
    } catch (const std::exception& exception) {
      std::snprintf(error, capacity, "%s", exception.what());
      return row + 1;
    }
  }
  return 0;
}
