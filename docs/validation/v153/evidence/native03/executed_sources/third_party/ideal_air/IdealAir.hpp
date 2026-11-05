// 固定组分理想空气热力学核。SI单位，显式60—2000K范围。
// 系数来源为Lemmon2000式(24)及固定CoolProp7.2.0 Air的OTH参考偏移。
// 本核只提供热力学和明确选择的工程输运，流场动能/SST k由调用方分离。
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace rah_air {
constexpr double gas_constant = 8.31451 / 0.02896546;
constexpr double reducing_temperature = 132.6312;
constexpr double reducing_density = 10447.7 * 0.02896546;
constexpr double minimum_temperature = 60.0;
constexpr double maximum_temperature = 2000.0;

struct Potential {
  double value, first, second, third;
};

struct State {
  double rho, temperature, pressure, energy, enthalpy, entropy, cp, cv, sound_speed;
  double dp_drho_e, dp_de_rho, dt_drho_e, dt_de_rho;
  double dh_drho_p, dh_dp_rho, ds_drho_p, ds_dp_rho, dcp_dt;
};

inline void require(bool condition, const char* message) {
  if (!condition) throw std::domain_error(message);
}

// 模块一，温度Helmholtz势及解析导数。密度项ln(delta)在组装状态时加入。
inline Potential potential(double temperature) {
  require(std::isfinite(temperature) && temperature >= minimum_temperature &&
          temperature <= maximum_temperature, "temperature outside 60..2000 K");
  const double tau = reducing_temperature / temperature;
  Potential result{10.3753039487406 + 3.31112445645577*tau,
                   3.31112445645577, 0.0, 0.0};
  constexpr std::array<double, 6> coefficients{{6.057194e-8, -2.10274769e-5,
      -0.000158860716, -13.841928076, 17.275266575, -0.00019536342}};
  constexpr std::array<double, 6> exponents{{-3.0, -2.0, -1.0, 0.0, 1.0, 1.5}};
  const double inverse_tau = 1/tau;
  const std::array<double, 6> powers{{inverse_tau*inverse_tau*inverse_tau,
      inverse_tau*inverse_tau, inverse_tau, 1.0, tau, tau*std::sqrt(tau)}};
  for (size_t i = 0; i < coefficients.size(); ++i) {
    const double exponent = exponents[i];
    const double term = coefficients[i]*powers[i];
    result.value += term;
    result.first += term*exponent/tau;
    result.second += term*exponent*(exponent-1)/(tau*tau);
    result.third += term*exponent*(exponent-1)*(exponent-2)/(tau*tau*tau);
  }
  constexpr double log_coefficient = 2.490888032;
  result.value += log_coefficient*std::log(tau);
  result.first += log_coefficient/tau;
  result.second -= log_coefficient/(tau*tau);
  result.third += 2*log_coefficient/(tau*tau*tau);
  constexpr std::array<double, 2> weights{{0.791309509, 0.212236768}};
  constexpr std::array<double, 2> frequencies{{25.36365, 16.90741}};
  for (size_t i = 0; i < weights.size(); ++i) {
    const double frequency = frequencies[i];
    const double decay = std::exp(-frequency*tau);
    const double complement = -std::expm1(-frequency*tau);
    result.value += weights[i]*std::log(complement);
    result.first += weights[i]*frequency*decay/complement;
    result.second -= weights[i]*frequency*frequency*decay/(complement*complement);
    result.third += weights[i]*frequency*frequency*frequency*decay*(1+decay)/
                    (complement*complement*complement);
  }
  constexpr double weight = -0.197938904;
  constexpr double frequency = 87.31279;
  const double scaled_decay = (2.0/3.0)*std::exp(-frequency*tau);
  const double denominator = 1+scaled_decay;
  result.value += weight*(frequency*tau + std::log1p(scaled_decay));
  result.first += weight*frequency/denominator;
  result.second += weight*frequency*frequency*scaled_decay/(denominator*denominator);
  result.third -= weight*frequency*frequency*frequency*scaled_decay*(1-scaled_decay)/
                  (denominator*denominator*denominator);
  return result;
}

// 模块二，一次组装完整状态，异常时调用者不会收到部分更新后的状态对象。
inline State from_rho_temperature(double rho, double temperature) {
  require(std::isfinite(rho) && rho > 0, "density must be finite and positive");
  const Potential alpha = potential(temperature);
  const double tau = reducing_temperature/temperature;
  State value{};
  value.rho = rho;
  value.temperature = temperature;
  value.pressure = rho*gas_constant*temperature;
  value.energy = gas_constant*reducing_temperature*alpha.first;
  value.enthalpy = value.energy + gas_constant*temperature;
  value.entropy = gas_constant*(tau*alpha.first-alpha.value-std::log(rho/reducing_density));
  value.cv = -gas_constant*tau*tau*alpha.second;
  value.cp = value.cv + gas_constant;
  value.sound_speed = std::sqrt(value.cp/value.cv*gas_constant*temperature);
  value.dp_drho_e = gas_constant*temperature;
  value.dp_de_rho = rho*gas_constant/value.cv;
  value.dt_drho_e = 0;
  value.dt_de_rho = 1/value.cv;
  value.dh_drho_p = -value.cp*temperature/rho;
  value.dh_dp_rho = value.cp/(rho*gas_constant);
  value.ds_drho_p = -value.cp/rho;
  value.ds_dp_rho = value.cv/value.pressure;
  value.dcp_dt = gas_constant/temperature*(2*tau*tau*alpha.second+tau*tau*tau*alpha.third);
  require(value.cv > 0 && value.pressure > 0 && std::isfinite(value.pressure) &&
          std::isfinite(value.entropy) && std::isfinite(value.dp_de_rho) &&
          std::isfinite(value.dh_drho_p) && std::isfinite(value.dh_dp_rho) &&
          std::isfinite(value.ds_dp_rho), "state exceeds numerical range");
  return value;
}

inline State from_pressure_temperature(double pressure, double temperature) {
  require(std::isfinite(pressure) && pressure > 0, "pressure must be finite and positive");
  return from_rho_temperature(pressure/(gas_constant*temperature), temperature);
}

// 模块三，单调热内能/焓/等压熵反解，区间外目标拒绝，最多60次有界迭代。
enum class InverseQuantity { Energy, Enthalpy, EntropyAtPressure };

inline double temperature_from(double target, InverseQuantity quantity, double pressure = 101325.0) {
  require(std::isfinite(target), "inverse target must be finite");
  require(std::isfinite(pressure) && pressure > 0, "pressure must be finite and positive");
  auto evaluate = [=](double temperature) {
    const Potential alpha = potential(temperature);
    const double tau = reducing_temperature/temperature;
    const double energy = gas_constant*reducing_temperature*alpha.first;
    const double cv = -gas_constant*tau*tau*alpha.second;
    if (quantity == InverseQuantity::Energy) return std::array<double, 2>{{energy, cv}};
    if (quantity == InverseQuantity::Enthalpy)
      return std::array<double, 2>{{energy+gas_constant*temperature, cv+gas_constant}};
    const double entropy = gas_constant*(tau*alpha.first-alpha.value-
        std::log(pressure/(gas_constant*temperature*reducing_density)));
    return std::array<double, 2>{{entropy, (cv+gas_constant)/temperature}};
  };
  double lower = minimum_temperature;
  double upper = maximum_temperature;
  // 单一模型的端点预计算只读缓存，熵按目标压力作解析修正。
  static const State low = from_pressure_temperature(101325.0, minimum_temperature);
  static const State high = from_pressure_temperature(101325.0, maximum_temperature);
  double low_value = low.energy;
  double high_value = high.energy;
  if (quantity == InverseQuantity::Enthalpy) {
    low_value = low.enthalpy;
    high_value = high.enthalpy;
  } else if (quantity == InverseQuantity::EntropyAtPressure) {
    // 直接使用同一熵表达式，避免压力修正改变端点最后一位。
    low_value = evaluate(lower)[0];
    high_value = evaluate(upper)[0];
  }
  require(target >= low_value && target <= high_value, "inverse target outside temperature range");
  if (target == low_value) return lower;
  if (target == high_value) return upper;
  double temperature = lower+(upper-lower)*(target-low_value)/(high_value-low_value);
  for (unsigned iteration = 0; iteration < 60; ++iteration) {
    const auto value = evaluate(temperature);
    const double residual = value[0]-target;
    if (std::abs(residual) <= 2e-13*std::max(std::abs(target), 1.0)) return temperature;
    if (residual < 0) lower = temperature;
    else upper = temperature;
    const double candidate = temperature-residual/value[1];
    temperature = candidate > lower && candidate < upper ? candidate : (lower+upper)/2;
  }
  throw std::runtime_error("temperature inversion failed");
}

inline State state(const std::string& mode, double first, double second) {
  require(std::isfinite(first) && std::isfinite(second), "nonfinite state input");
  if (mode == "rhoT") return from_rho_temperature(first, second);
  if (mode == "PT") return from_pressure_temperature(first, second);
  if (mode == "rhoe") return from_rho_temperature(first, temperature_from(second, InverseQuantity::Energy));
  if (mode == "Prho") {
    require(first > 0 && second > 0, "pressure and density must be positive");
    double temperature = first/(second*gas_constant);
    const double tolerance = 8*std::numeric_limits<double>::epsilon()*maximum_temperature;
    require(temperature >= minimum_temperature-tolerance &&
            temperature <= maximum_temperature+tolerance, "derived temperature outside 60..2000 K");
    temperature = std::max(minimum_temperature, std::min(maximum_temperature, temperature));
    return from_rho_temperature(second, temperature);
  }
  if (mode == "hs") {
    const double temperature = temperature_from(first, InverseQuantity::Enthalpy);
    const State reference = from_rho_temperature(reducing_density, temperature);
    const double rho = reducing_density*std::exp((reference.entropy-second)/gas_constant);
    return from_rho_temperature(rho, temperature);
  }
  if (mode == "Ps") return from_pressure_temperature(first, temperature_from(second, InverseQuantity::EntropyAtPressure, first));
  throw std::invalid_argument("unsupported state pair");
}

struct Transport {
  double viscosity, conductivity, dviscosity_dt, dconductivity_dt;
};

// 模块四，显式工程输运。湍黏度偏导以给定mu_t冻结，lambda导数含dcp/dT。
inline Transport transport(const State& state, double eddy_viscosity = 0.0) {
  require(std::isfinite(eddy_viscosity) && eddy_viscosity >= 0, "invalid eddy viscosity");
  const double temperature = state.temperature;
  const double viscosity = 1.716e-5*std::pow(temperature/273.15, 1.5)*383.55/(temperature+110.4);
  const double derivative = viscosity*(1.5/temperature-1/(temperature+110.4));
  const double factor = viscosity/.72+eddy_viscosity/.9;
  const double conductivity = state.cp*factor;
  const double conductivity_derivative = state.dcp_dt*factor+state.cp*derivative/.72;
  require(std::isfinite(conductivity) && std::isfinite(conductivity_derivative), "transport exceeds numerical range");
  return {viscosity, conductivity, derivative, conductivity_derivative};
}
}  // namespace rah_air
