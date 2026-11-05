// 独立原生物性参照和有限性能检查，链接现有固定CoolProp静态库。
#include "../third_party/ideal_air/IdealAir.hpp"
#include "AbstractState.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

using Clock = std::chrono::steady_clock;

int main(int argc, char** argv) {
  try {
    if (argc != 4) throw std::runtime_error("usage: probe --core|--oracle|--benchmark input output");
    std::unique_ptr<CoolProp::AbstractState> fluid(CoolProp::AbstractState::factory("HEOS", "Air"));
    fluid->specify_phase(CoolProp::iphase_gas);
    std::ifstream input(argv[2]);
    std::ofstream output(argv[3]);
    if (!input || !output) throw std::runtime_error("cannot open input/output");
    output << std::setprecision(17);
    const std::string kind = argv[1];
    if (kind == "--benchmark") {
      std::vector<rah_air::State> states;
      std::vector<double> heos_energies;
      for (unsigned i = 0; i < 1000; ++i) {
        states.push_back(rah_air::from_rho_temperature(.02, 80.+1.9*i));
        fluid->update(CoolProp::DmassT_INPUTS, .02, 80.+1.9*i);
        heos_energies.push_back(fluid->umass());
      }
      volatile double checksum = 0;
      output << "repeat,analytic_rhoT_s,analytic_rhoe_s,HEOS_rhoT_s,HEOS_rhoe_s,checksum\n";
      for (unsigned repeat = 0; repeat < 4; ++repeat) {
        const auto start = Clock::now();
        for (unsigned i = 0; i < 20000; ++i) {
          const auto& source = states[i%states.size()];
          auto value = rah_air::from_rho_temperature(source.rho, source.temperature);
          checksum += value.energy + value.cp + value.entropy;
        }
        const auto middle = Clock::now();
        for (unsigned i = 0; i < 20000; ++i) {
          const auto& source = states[i%states.size()];
          auto value = rah_air::state("rhoe", source.rho, source.energy);
          checksum += value.energy + value.cp + value.entropy;
        }
        const auto second = Clock::now();
        for (unsigned i = 0; i < 20000; ++i) {
          const auto& source = states[i%states.size()];
          fluid->update(CoolProp::DmassT_INPUTS, source.rho, source.temperature);
          checksum += fluid->umass() + fluid->cpmass() + fluid->smass();
        }
        const auto end = Clock::now();
        for (unsigned i = 0; i < 20000; ++i) {
          const auto& source = states[i%states.size()];
          fluid->update(CoolProp::DmassUmass_INPUTS, source.rho, heos_energies[i%states.size()]);
          checksum += fluid->umass() + fluid->cpmass() + fluid->smass();
        }
        const auto inverse_end = Clock::now();
        output << repeat << ',' << std::chrono::duration<double>(middle-start).count() << ','
               << std::chrono::duration<double>(second-middle).count() << ','
               << std::chrono::duration<double>(end-second).count() << ','
               << std::chrono::duration<double>(inverse_end-end).count() << ',' << checksum << '\n';
      }
      return 0;
    }
    if (kind != "--core" && kind != "--oracle") throw std::runtime_error("unsupported probe mode");
    std::string line;
    std::getline(input, line);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != "mode,a,b") throw std::runtime_error("invalid input header");
    output << "row,accepted,T,rho,p,e,h,s,cp,cv,sound,dp_drho_e,dp_de_rho,dt_drho_e,dt_de_rho,"
              "dh_drho_p,dh_dp_rho,ds_drho_p,ds_dp_rho,dcp_dt,mu,lambda,dmu_dt,dlambda_dt\n";
    unsigned row = 0;
    while (std::getline(input, line)) {
      if (line.empty()) continue;
      std::replace(line.begin(), line.end(), ',', ' ');
      std::istringstream stream(line);
      std::string mode, first_text, second_text, extra;
      if (!(stream >> mode >> first_text >> second_text) || (stream >> extra)) throw std::runtime_error("invalid input row");
      const double first = std::stod(first_text);
      const double second = std::stod(second_text);
      try {
        rah_air::State value{};
        if (kind == "--core") value = rah_air::state(mode, first, second);
        else {
          if (mode != "rhoT") throw std::domain_error("oracle accepts rhoT only");
          fluid->update(CoolProp::DmassT_INPUTS, first, second);
          value.rho = first;
          value.temperature = second;
          value.pressure = first*fluid->gas_constant()/fluid->molar_mass()*second;
          value.energy = fluid->umass_idealgas();
          value.enthalpy = fluid->hmass_idealgas();
          value.entropy = fluid->smass_idealgas();
          value.cp = fluid->cp0mass();
          value.cv = value.cp-fluid->gas_constant()/fluid->molar_mass();
          value.sound_speed = std::sqrt(value.cp/value.cv*value.pressure/value.rho);
        }
        const auto transport = rah_air::transport(value);
        output << row << ",1," << value.temperature << ',' << value.rho << ',' << value.pressure << ','
               << value.energy << ',' << value.enthalpy << ',' << value.entropy << ',' << value.cp << ','
               << value.cv << ',' << value.sound_speed << ',' << value.dp_drho_e << ',' << value.dp_de_rho << ','
               << value.dt_drho_e << ',' << value.dt_de_rho << ',' << value.dh_drho_p << ',' << value.dh_dp_rho << ','
               << value.ds_drho_p << ',' << value.ds_dp_rho << ',' << value.dcp_dt << ',' << transport.viscosity << ','
               << transport.conductivity << ',' << transport.dviscosity_dt << ',' << transport.dconductivity_dt << '\n';
      } catch (const std::exception& error) {
        output << row << ",0";
        for (unsigned i = 0; i < 22; ++i) output << ",nan";
        output << '\n';
        std::cerr << "rejected row " << row << ": " << error.what() << '\n';
      }
      ++row;
    }
    if (!row || !output) throw std::runtime_error("empty or failed output");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
