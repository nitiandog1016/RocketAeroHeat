// 模块职责：从原生单区算例文件完成初始化，读回工作对象与节点状态。
#include "drivers/CSinglezoneDriver.hpp"
#include "variables/CNSVariable.hpp"
#include "CaloricFactoryV261.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>

class StaticAirDriver final : public CSinglezoneDriver {
public:
  using CSinglezoneDriver::CSinglezoneDriver;

  void inspect() {
    auto* cfg = config_container[0];
    auto* geometry = geometry_container[0][0][0];
    auto* solver = solver_container[0][0][0][FLOW_SOL];
    auto* nodes = dynamic_cast<CNSVariable*>(solver->GetNodes());
    auto* gas = solver->GetFluidModel();
    if (!nodes || !gas || geometry->GetnDim()!=3 || geometry->GetnPointDomain()!=125 ||
        cfg->GetKind_FluidModel()!=RAH_PHASE_AIR || omp_get_max_threads()!=1)
      throw std::runtime_error("STATIC_AIR_CASE_IDENTITY");
    auto* stable = dynamic_cast<CSU2PhaseAwareAirStableV260*>(gas);
    const char* name = rah_caloric_factory::method(gas);
    const auto startCalls = stable ? stable->calls : 0;

    // 两帧之间只恢复原生原始量，各节点的全部守恒量均需逐位一致。
    for (const auto* filename : {"initial.csv","terminal.csv"}) {
      std::ofstream output(filename);
      output << std::setprecision(17) << "id,x,y,z,rho,rhoUx,rhoUy,rhoUz,rhoE_internal,T_K,p_Pa,method\n";
      for (unsigned long i=0; i<geometry->GetnPointDomain(); ++i) {
        double initial[5];
        for (unsigned j=0; j<5; ++j) initial[j]=nodes->GetSolution(i,j);
        if (!nodes->SetPrimVar(i,0.,0.,gas)) throw std::runtime_error("INITIAL_STATE_REFUSED");
        output << geometry->nodes->GetGlobalIndex(i);
        for (unsigned j=0; j<3; ++j) output << ',' << geometry->nodes->GetCoord(i,j);
        for (unsigned j=0; j<5; ++j) {
          double value=nodes->GetSolution(i,j);
          if (std::memcmp(&value,&initial[j],sizeof(value))) throw std::runtime_error("CHANGED_CONSERVED_STATE");
          output << ',' << value;
        }
        output << ',' << nodes->GetTemperature(i) << ',' << nodes->GetPressure(i) << ',' << name << '\n';
      }
      output.close();
      if (!output) throw std::runtime_error("FRAME_WRITE_FAILED");
    }
    const auto frameCalls = stable ? stable->calls - startCalls : 0;
    std::string inversePath = stable ? stable->lastCaloricInverse.method : "legacy-TOMS748";

    // 在独立节点上检验新类型的原生异常处理。保留拒绝候选和恢复结果供复核。
    double rho=nodes->GetSolution(0,0), rhoE=nodes->GetSolution(0,4), velocity[3]{};
    CNSVariable probe(rho,velocity,rhoE/rho,1,3,5,cfg);
    for (unsigned j=0; j<5; ++j) {
      probe.SetSolution_Old(0,j,nodes->GetSolution(0,j));
      probe.SetSolution(0,j,nodes->GetSolution(0,j));
    }
    probe.SetSolution(0,4,-rho*1e6);
    const bool accepted = probe.SetPrimVar(0,0.,0.,gas);
    bool restored = true;
    for (unsigned j=0; j<5; ++j) restored &= probe.GetSolution(0,j)==nodes->GetSolution(0,j);
    if (accepted || !restored) throw std::runtime_error("DOMAIN_RECOVERY_FAILED");
    std::ofstream report("initialization.json");
    report << std::setprecision(17) << std::boolalpha
           << "{\"workerMethod\":\"" << name << "\",\"inversePath\":\"" << inversePath
           << "\",\"stableRecoveryCallsInFrames\":" << frameCalls << ",\"threadsObserved\":1,"
           << "\"nodes\":125,\"nonlinearUpdates\":0,\"timeUpdates\":0,"
           << "\"rejectedCandidateSpecificInternalEnergy\":-1000000,"
           << "\"invalidCandidateAccepted\":" << accepted << ",\"oldConservedStateRestored\":" << restored
           << ",\"restoredTemperatureK\":" << probe.GetTemperature(0) << "}\n";
    report.close();
    if (!report) throw std::runtime_error("RECEIPT_WRITE_FAILED");
  }
};

int main(int argc,char** argv) {
  if (argc!=2) return 2;
  try {
    SU2_MPI::Init(&argc,&argv);
    StaticAirDriver driver(argv[1],1,SU2_MPI::GetComm());
    driver.inspect();
    driver.Finalize();
    SU2_MPI::Finalize();
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
