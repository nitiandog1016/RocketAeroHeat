"""模块用途。为稀疏同阶段输出构建一次独立运行包，复用已审核构建器。"""
from buildSU2IdealAirRuntime import build, ROOT


if __name__ == '__main__':
    build(ROOT/'outputs/v156_validation/build_reproduction', ROOT/'outputs/v156_validation/phase_patch',
          'rah-ideal-air-phase-v2', ROOT/'docs/validation/v156/RECOVERY_PLAN.md',
          ROOT/'scripts/prepareSU2IdealAirPhasePatch.py',
          ('SU2_CFD/src/variables/CNSVariable.cpp', 'SU2_CFD/src/output/CFlowCompOutput.cpp',
           'SU2_CFD/src/output/COutput.cpp', 'SU2_CFD/src/output/CFlowOutput.cpp'))
