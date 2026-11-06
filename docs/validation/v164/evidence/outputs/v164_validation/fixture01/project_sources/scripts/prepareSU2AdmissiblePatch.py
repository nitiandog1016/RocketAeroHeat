"""模块用途。准确父源上生成限定物性回退补丁，供独立实际类测试和单对象构建复用。"""
from contextlib import contextmanager
from pathlib import Path

from buildSU2IdealAirRuntime import ROOT, SOURCE, BUILD, sha, save
from prepareSU2IdealAirPatch import replace
from su2ThermalRuntime import read

PARENT = ROOT/'outputs/v163_validation/build03'
FILE = 'SU2_CFD/src/variables/CNSVariable.cpp'
HEADER = 'SU2_CFD/include/variables/AdmissibilityTrace.hpp'
FLAVOR = 'rah-ideal-air-admissible-v1'


def patch(text):
    """仅在物性拒绝时触发原生回退，成功态不读取陈旧物性缓存。"""
    text = replace(text, '#include "../../include/fluid/CFluidModel.hpp"',
        '#include "../../include/fluid/CFluidModel.hpp"\n'
        '#include "../../include/fluid/CSU2IdealAir.hpp"\n'
        '#include "../../include/variables/AdmissibilityTrace.hpp"\n#include <array>\n#include <string>')
    start = text.index('bool CNSVariable::SetPrimVar(')
    end = text.index('\nvoid CNSVariable::SetSecondaryVar(', start)
    original = text[start:end]
    changed = replace(original, '''  FluidModel->SetTDState_rhoe(density, staticEnergy);

  bool check_dens  = SetDensity(iPoint);
  bool check_press = SetPressure(iPoint, FluidModel->GetPressure());
  bool check_sos   = SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
  bool check_temp  = SetTemperature(iPoint, FluidModel->GetTemperature());''', '''  // 模块一。新物性拒绝候选态时，进入原生旧态恢复，避免读取旧物性字段。
  bool air_state_rejected = false;
  std::string rejected_reason;
  try {
    FluidModel->SetTDState_rhoe(density, staticEnergy);
  } catch (const std::domain_error& error) {
    if (dynamic_cast<CSU2IdealAir*>(FluidModel) == nullptr) throw;
    air_state_rejected = true;
    rejected_reason = error.what();
  }

  bool check_dens  = air_state_rejected || SetDensity(iPoint);
  bool check_press = air_state_rejected || SetPressure(iPoint, FluidModel->GetPressure());
  bool check_sos   = air_state_rejected || SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
  bool check_temp  = air_state_rejected || SetTemperature(iPoint, FluidModel->GetTemperature());''')
    changed = replace(changed, '    /*--- Copy the old solution ---*/', '''    // 模块二。失败证据使用明确的局部编号，保留原始全部守恒量与当前k。
    std::array<double, 5> trial{}, restored{};
    for (unsigned long iVar = 0; iVar < nVar; ++iVar) {
      trial[iVar] = Solution(iPoint,iVar);
      restored[iVar] = Solution_Old(iPoint,iVar);
    }
    if (nDim == 2) {
      trial[4] = trial[3]; trial[3] = 0.0;
      restored[4] = restored[3]; restored[3] = 0.0;
    }
    const double trial_energy = staticEnergy;

    /*--- Copy the old solution ---*/''')
    changed = replace(changed, '    FluidModel->SetTDState_rhoe(density, staticEnergy);\n\n    SetDensity(iPoint);', '''    // 模块三。旧守恒态必须与本次实际k共同通过完整物性查询，否则明确终止。
    try {
      FluidModel->SetTDState_rhoe(density, staticEnergy);
    } catch (const std::domain_error& error) {
      if (dynamic_cast<CSU2IdealAir*>(FluidModel) != nullptr)
        rah_admissibility::record(iPoint, nDim, trial.data(), restored.data(),
                                  trial_energy, staticEnergy, turb_ke, false, error.what());
      throw;
    }
    if (air_state_rejected)
      rah_admissibility::record(iPoint, nDim, trial.data(), restored.data(),
                                trial_energy, staticEnergy, turb_ke, true, rejected_reason.c_str());

    SetDensity(iPoint);''')
    return text[:start]+changed+text[end:]


@contextmanager
def overlay(folder):
    """临时恢复父包准确源树，保留本次执行副本，并在所有终态恢复共享文件。"""
    old_tree = read(PARENT/'source_tree.json')
    snapshots = {SOURCE/p.relative_to(PARENT/'executed_sources'): p.read_bytes()
                 for p in (PARENT/'executed_sources').rglob('*') if p.is_file()}
    snapshots[SOURCE/FILE] = patch(snapshots[SOURCE/FILE].decode('utf8')).encode('utf8')
    snapshots[SOURCE/HEADER] = (ROOT/'third_party/ideal_air/AdmissibilityTrace.hpp').read_bytes()
    originals = {path: path.read_bytes() if path.exists() else None for path in snapshots}
    try:
        for path, raw in snapshots.items():
            path.write_bytes(raw)
        assert all(sha(Path(name)) == digest for name, digest in old_tree.items() if name != str(SOURCE/FILE))
        tree = {name: sha(Path(name)) for name in old_tree}
        tree[str(SOURCE/HEADER)] = sha(SOURCE/HEADER)
        save(folder/'source_tree.json', tree)
        for name in ('docs/validation/v164/PLAN.md', 'scripts/prepareSU2AdmissiblePatch.py',
                     'third_party/ideal_air/AdmissibilityTrace.hpp'):
            target = folder/'project_sources'/name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((ROOT/name).read_bytes())
        yield snapshots
        assert all(sha(Path(name)) == digest for name, digest in tree.items())
    finally:
        restored = {}
        for path, raw in originals.items():
            target = folder/'executed_sources'/path.relative_to(SOURCE)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(path.read_bytes())
            if raw is None:
                path.resolve().relative_to(SOURCE.resolve())
                path.unlink()
                restored[str(path)] = not path.exists()
            else:
                path.write_bytes(raw)
                restored[str(path)] = path.read_bytes() == raw
        save(folder/'source_restore.json', restored)
        assert all(restored.values())
