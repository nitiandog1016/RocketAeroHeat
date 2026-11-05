"""模块用途。复用物性补丁，增加与原始热状态绑定的k字段和首更新保存。"""
from pathlib import Path
import difflib
import hashlib
import json
import zipfile

from prepareSU2IdealAirPatch import BASE, ROOT, patches, replace

OUT = ROOT/'outputs/v156_validation/patch01'


def phase_patches(original):
    changed = patches(original)
    name = 'SU2_CFD/include/variables/CNSVariable.hpp'
    text = replace(original[name], '  VectorType Tau_Wall;',
        '  VectorType FlowPrimitiveK; /*!< 对应最近一次原始热状态计算的独立k输入。 */\n\n  VectorType Tau_Wall;')
    changed[name] = replace(text, 'public:\n',
        'public:\n  su2double GetFlowPrimitiveK(unsigned long iPoint) const { return FlowPrimitiveK(iPoint); }\n\n')
    name = 'SU2_CFD/src/variables/CNSVariable.cpp'
    text = replace(original[name], '  Tau_Wall.resize(nPoint)',
        '  FlowPrimitiveK.resize(nPoint) = su2double(0.0);\n  Tau_Wall.resize(nPoint)')
    changed[name] = replace(text,
        'bool CNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {',
        'bool CNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {\n\n'
        '  // 直接捕获热状态计算输入，包含原生回退后的同一次调用。\n'
        '  FlowPrimitiveK(iPoint) = turb_ke;')
    name = 'SU2_CFD/src/output/CFlowCompOutput.cpp'
    text = replace(original[name], '#include "../../include/solvers/CSolver.hpp"',
        '#include "../../include/solvers/CSolver.hpp"\n#include "../../include/variables/CNSVariable.hpp"')
    text = replace(text, '  AddCommonFVMOutputs(config);',
        '  if (config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP)\n'
        '    AddVolumeOutput("FLOW_PRIMITIVE_K", "Flow_Primitive_K", "FLOW_PHASE", "Independent k used by saved flow primitive state");\n\n'
        '  AddCommonFVMOutputs(config);')
    changed[name] = replace(text, '  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);',
        '  if (config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP)\n'
        '    SetVolumeOutputValue("FLOW_PRIMITIVE_K", iPoint, static_cast<const CNSVariable*>(Node_Flow)->GetFlowPrimitiveK(iPoint));\n\n'
        '  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);')
    name = 'SU2_CFD/src/output/COutput.cpp'
    changed[name] = replace(original[name],
        '      return ((Iter > 0) && (Iter % config->GetVolumeOutputFrequency(iFile) == 0)) || force_writing;',
        '      // 新模型首个完成更新也保存，后续沿用原生局部索引和强制末次规则。\n'
        '      return (((Iter > 0 || config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP) &&\n'
        '                (Iter % config->GetVolumeOutputFrequency(iFile) == 0))) || force_writing;')
    return changed


def main():
    assert not OUT.exists()
    names = ['Common/include/option_structure.hpp', 'Common/src/CConfig.cpp',
        'SU2_CFD/src/solvers/CEulerSolver.cpp', 'SU2_CFD/include/solvers/CEulerSolver.hpp',
        'SU2_CFD/src/solvers/CNSSolver.cpp', 'SU2_CFD/include/variables/CNSVariable.hpp',
        'SU2_CFD/src/variables/CNSVariable.cpp', 'SU2_CFD/src/output/CFlowCompOutput.cpp',
        'SU2_CFD/src/output/COutput.cpp']
    with zipfile.ZipFile(str(BASE/'su2.zip')) as bundle:
        original = {name: bundle.read(next(p for p in bundle.namelist() if p.endswith('/'+name))).decode('utf8')
                    for name in names}
    changed = phase_patches(original)
    OUT.mkdir(parents=True)
    diff = []
    manifest = {}
    for name, text in changed.items():
        target = OUT/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(text.encode('utf8'))
        manifest[name] = dict(originalSHA256=hashlib.sha256(original[name].encode()).hexdigest(),
                             patchedSHA256=hashlib.sha256(text.encode()).hexdigest())
        diff.extend(difflib.unified_diff(original[name].splitlines(True), text.splitlines(True),
                                       fromfile='a/'+name, tofile='b/'+name))
    (OUT/'ideal_air.diff').write_bytes(''.join(diff).encode('utf8'))
    (OUT/'patch_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('Prepared stage-aware patch', len(changed))


if __name__ == '__main__':
    main()
