"""生成固定SU2的显式变比热补丁，所有替换均限定函数和唯一原语句。"""
from pathlib import Path
import difflib
import hashlib
import json
import zipfile

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'tools/runtime/su2-thermal-build'
OUT=ROOT/'outputs/v154_validation/patch01'


def replace(text,old,new,count=1):
    assert text.count(old)==count,(old[:120],text.count(old),count)
    return text.replace(old,new)


def patches(original):
    changed={}
    name='Common/include/option_structure.hpp'
    text=replace(original[name],'  FLUID_CANTERA = 13,','  IDEAL_AIR_VARIABLE_CP = 14, // RocketAeroHeat固定组分变比热空气\n  FLUID_CANTERA = 13,')
    changed[name]=replace(text,'  MakePair("COOLPROP", COOLPROP)','  MakePair("IDEAL_AIR_VARIABLE_CP", IDEAL_AIR_VARIABLE_CP)\n  MakePair("COOLPROP", COOLPROP)')
    name='Common/src/CConfig.cpp'
    guard='''  // 模块用途。只接受已明确接入的直接SI一般EOS路径。
  if (Kind_FluidModel == IDEAL_AIR_VARIABLE_CP) {
    if (SystemMeasurements != SI || Ref_NonDim != DIMENSIONAL || Time_Domain ||
        Multizone_Problem || DiscreteAdjoint || ContinuousAdjoint ||
        (Kind_Solver != MAIN_SOLVER::NAVIER_STOKES && Kind_Solver != MAIN_SOLVER::RANS) ||
        (Kind_Turb_Model != TURB_MODEL::NONE && Kind_Turb_Model != TURB_MODEL::SST) ||
        Kind_Upwind_Flow != UPWIND::ROE || nMarker_Wall_Functions != 0)
      SU2_MPI::Error("IDEAL_AIR_VARIABLE_CP requires direct SI steady single-zone NS/SST, ROE, resolved wall", CURRENT_FUNCTION);
  }

'''
    changed[name]=replace(original[name],'  /*--- Check if CoolProp is used with non-dimensionalization. ---*/',guard+'  /*--- Check if CoolProp is used with non-dimensionalization. ---*/')
    name='SU2_CFD/src/solvers/CEulerSolver.cpp'
    text=replace(original[name],'#include "../../include/fluid/CCoolProp.hpp"',
        '#include "../../include/fluid/CCoolProp.hpp"\n#include "../../include/fluid/CSU2IdealAir.hpp"\n#include "../../include/fluid/SU2EnergyState.hpp"')
    text=replace(text,'    case COOLPROP:\n\n      auxFluidModel = new CCoolProp(config->GetFluid_Name());',
        '    case IDEAL_AIR_VARIABLE_CP:\n      config->SetGas_Constant(rah_air::gas_constant);\n      auxFluidModel = new CSU2IdealAir();\n      break;\n\n    case COOLPROP:\n\n      auxFluidModel = new CCoolProp(config->GetFluid_Name());')
    text=replace(text,'      case COOLPROP:\n        FluidModel[thread] = new CCoolProp(config->GetFluid_Name());',
        '      case IDEAL_AIR_VARIABLE_CP:\n        FluidModel[thread] = new CSU2IdealAir();\n        break;\n\n      case COOLPROP:\n        FluidModel[thread] = new CCoolProp(config->GetFluid_Name());')
    text=replace(text,'    case COOLPROP:\n      ModelTable << "COOLPROP";',
        '    case IDEAL_AIR_VARIABLE_CP:\n      ModelTable << "IDEAL_AIR_VARIABLE_CP Lemmon2000 OTH 60..2000K";\n      break;\n    case COOLPROP:\n      ModelTable << "COOLPROP";')
    old='''      if (!ideal_gas || low_mach_corr) {
        ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_i, Secondary_i);
        ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_j, Secondary_j);
      }'''
    new='''      bool rah_invalid_reconstruction = false;
      if (!ideal_gas || low_mach_corr) {
        const bool variable_air = config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP;
        const bool with_k = variable_air && config->GetKind_Turb_Model() == TURB_MODEL::SST;
        const su2double k_i = with_k ? solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0) : 0.0;
        const su2double k_j = with_k ? solver_container[TURB_SOL]->GetNodes()->GetSolution(jPoint,0) : 0.0;
        try {
          ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_i, Secondary_i, k_i);
          ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_j, Secondary_j, k_j);
        } catch (const std::domain_error&) {
          if (!variable_air) throw;
          // 拒绝超温域重构，复用原生单元值降阶路径，避免半更新状态。
          rah_invalid_reconstruction = true;
          std::copy_n(V_i,nPrimVar,Primitive_i); std::copy_n(V_j,nPrimVar,Primitive_j);
          std::copy_n(S_i,nSecondaryVar,Secondary_i); std::copy_n(S_j,nSecondaryVar,Secondary_j);
        }
      }'''
    text=replace(text,old,new)
    text=replace(text,'bool bad_recon = neg_sound_speed || neg_pres_or_rho_i || neg_pres_or_rho_j;',
        'bool bad_recon = rah_invalid_reconstruction || neg_sound_speed || neg_pres_or_rho_i || neg_pres_or_rho_j;')
    text=replace(text,'su2double *primitive, su2double *secondary) {',
        'su2double *primitive, su2double *secondary, su2double turbulent_k) {')
    text=replace(text,'primitive[prim_idx.Enthalpy()] = fluidModel->GetStaticEnergy() + pressure / density + 0.5*velocity2;',
        'primitive[prim_idx.Enthalpy()] = fluidModel->GetStaticEnergy() + pressure / density + 0.5*velocity2 + turbulent_k;')

    # 仅修改新模型Riemann路径，热状态扣k，特征矩阵使用含k守恒χ。
    start=text.index('void CEulerSolver::BC_Riemann(')
    end=text.index('\nvoid CEulerSolver::',start+1)
    boundary=text[start:end]
    boundary=replace(boundary,'      const su2double StaticEnergy_i = Energy_i - 0.5*Velocity2_i;',
        '''      const bool variable_air = config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP;
      const su2double k_i = variable_air && tkeNeeded ? solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0) : 0.0;
      const su2double StaticEnergy_i = rah_air::thermal_energy(Energy_i,Velocity2_i,k_i);''')
    boundary=replace(boundary,'GetFluidModel()->GetdPdrho_e() - Kappa_i * StaticEnergy_i;',
        'GetFluidModel()->GetdPdrho_e() - Kappa_i * (StaticEnergy_i+k_i);')
    anchor='''        Energy_e = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_e;
        break;'''
    boundary=replace(boundary,anchor,'''        Energy_e = GetFluidModel()->GetStaticEnergy() + 0.5*Velocity2_e + k_i;
        break;''')
    boundary=replace(boundary,'      const su2double StaticEnergy_b = Energy_b - 0.5*Velocity2_b;',
        '''      const bool outgoing = GeometryToolbox::DotProduct(nDim,Velocity_b,UnitNormal) >= 0;
      const bool pressure_boundary = config->GetKind_Data_Riemann(Marker_Tag) == STATIC_PRESSURE;
      const su2double k_b = variable_air && tkeNeeded ? ((outgoing || pressure_boundary) ? k_i : GetTke_Inf()) : 0.0;
      const su2double StaticEnergy_b = rah_air::thermal_energy(Energy_b,Velocity2_b,k_b);''')
    boundary=replace(boundary,'GetFluidModel()->GetdPdrho_e() - Kappa_b * StaticEnergy_b;',
        'GetFluidModel()->GetdPdrho_e() - Kappa_b * (StaticEnergy_b+k_b);')
    text=text[:start]+boundary+text[end:]
    changed[name]=text
    name='SU2_CFD/include/solvers/CEulerSolver.hpp'
    changed[name]=replace(original[name],'                                             su2double *secondary);',
        '                                             su2double *secondary, su2double turbulent_k = 0.0);')
    name='SU2_CFD/src/solvers/CNSSolver.cpp'
    text=original[name]
    start=text.index('void CNSSolver::BC_Isothermal_Wall_Generic(')
    end=text.index('\nvoid CNSSolver::',start+1)
    wall=text[start:end]
    wall=replace(wall,'      su2double dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );',
        '''      su2double dTdrho = 1.0/Density * ( -Twall + (Gamma-1.0)/Gas_Constant*(Vel2/2.0) );
      su2double dTdrhoE = (Gamma-1.0)/(Gas_Constant*Density);
      if (config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP) {
        // 原生弱壁面近似对角预条件，采用同阶段实际状态和冻结k的热力学链式导数。
        const su2double* secondary = nodes->GetSecondary(iPoint);
        dTdrho = secondary[2] + secondary[3]*(Vel2-nodes->GetEnergy(iPoint))/Density;
        dTdrhoE = secondary[3]/Density;
      }''')
    wall=replace(wall,'thermal_conductivity/dist_ij * (Gamma-1.0)/(Gas_Constant*Density) * Area;',
        'thermal_conductivity/dist_ij * dTdrhoE * Area;')
    changed[name]=text[:start]+wall+text[end:]
    return changed


def main():
    assert not OUT.exists(); OUT.mkdir(parents=True)
    archive=BASE/'su2.zip'
    with zipfile.ZipFile(str(archive)) as bundle:
        originals={'/'.join(name.split('/')[1:]):bundle.read(name).decode('utf8') for name in bundle.namelist()
            if not name.endswith('/') and name.endswith(('.cpp','.hpp'))}
    changed=patches(originals)
    diff=[]
    for name,text in changed.items():
        path=OUT/name; path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(text.encode('utf8'))
        diff.extend(difflib.unified_diff(originals[name].splitlines(True),text.splitlines(True),fromfile='a/'+name,tofile='b/'+name))
    (OUT/'ideal_air.diff').write_bytes(''.join(diff).encode('utf8'))
    manifest={name:dict(originalSHA256=hashlib.sha256(originals[name].encode()).hexdigest(),
        patchedSHA256=hashlib.sha256(text.encode()).hexdigest()) for name,text in changed.items()}
    (OUT/'patch_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Prepared exact patches',len(changed))


if __name__=='__main__':
    main()
