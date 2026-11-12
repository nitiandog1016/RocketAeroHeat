"""模块职责：暂存独立Euler身份的最小适配补丁；旧源码/对象/引擎保持。
只扩注册和既有相态气体分派，不修改通量公式、推进、线搜索或滑移内核。
源码静态审查不能替代编译及原生Euler兼容控制。
"""
from pathlib import Path
import difflib,json
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from preparePhaseEulerV228 import BASE,bound_sources,IDENTITY

GUARD='''  // RocketAeroHeat: 独立无黏身份；旧RAH_PHASE_AIR的NS/SA合同保持。
  if (Kind_FluidModel == RAH_PHASE_AIR_EULER) {
    if (SystemMeasurements != SI || Ref_NonDim != DIMENSIONAL || Time_Domain ||
        Multizone_Problem || DiscreteAdjoint || ContinuousAdjoint ||
        Kind_Solver != MAIN_SOLVER::EULER || Kind_Turb_Model != TURB_MODEL::NONE ||
        Kind_Upwind_Flow != UPWIND::ROE || nMarker_WallFunctions != 0 ||
        nMarker_Euler == 0 || nMarker_Riemann == 0 ||
        nMarker_Isothermal || nMarker_HeatFlux || nMarker_HeatTransfer ||
        nMarker_FarField || nMarker_Supersonic_Inlet || nMarker_Supersonic_Outlet ||
        nMarker_Inlet || nMarker_Outlet || Low_Mach_Corr || nMarker_Giles || GetDynamic_Grid())
      SU2_MPI::Error("RAH_PHASE_AIR_EULER requires direct SI steady EULER/NONE, ROE, slip wall and ambient PT characteristics", CURRENT_FUNCTION);
    for (unsigned short i=0; i<nMarker_Riemann; ++i)
      if (Kind_Data_Riemann[i] != STATIC_SUPERSONIC_INFLOW_PT)
        SU2_MPI::Error("RAH_PHASE_AIR_EULER supports ambient PT characteristic boundary only", CURRENT_FUNCTION);
  }

'''

def transform(rel,s):
    if rel=='Common/include/option_structure.hpp':
        a='  RAH_PHASE_AIR = 15, // 相态标准Air，内部能量+300000J/kg'
        assert s.count(a)==1
        s=s.replace(a,'  RAH_PHASE_AIR_EULER = 16, // 独立Euler分派，同相态/能量坐标\n'+a)
        a='  MakePair("RAH_PHASE_AIR", RAH_PHASE_AIR)';assert s.count(a)==1
        return s.replace(a,'  MakePair("RAH_PHASE_AIR_EULER", RAH_PHASE_AIR_EULER)\n'+a)
    if rel=='Common/src/CConfig.cpp':
        marker='  /*--- Check if CoolProp is used with non-dimensionalization. ---*/';assert s.count(marker)==1
        return s.replace(marker,GUARD+marker)
    if rel=='SU2_CFD/src/solvers/CEulerSolver.cpp':
        # 只让新身份沿现有气体/重构/一般EOS特征分支，不改它们的数学正文。
        assert s.count('case RAH_PHASE_AIR:')==3
        s=s.replace('case RAH_PHASE_AIR:','case RAH_PHASE_AIR_EULER:\n    case RAH_PHASE_AIR:')
        a='config->GetKind_FluidModel() == RAH_PHASE_AIR';assert s.count(a)==3
        s=s.replace(a,'('+a+' || config->GetKind_FluidModel() == RAH_PHASE_AIR_EULER)')
        return s
    if rel=='SU2_CFD/src/numerics/flow/convection/roe.cpp':
        a='config->GetKind_FluidModel() == RAH_PHASE_AIR';assert s.count(a)==1
        return s.replace(a,'('+a+' || config->GetKind_FluidModel() == RAH_PHASE_AIR_EULER)')
    raise ValueError(rel)

def main():
    out=BASE/'stage01';out.mkdir(exist_ok=False);sources,text=bound_sources();patch=[];manifest={}
    changed=['Common/include/option_structure.hpp','Common/src/CConfig.cpp',
             'SU2_CFD/src/solvers/CEulerSolver.cpp','SU2_CFD/src/numerics/flow/convection/roe.cpp']
    for rel in changed:
        before=text[rel];after=transform(rel,before)
        path=out/'isolated_source'/rel;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(after,encoding='utf8')
        patch.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='a/'+rel,tofile='b/'+rel))
        manifest[rel]=dict(original=sources[rel],isolatedPath=str(path),patchedSHA256=sha(path))
    # 严格可逆变换保护旧NS配置；无语义授予，无当前引擎替换。
    assert (out/'isolated_source/Common/src/CConfig.cpp').read_text(encoding='utf8').replace(GUARD,'')==text['Common/src/CConfig.cpp']
    assert all(sha(v['path'])==v['sha256'] for v in sources.values())
    (out/'phase_euler_identity.patch').write_text(''.join(patch),encoding='utf8')
    save(out/'manifest.json',manifest)
    save(out/'report.json',dict(passed=True,identity=IDENTITY,staticSourceOnly=True,changedSourceFiles=4,
      oldNSGuardTextUnchanged=True,oldSourceUnchanged=True,oldEngineUnchanged=True,
      nativeGasClassUnchanged=True,slipBoundaryKernelUnchanged=True,primitiveUpdateKernelUnchanged=True,
      numericalFormulaChanges=0,builds=0,nativeQueries=0,CFDRuns=0,compiled=False,
      actualEulerCompatibilityPassed=False,runAuthorized=False,newQualifiedAreaM2=0,
      pending='one isolated build and whole Euler slip/primitive/characteristic/save-restart control under a new budget'))
    print(json.dumps(read(out/'report.json')))
if __name__=='__main__':main()
