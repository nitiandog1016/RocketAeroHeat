"""模块用途：有限只读确认v171空间Tw包的原生IDEAL_GAS路径、准确源码/对象与输出合同。"""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
BUILD=ROOT/'outputs/v171_validation/build01'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
OBJECTS=ROOT/'tools/runtime/su2-thermal-build/build-v117'
PACKAGE=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
PLAN=ROOT/'docs/validation/v174/PLAN.md'
TREE_SHA='c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8'
EXE_SHA='3309e77ed0fdfdb465e6d774f2157b7ef8d2cb2252a3d4abbfab7ee49fd7a217'
SOURCES={
    'Common/src/CConfig.cpp':[(1326,1329),(3793,3805),(5772,5781)],
    'Common/src/geometry/CGeometry.cpp':[(2752,2817)],
    'Common/src/geometry/RAHAxialWallTemperature.hpp':[(16,77)],
    'SU2_CFD/src/python_wrapper_structure.cpp':[(1,20)],
    'SU2_CFD/src/fluid/CIdealGas.cpp':[(31,101)],
    'SU2_CFD/src/fluid/CFluidModel.cpp':[(79,108),(131,134)],
    'SU2_CFD/include/fluid/CConstantPrandtl.hpp':[(42,60)],
    'SU2_CFD/include/fluid/SU2EnergyState.hpp':[(1,10)],
    'SU2_CFD/src/variables/CEulerVariable.cpp':[(31,64)],
    'SU2_CFD/include/variables/CEulerVariable.hpp':[(48,65),(124,148),(201,214)],
    'SU2_CFD/src/variables/CNSVariable.cpp':[(140,231)],
    'SU2_CFD/src/solvers/CEulerSolver.cpp':[(873,891),(1140,1179),(1838,1850),(1964,1982),
        (5290,5296),(5331,5348),(5399,5427),(5469,5488),(5527,5558)],
    'SU2_CFD/src/solvers/CNSSolver.cpp':[(690,750)],
    'SU2_CFD/src/drivers/CDriver.cpp':[(1667,1682),(1804,1828)],
    'SU2_CFD/src/numerics/flow/flow_diffusion.cpp':[(393,446),(486,510)],
    'SU2_CFD/include/numerics_simd/CNumericsSIMD.cpp':[(38,55),(97,119)],
    'SU2_CFD/include/numerics_simd/flow/variables.hpp':[(40,95)],
    'SU2_CFD/include/numerics_simd/flow/convection/common.hpp':[(74,130)],
    'SU2_CFD/include/numerics_simd/flow/convection/upwind.hpp':[(51,76),(112,148)],
    'SU2_CFD/include/numerics_simd/flow/diffusion/viscous_fluxes.hpp':[(252,300)],
    'SU2_CFD/src/output/CFlowCompOutput.cpp':[(220,276),(309,327)],
    'SU2_CFD/src/output/filewriter/CSU2FileWriter.cpp':[(44,68)],
}
OBJECT_ENDINGS=('CConfig.cpp.obj','geometry_CGeometry.cpp.obj','python_wrapper_structure.cpp.obj',
    'fluid_CIdealGas.cpp.obj','fluid_CFluidModel.cpp.obj','variables_CEulerVariable.cpp.obj',
    'variables_CNSVariable.cpp.obj','solvers_CEulerSolver.cpp.obj','solvers_CNSSolver.cpp.obj',
    'drivers_CDriver.cpp.obj','flow_flow_diffusion.cpp.obj','numerics_simd_CNumericsSIMD.cpp.obj',
    'output_CFlowCompOutput.cpp.obj','filewriter_CSU2FileWriter.cpp.obj')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')


def main(output):
    """模块一：准确包/源/对象绑定；模块二：具名分派断言；模块三：独立接口要求。"""
    started=time.monotonic();out=Path(output).resolve();out.relative_to(ROOT/'outputs')
    out.mkdir(parents=True,exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    locks={};report=dict(passed=False,CFDLaunches=0,builds=0,largeFieldsRead=0,sharedSourceMutations=0)
    def bind(path,expected=None):
        digest=sha(path)
        if expected is not None:assert digest==expected,str(path)
        locks[str(path.resolve())]=digest
        return digest
    try:
        for path in (Path(__file__),PLAN,PACKAGE/'runtime.json',BUILD/'runtime_installed.json',BUILD/'source_tree.json',
                     BUILD/'objects_after.json',BUILD/'archive_members.json',BUILD/'commands.json',BUILD/'link_process.json'):
            bind(path)
        plan=PLAN.read_text('utf8')
        assert 'v174' in plan and 'IDEAL_GAS' in plan and '40秒' in plan
        (out/'plan_at_review.md').write_bytes(PLAN.read_bytes())
        runtime=json.loads((PACKAGE/'runtime.json').read_text('utf8'))
        assert runtime==json.loads((BUILD/'runtime_installed.json').read_text('utf8'))
        assert runtime['flavor']=='rah-ideal-air-spatial-wall-v1' and runtime['asciiSignificantDigits']==17
        assert runtime['sourceTreeSHA256']==TREE_SHA==bind(BUILD/'source_tree.json')
        exe=PACKAGE/'bin/SU2_CFD.exe';assert bind(exe)==EXE_SHA
        assert EXE_SHA==next(row['sha256'] for row in runtime['files'] if row['path']=='bin/SU2_CFD.exe')
        assert json.loads((BUILD/'link_process.json').read_text('utf8'))['exitCode']==0
        tree=json.loads((BUILD/'source_tree.json').read_text('utf8'))
        texts={};source_records=[];excerpts=[]
        for relative,intervals in SOURCES.items():
            shared=SOURCE/relative;expected=tree[str(shared)]
            candidates=[BUILD/'executed_sources'/relative,
                ROOT/'outputs/v163_validation/build03/executed_sources'/relative,shared]
            selected=next((p for p in candidates if p.is_file() and sha(p)==expected),None)
            assert selected is not None,'缺准确源: '+relative
            bind(selected,expected);texts[relative]=selected.read_text('utf8')
            copy=out/'executed_sources'/relative;copy.parent.mkdir(parents=True,exist_ok=True);copy.write_bytes(selected.read_bytes())
            source_records.append(dict(relativePath=relative,selectedPath=str(selected),sha256=expected,
                exactCopy=str(copy),sharedSourceMatches=shared.is_file() and sha(shared)==expected))
            lines=texts[relative].splitlines()
            for first,last in intervals:
                excerpts.append('\n'+relative+':'+str(first)+'-'+str(last)+' SHA256='+expected+'\n')
                excerpts.extend(str(i+1)+': '+lines[i]+'\n' for i in range(first-1,min(last,len(lines))))
        (out/'exact_source_excerpts.txt').write_text(''.join(excerpts),encoding='utf8')
        objects=json.loads((BUILD/'objects_after.json').read_text('utf8'));assert len(objects)==298
        object_records=[]
        for ending in OBJECT_ENDINGS:
            found=[name for name in objects if name.endswith(ending)];assert len(found)==1,(ending,found)
            name=found[0];expected=objects[name]
            choices=[BUILD/'produced_objects'/name,OBJECTS/name]
            selected=next((p for p in choices if p.is_file() and sha(p)==expected),None)
            assert selected is not None,'缺准确对象: '+name
            bind(selected,expected);object_records.append(dict(object=name,selectedPath=str(selected),sha256=expected))

        cfg=texts['Common/src/CConfig.cpp'];gas=texts['SU2_CFD/src/fluid/CIdealGas.cpp']
        euler=texts['SU2_CFD/src/solvers/CEulerSolver.cpp'];node=texts['SU2_CFD/src/variables/CNSVariable.cpp']
        wall=texts['SU2_CFD/src/solvers/CNSSolver.cpp'];geom=texts['Common/src/geometry/CGeometry.cpp']
        output_text=texts['SU2_CFD/src/output/CFlowCompOutput.cpp']
        checks={}
        def check(name,value):
            checks[name]=bool(value);assert checks[name],name
        check('nativeGasDispatch', 'auxFluidModel = new CIdealGas(Gamma, config->GetGas_Constant());' in euler and
              'FluidModel[thread] = new CIdealGas(Gamma, Gas_ConstantND);' in euler)
        check('gammaAndRExplicitOptions','addDoubleOption("GAMMA_VALUE", Gamma, 1.4)' in cfg and
              'addDoubleOption("GAS_CONSTANT", Gas_Constant, 287.058)' in cfg)
        for token in ('Cp = Gamma / Gamma_Minus_One * Gas_Constant;','Cv = Cp - R;',
                      'Temperature = dTde_rho * StaticEnergy;',
                      'su2double e = T * Gas_Constant / Gamma_Minus_One;',
                      'Pressure = dPde_rho * StaticEnergy;', 'SoundSpeed2 = Gamma * Gas_Constant * Temperature;'):
            check('gas:'+token,token in gas)
        check('nativeCNSNoProjectBacktracking','rah_step' not in node and 'catch' not in node)
        check('nativeOldStateRecoveryRetained','Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);' in node)
        check('nativeGasSIMDFlowPath','if (ideal_gas && !low_mach_corr)' in euler and
              'EdgeFluxResidual(geometry, solver_container, config);\n      return;' in euler)
        factory=texts['SU2_CFD/include/numerics_simd/CNumericsSIMD.cpp']
        check('nativeGasRoeViscousFactory','new CRoeScheme<ViscousDecorator>' in factory and
              'createUpwindIdealNumerics<CCompressibleViscousFlux<nDim> >' in factory)
        recon=texts['SU2_CFD/include/numerics_simd/flow/convection/common.hpp']
        check('nativeGasReconstructionAndFallback','V.i.density() = V.i.pressure() / (gasConst * V.i.temperature());' in recon and
              'solution.UpdateNonPhysicalEdgeCounter' in recon and 'V.i.all(iVar) = nonPhysical * V1st.i.all(iVar)' in recon)
        secondary=texts['SU2_CFD/src/variables/CEulerVariable.cpp']
        check('nativeGasSecondaryZero','if (ideal_gas) return 0;' in secondary and 'ndim + 9, EulerNPrimVarGrad' in secondary)
        idx=texts['SU2_CFD/include/variables/CEulerVariable.hpp']
        simdidx=texts['SU2_CFD/include/numerics_simd/flow/variables.hpp']
        check('conductivityCpIndices','ThermalConductivity() const { return nDim+7; }' in idx and
              'CpTotal() const { return nDim+8; }' in idx and 'thermalCond() { return all(nDim+7); }' in simdidx and
              'cp() { return all(nDim+8); }' in simdidx)
        check('nativeConstantPrandtl','kt_ = mu_lam * cp / pr_lam_;' in texts['SU2_CFD/include/fluid/CConstantPrandtl.hpp'])
        check('noRANSConductivityForNS','Kind_ConductivityModel_Turb = CONDUCTIVITYMODEL_TURB::NONE;' in cfg)
        check('weakWallFormulaPreserved','su2double dTdn = -(There - Twall)/dist_ij;' in wall and
              'su2double Res_Visc = thermal_conductivity * dTdn * Area;' in wall)
        check('variableCpWallDerivativeGuarded','if (config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP) {' in wall and
              'dTdrhoE = secondary[3]/Density;' in wall)
        check('riemannOriginalMarkerBranches','case STATIC_SUPERSONIC_INFLOW_PT:' in euler and
              'Velocity_e[iDim] = Mach[iDim]*GetFluidModel()->GetSoundSpeed();' in euler and
              'case STATIC_PRESSURE:' in euler and 'GetFluidModel()->SetTDState_Prho(Pressure_e, Density_e);' in euler)
        check('riemannEnergyHelperOnlyKineticSubtraction',
              'return total_energy-.5*velocity_squared-turbulent_k;' in texts['SU2_CFD/include/fluid/SU2EnergyState.hpp'])
        begin=geom.index('void CGeometry::SetCustomBoundary(CConfig* config) {')
        end=geom.index('void CGeometry::UpdateCustomBoundaryConditions',begin)
        thermal=geom[begin:end]
        check('spatialTwIndependentOfFluidModel','RAHAxialWallTemperature rah;' in thermal and 'GetKind_FluidModel' not in thermal)
        check('spatialTwStrictSIAndSingleWall','config->GetnMGLevels() != 0' in thermal and
              'Marker_Tag != "wall"' in thermal and 'std::setprecision(17)' in thermal)
        check('restart17DigitsForAllModels','data.precision(std::numeric_limits<double>::max_digits10 - 1);' in
              texts['SU2_CFD/src/output/filewriter/CSU2FileWriter.cpp'])
        check('flowPrimitiveKVariableCpOnly','if (config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP)\n    AddVolumeOutput("FLOW_PRIMITIVE_K"' in output_text)
        check('noNewStepOrQueryPatchInV171',not any('PropertyQueryContext' in value or 'ConservativeStepTrace' in value for value in texts.values()))
        primitive=dict(T=0,u=1,v=2,w=3,p=4,rho=5,H=6,a=7,mu=8,mu_t=9,kappa=10,cp=11)
        columns=['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy',
            'Pressure','Temperature','Mach','Pressure_Coefficient','Velocity_x','Velocity_y','Velocity_z',
            'Laminar_Viscosity','Skin_Friction_Coefficient_x','Skin_Friction_Coefficient_y','Skin_Friction_Coefficient_z',
            'Heat_Flux','Y_Plus','Residual_Density','Residual_Momentum_x','Residual_Momentum_y','Residual_Momentum_z','Residual_Energy']
        report.update(passed=True,reviewType='exact-source-native-ideal-gas-path',runtimeFlavor=runtime['flavor'],
            runtimeExecutableSHA256=EXE_SHA,sourceTreeSHA256=TREE_SHA,runtimeMetadataChanged=False,
            packageIsUnmodifiedOfficialDistribution=False,originalPackageModelMetadata=runtime['thermodynamicModel'],
            proposedCaseModel='SU2-native-IDEAL_GAS-gamma1.4-R8.31451/0.02896546',
            sourceBindings=source_records,objectBindings=object_records,allRecordedObjectCount=298,checks=checks,
            primitive3DZeroBased=primitive,primitiveCount3D=12,secondaryCount=0,reconstructedGradientCount3D=5,
            expectedRestartColumnsForExplicitGroups=columns,Flow_Primitive_K_expected=False,
            internalEnergyReference='e=cv*T, zero at T=0 in calorically-perfect model; not Lemmon OTH',
            oldOTHConservativeRestartAllowed=False,proposedBoundary='retain v171 outer RIEMANN static PT and outlet static pressure',
            riemannFlowInputIsMachVector=True,farfieldChangeRecommended=False,
            sourceIndicatesBinaryReusableForSmallGasTrial=True,newNativeGasRunVerified=False,
            nativeIdealGasHas60To2000KGuard=False,proposedDevelopmentEnvelopeK=[150,360],envelopeIsExternalAuditNotClipping=True,
            originalResidualGate=-8,strictHeatGradientScaledGate=2e-10,
            finitePositiveStateMustBeIndependentlyAudited=True,convergenceOrPhysicalQualificationProven=False,
            requiredIndependentContract=['exact v171 binary/source/object identity','IDEAL_GAS and explicit gamma/R',
                'new same-geometry conservative state in cv*T reference','single nonuniform Tw table and applied receipt',
                'same RIEMANN markers with Mach vector [-1.97,0,0]','model-specific 27-column saved state',
                'separate runtime metadata from actual case physics','phase/capacity/inventory dispatch before launch',
                'normal/abnormal saved-frame separation','independent EOS/energy/heat/method labels'])
    finally:
        stable={path:Path(path).is_file() and sha(Path(path))==digest for path,digest in locks.items()}
        if not all(stable.values()):report['passed']=False
        report.update(inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started)
        save(out/'report.json',report)
    assert report['passed'] and all(report['inputsStable'].values())
    print(json.dumps(dict(passed=True,elapsedSeconds=report['elapsedSeconds'],sources=len(source_records),
        objects=len(object_records),checks=len(checks),CFD=0,builds=0,gasRunVerified=False)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    args=parser.parse_args();main(args.output)
