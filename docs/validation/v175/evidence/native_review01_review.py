"""模块用途：有限核对准确热壁与外边界原生实现；不运行CFD、不构建、不读取大场。"""
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/v175_validation/native_review01'
PLAN=ROOT/'docs/validation/v175/PLAN.md'
BUILD=ROOT/'outputs/v174_validation/build01'
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
RUNTIME=ROOT/'tools/runtime/su2/8.6.0-rah-native-ideal-gas-spatial-wall-v1/runtime.json'
TREE_SHA='7f5df95b0c5213266f21d8a86abf23f0d4ebc8f1b21ff12a169d04f72087720b'
RUNTIME_SHA='faaed7f318876335f5309fb79bbbc0facb2d50a9972bfaaa94b97220f9457ebb'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text('utf8'))
def save(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf8')


def main():
    """模块一：仅准确源与小配置绑定；模块二：保留可复核行号和完整原文字节。"""
    started=time.monotonic();OUT.mkdir(parents=True,exist_ok=False)
    locked={str(p):sha(p) for p in (Path(__file__),PLAN,RUNTIME,BUILD/'source_tree.json',BUILD/'objects_after.json')}
    assert sha(RUNTIME)==RUNTIME_SHA and sha(BUILD/'source_tree.json')==TREE_SHA
    runtime=read(RUNTIME);tree=read(BUILD/'source_tree.json')
    assert runtime['sourceTreeSHA256']==TREE_SHA and runtime['supportedModel']=='IDEAL_GAS'
    bindings={};texts={}
    names=(
        'SU2_CFD/src/solvers/CNSSolver.cpp',
        'SU2_CFD/src/solvers/CEulerSolver.cpp',
        'Common/src/CConfig.cpp',
        'SU2_CFD/src/integration/CIntegration.cpp',
        'config_template.cfg',
        'TestCases/navierstokes/flatplate/flatplate_supersonic.cfg',
        'TestCases/navierstokes/flatplate/lam_flatplate.cfg',
        'TestCases/serial_regression.py')
    for name in names:
        accurate=BUILD/'executed_sources'/name
        if not accurate.is_file():accurate=SOURCE/name
        digest=sha(accurate);assert digest==tree[str(SOURCE/name)]
        locked[str(accurate)]=digest
        destination=OUT/'executed_sources'/name;destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_bytes(accurate.read_bytes())
        texts[name]=accurate.read_text('utf8')
        bindings[name]=dict(accurateSource=str(accurate),sourceTreeKey=str(SOURCE/name),sha256=digest,
            savedSource=str(destination),matchesBuiltTree=True)
    ns=texts[names[0]];euler=texts[names[1]];config=texts[names[2]];integration=texts[names[3]]
    wall=ns[ns.index('void CNSSolver::BC_Isothermal_Wall_Generic('):ns.index('void CNSSolver::BC_Isothermal_Wall(')]
    far=euler[euler.index('void CEulerSolver::BC_Far_Field('):euler.index('void CEulerSolver::BC_Riemann(')]
    riemann=euler[euler.index('void CEulerSolver::BC_Riemann('):euler.index('\nvoid CEulerSolver::',euler.index('void CEulerSolver::BC_Riemann(')+1)]
    checks={
        'weakWallUsesPrescribedOrCustomTw':'GetCustomBoundaryTemperature' in wall and 'GetIsothermal_Temperature' in wall,
        'weakWallNeighborGradient':'dTdn = -(There - Twall)/dist_ij' in wall,
        'weakWallEnergyResidual':'LinSysRes(iPoint, nDim+1) += Res_Conv - Res_Visc' in wall,
        'strongRowsOnlyMomentum':'for (auto iVar = 1u; iVar <= nDim; iVar++)' in wall,
        'chtUsesSameWeakWallGeneric':'config, val_marker, true)' in ns,
        'riemannPTBuildsMachScaledExternalState':'Velocity_e[iDim] = Mach[iDim]*GetFluidModel()->GetSoundSpeed()' in riemann,
        'riemannSelectsIncomingInteriorCharacteristics':'if (Lambda_i[jVar] < 0)' in riemann,
        'riemannUsesInteriorAcousticEigenvalues':'Lambda_i[nVar-2] = ProjVelocity_i + SoundSpeed_i' in riemann,
        'riemannUsesProjectedPhysicalFlux':'conv_numerics->GetInviscidProjFlux' in riemann,
        'farfieldUsesFreestreamAcousticSign':'if (Qn_Infty > -SoundSpeed_Infty)' in far and 'if (Qn_Infty > SoundSpeed_Infty)' in far,
        'farfieldUsesIncomingOutgoingEntropy':'Entropy = Entropy_Bound' in far and 'Entropy = Entropy_Infty' in far,
        'farfieldUsesNativeUpwindFlux':'conv_numerics->SetPrimitive(V_domain, V_infty)' in far and 'conv_numerics->ComputeResidual(config)' in far,
        'idealGasRecognized':'(Kind_FluidModel == IDEAL_GAS)' in config,
        'nonIdealFarfieldRejected':'if (nMarker_FarField != 0)' in config,
        'farfieldAndRiemannHaveSeparateNativeDispatch':'case FAR_FIELD:' in integration and 'case RIEMANN_BOUNDARY:' in integration,
        'officialSupersonicThermalFarfield':'MARKER_ISOTHERMAL= ( WALL, 300.0 )' in texts[names[5]] and 'MARKER_FAR= ( INLET )' in texts[names[5]],
        'officialSupersonicNeedsRestart':'RESTART_SOL= YES' in texts[names[5]],
        'standardFlatplateIsAdiabatic':'MARKER_HEATFLUX= ( wall, 0.0 )' in texts[names[6]],
    }
    assert all(checks.values())
    excerpts={}
    for name,begin,end in ((names[0],594,804),(names[1],5028,5282),(names[1],5285,5570),
            (names[2],3690,3703),(names[2],4620,4637),(names[3],120,144),(names[4],1175,1195)):
        excerpts[name+':'+str(begin)]=dict(firstLine=begin,lastLine=end,
            text='\n'.join(texts[name].splitlines()[begin-1:end]))
    save(OUT/'source_bindings.json',bindings);save(OUT/'source_excerpts.json',excerpts)
    (OUT/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    (OUT/'plan_at_review.md').write_bytes(PLAN.read_bytes())
    stable={name:sha(Path(name))==digest for name,digest in locked.items()}
    assert all(stable.values())
    report=dict(passed=True,scope='exact-source-and-configuration-only',checks=checks,
        sourceFiles=len(names),sourceTreeSHA256=TREE_SHA,runtimeSHA256=RUNTIME_SHA,
        sourceCommit=runtime['sourceCommit'],inputSHA256=locked,inputsStable=stable,
        existingRiemannIsAllStateForced=False,farfieldAndRiemannEquivalent=False,
        nativeIdealGasFarfieldSupported=True,variableCpFarfieldSupported=False,
        weakWallReplacedByStrongTw=False,chtEliminatesWeakWall=False,
        selectedMatureAlternative='native-IDEAL_GAS-external-characteristic-farfield-on-outer-and-outlet',
        alternativeChangesBoundaryStateSelectionAndFlux=True,
        completeOfficialSupersonicCasePresentLocally=False,
        officialConfigAloneIsIndependentHeatValidation=False,
        nacaHeatReferenceOrConvergenceClosed=False,
        candidateRequiresSeparatePredeclaredPhysicalContract=True,
        nacaTrialAuthorizedByThisReview=False,oldNacaFailureOr360KLimitChanged=False,
        physicalQualified=False,newQualifiedArea_m2=0,newCFDTasks=0,newBuilds=0,
        largeFieldFilesRead=0,sharedSourceWrites=0,elapsedSeconds=time.monotonic()-started)
    save(OUT/'report.json',report)
    print(json.dumps(dict(passed=True,checks=len(checks),sourceFiles=len(names),newCFDTasks=0,newBuilds=0,
        elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
