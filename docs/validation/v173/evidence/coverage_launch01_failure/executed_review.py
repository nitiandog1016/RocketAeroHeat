"""模块用途：只读核对准确空间壁温运行包的物性查询覆盖，不编译、不运行CFD。"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import time

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'tools/runtime/su2-thermal-build/source'
BUILD = ROOT / 'outputs/v171_validation/build01'
PACKAGE = ROOT / 'tools/runtime/su2/8.6.0-rah-ideal-air-spatial-wall-v1'
TREE_SHA = 'c600dc43dc777e40a344f253b1adf89e7c51aa95caa98c6a58f57db1c65f0ee8'
EXE_SHA = '3309e77ed0fdfdb465e6d774f2157b7ef8d2cb2252a3d4abbfab7ee49fd7a217'
EULER = 'SU2_CFD/src/solvers/CEulerSolver.cpp'
CNS = 'SU2_CFD/src/variables/CNSVariable.cpp'
WALL = 'SU2_CFD/src/solvers/CNSSolver.cpp'
BASE = 'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl'
ROE = 'SU2_CFD/src/numerics/flow/convection/roe.cpp'
CONFIG = 'Common/src/CConfig.cpp'
SOURCES = {
    CONFIG: [(1451,1458),(3793,3805)],
    EULER: [(310,326),(835,841),(925,934),(960,995),(1170,1178),
            (1964,2013),(2065,2082),(2567,2602),(5285,5562)],
    CNS: [(140,187)],
    WALL: [(630,793)],
    BASE: [(1150,1315)],
    ROE: [(788,824),(985,1021)],
    'SU2_CFD/src/drivers/CDriver.cpp': [(1666,1683)],
    'SU2_CFD/include/fluid/CSU2IdealAir.hpp': [(14,41)],
    'SU2_CFD/include/fluid/IdealAir.hpp': [(13,17),(140,189)],
    'SU2_CFD/src/output/CFlowCompOutput.cpp': [],
}
QUERY = re.compile(r'(?:SetTDState_(?:rhoe|PT|Prho|rhoT|hs|Ps)|SetEnergy_Prho|ComputeDerivativeNRBC_Prho)\s*\(')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')


def section(text, start, end):
    """模块一：在具名函数之间提取边界，防止把相邻函数的查询算入当前路径。"""
    begin = text.index(start)
    return text[begin:text.index(end, begin+len(start))]


def main(output):
    """模块二：绑定二进制、准确源码和实际配置，并逐项记录保护及拒绝范围。"""
    started = time.monotonic()
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    runtime_path = PACKAGE/'runtime.json'
    tree_path = BUILD/'source_tree.json'
    binary = PACKAGE/'bin/SU2_CFD.exe'
    config_path = ROOT/'outputs/v172_validation/real01/case.cfg'
    runtime = json.loads(runtime_path.read_text('utf8'))
    tree = json.loads(tree_path.read_text('utf8'))
    assert sha(tree_path) == runtime['sourceTreeSHA256'] == TREE_SHA
    assert sha(binary) == EXE_SHA == next(x['sha256'] for x in runtime['files'] if x['path']=='bin/SU2_CFD.exe')
    assert runtime['flavor'] == 'rah-ideal-air-spatial-wall-v1'
    locks = {str(p):sha(p) for p in (runtime_path,tree_path,binary,config_path,Path(__file__).resolve())}
    sources, texts, excerpts = [], {}, []
    for relative, intervals in SOURCES.items():
        shared = SOURCE/relative
        expected = tree[str(shared)]
        candidates = [BUILD/'executed_sources'/relative,
                      ROOT/'outputs/v163_validation/build03/executed_sources'/relative,shared]
        selected = next((p for p in candidates if p.is_file() and sha(p)==expected),None)
        assert selected is not None, '没有准确源：'+relative
        locks[str(selected)] = expected
        texts[relative] = selected.read_text('utf8')
        sources.append(dict(relativePath=relative,selectedSource=str(selected),sha256=expected,
                            sharedSourceMatches=shared.is_file() and sha(shared)==expected))
        lines = texts[relative].splitlines()
        for first,last in intervals:
            excerpts.append('\n'+relative+':'+str(first)+'-'+str(last)+' SHA256='+expected+'\n')
            excerpts.extend(str(i+1)+': '+lines[i]+'\n' for i in range(first-1,min(last,len(lines))))
    (out/'exact_source_excerpts.txt').write_text(''.join(excerpts),encoding='utf8')
    config = {}
    for line in config_path.read_text('utf8').splitlines():
        line = line.split('%',1)[0].strip()
        if '=' in line:
            key,value = line.split('=',1)
            config[key.strip()] = value.strip()
    required = dict(SOLVER='NAVIER_STOKES',KIND_TURB_MODEL='NONE',FLUID_MODEL='IDEAL_AIR_VARIABLE_CP',
                    INIT_OPTION='TD_CONDITIONS',NUM_METHOD_GRAD='WEIGHTED_LEAST_SQUARES',
                    CONV_NUM_METHOD_FLOW='ROE',MUSCL_FLOW='YES',MGLEVEL='0',AXISYMMETRIC='NO')
    for key,value in required.items():
        assert config[key] == value, (key,config.get(key))
    assert config['MARKER_SYM'] == '(sym_y, sym_z)'
    assert config['MARKER_ISOTHERMAL'] == '(wall, 335)'
    assert config['MARKER_PYTHON_CUSTOM'] == '(wall)'
    assert 'STATIC_SUPERSONIC_INFLOW_PT' in config['MARKER_RIEMANN'] and 'STATIC_PRESSURE' in config['MARKER_RIEMANN']
    assert set(k for k in config if k.startswith('MARKER_')) == {
        'MARKER_ISOTHERMAL','MARKER_SYM','MARKER_PLOTTING','MARKER_MONITORING','MARKER_RIEMANN','MARKER_PYTHON_CUSTOM'}
    assert config.get('FREESTREAM_OPTION','TEMPERATURE_FS') == 'TEMPERATURE_FS'
    assert 'FreeStreamOption_Map, FREESTREAM_OPTION::TEMPERATURE_FS' in texts[CONFIG]

    # 模块三：查询位置采用固定行号并验证原语句；这是准确源静态覆盖，不是旧异常调用栈。
    queries = []
    def add(relative,line,mode,active,role,protection,token):
        actual = texts[relative].splitlines()[line-1].strip()
        assert token in actual, (relative,line,actual)
        queries.append(dict(source=relative,line=line,mode=mode,active=active,role=role,
                            currentProtection=protection,exactStatement=actual))
    add(EULER,321,'rhoe',True,'构造期节点检查','物性域严格拒绝，无本地回退','SetTDState_rhoe(Density, StaticEnergy)')
    add(EULER,926,'PT',True,'给定来流初始化','固定输入，物性域严格拒绝','SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream)')
    add(EULER,931,'Prho',False,'密度形式来流初始化','当前TEMPERATURE_FS未走此分支','SetTDState_Prho')
    add(EULER,989,'rhoT',False,'Reynolds形式初始化','当前TD_CONDITIONS未走此分支','SetTDState_rhoT')
    add(EULER,1175,'Prho',True,'无量纲参考态能量','SetEnergy_Prho转发Prho；当前SI量纲保持','SetEnergy_Prho')
    add(CNS,153,'rhoe',True,'更新后当前单元状态','域拒绝先于原生check/旧态回退；无catch','SetTDState_rhoe(density, staticEnergy)')
    add(CNS,177,'rhoe',True,'原生非法状态回退后的旧态','仅其他检查能先到此处；旧态再次严格拒绝','SetTDState_rhoe(density, staticEnergy)')
    add(EULER,2073,'Prho',True,'MUSCL两侧重构一致性','当前模型domain_error触发两侧原始primitive+secondary配对回退','SetTDState_Prho(pressure, density)')
    add(EULER,2114,'Prho',False,'LowMach修正左侧能量','当前模型配置门禁止Low_Mach_Corr','SetEnergy_Prho')
    add(EULER,2117,'Prho',False,'LowMach修正右侧能量','当前模型配置门禁止Low_Mach_Corr','SetEnergy_Prho')
    add(EULER,5346,'rhoe',True,'RIEMANN内部单元','无本地保护；合法接受单元可间接保护同U','SetTDState_rhoe(Density_i, StaticEnergy_i)')
    add(EULER,5378,'PT',False,'RIEMANN总态输入','当前模型配置门禁止TOTAL_CONDITIONS_PT','SetTDState_PT')
    add(EULER,5394,'hs',False,'RIEMANN总态到静焓熵','当前模型配置门禁止TOTAL_CONDITIONS_PT','SetTDState_hs')
    add(EULER,5414,'PT',True,'RIEMANN固定静PT外场','固定输入严格拒绝，无边界回退','SetTDState_PT(P_static, T_static)')
    add(EULER,5441,'Prho',False,'RIEMANN固定静p-rho外场','当前模型配置门禁止STATIC_SUPERSONIC_INFLOW_PD','SetTDState_Prho')
    add(EULER,5478,'Prho',True,'RIEMANN给定p与内部rho组合','不由合法单元U保证；严格拒绝，保持给定压力','SetTDState_Prho(Pressure_e, Density_e)')
    add(EULER,5548,'rhoe',True,'RIEMANN特征混合边界U','两侧分别合法不保证混合合法；无本地保护','SetTDState_rhoe(Density_b, StaticEnergy_b)')

    state = texts[CNS]
    assert 'catch' not in state
    assert state.index('FluidModel->SetTDState_rhoe(density, staticEnergy);') < state.index('if (check_dens || check_press')
    muscl = section(texts[EULER],'void CEulerSolver::Upwind_Residual','void CEulerSolver::ComputeConsistentExtrapolation')
    for token in ('catch (const std::domain_error&)', 'if (!variable_air) throw;',
                  'std::copy_n(V_i,nPrimVar,Primitive_i); std::copy_n(V_j,nPrimVar,Primitive_j);',
                  'std::copy_n(S_i,nSecondaryVar,Secondary_i); std::copy_n(S_j,nSecondaryVar,Secondary_j);',
                  'bad_recon = nodes->UpdateNonPhysicalEdgeCounter(iEdge, bad_recon);',
                  'numerics->SetPrimitive(bad_recon? V_i : Primitive_i,  bad_recon? V_j : Primitive_j);',
                  'numerics->SetSecondary(bad_recon? S_i : Secondary_i,  bad_recon? S_j : Secondary_j);'):
        assert token in muscl, token
    boundary = section(texts[EULER],'void CEulerSolver::BC_Riemann','void CEulerSolver::BC_TurboRiemann')
    assert 'catch' not in boundary
    assert 'u_b[iVar] += P_Tensor[iVar][jVar]*dw[jVar];' in boundary
    assert 'Low_Mach_Corr || nMarker_Giles != 0 || GetDynamic_Grid()' in texts[CONFIG]
    assert 'Kind_Data_Riemann[i] != STATIC_SUPERSONIC_INFLOW_PT && Kind_Data_Riemann[i] != STATIC_PRESSURE' in texts[CONFIG]
    body_without_queries = {
        'spatialIsothermalWall': section(texts[WALL],'void CNSSolver::BC_Isothermal_Wall_Generic(',
                                        'void CNSSolver::BC_Isothermal_Wall('),
        'symmetry': section(texts[BASE],'::BC_Sym_Plane(', '::BC_Periodic('),
        'generalRoe': texts[ROE][texts[ROE].index('CNumerics::ResidualType<> CUpwGeneralRoe_Flow::ComputeResidual'):],
        'primitiveOutput': texts['SU2_CFD/src/output/CFlowCompOutput.cpp'],
    }
    for name,body in body_without_queries.items():
        assert not QUERY.search(body), name
    assert 'nodes->SetVelocity_Old(iPoint, zero);' in body_without_queries['spatialIsothermalWall']
    assert 'solutionOld[iVel + iDim] -= vp * UnitNormal[iDim];' in body_without_queries['symmetry']
    assert 'GetCustomBoundaryTemperature(val_marker, iVertex)' in body_without_queries['spatialIsothermalWall']
    assert 'su2double dTdn = -(There - Twall)/dist_ij;' in body_without_queries['spatialIsothermalWall']
    assert 'if (RoeSoundSpeed2 <= 0.0)' in body_without_queries['generalRoe']
    assert 'new CUpwGeneralRoe_Flow' in texts['SU2_CFD/src/drivers/CDriver.cpp']
    model = texts['SU2_CFD/include/fluid/CSU2IdealAir.hpp']
    assert 'assign(rah_air::state("rhoe",rho,e))' in model
    assert 'assign(rah_air::state("Ps",p,s))' in model
    assert 'minimum_temperature = 60.0' in texts['SU2_CFD/include/fluid/IdealAir.hpp']
    assert 'maximum_temperature = 2000.0' in texts['SU2_CFD/include/fluid/IdealAir.hpp']
    stable = {path:sha(Path(path))==value for path,value in locks.items()}
    assert all(stable.values())
    report = dict(passed=True,reviewKind='finite-static-exact-source-query-coverage',
                  runtimeFlavor=runtime['flavor'],sourceCommit=runtime['sourceCommit'],runtimeSHA256=EXE_SHA,
                  sourceTreeSHA256=TREE_SHA,sourceBindings=sources,actualConfiguration=config,queries=queries,
                  activeQueryModes=sorted({q['mode'] for q in queries if q['active']}),
                  activeHsOrPsIdentified=False,queryFreeReviewedBodies=list(body_without_queries),
                  existingMusclPairedFallbackVerified=True,cellDomainExceptionPreemptsNativeFallback=True,
                  boundaryPressureDensityCombinationProtected=False,characteristicMixedStateProtected=False,
                  oldStateMayBeModifiedByNativeNoSlipOrSymmetry=True,
                  selectedCombination=dict(cellPath='same-conservative-direction admissible step or strict rejection',
                      reconstruction='retain exact paired fallback',boundary='strict rejection with context; prescribed pressure unchanged',
                      queryContext='mode and two inputs; phase and correctly identified node if available',
                      completeAdmissibilityCoverage=False,finitePossiblyRejectingCandidateSupported=True),
                  previousV172ExceptionActualCallerRecovered=False,allProgramPathsFormallyProven=False,
                  staleModelFieldsMayBeUsedAfterFailedQuery=False,sourceSHA256=sha(Path(__file__).resolve()),
                  inputSHA256=locks,inputsStable=stable,nativeCFDRuns=0,builds=0,sharedSourceMutations=0,
                  newPhysicalQualification=False,elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report)
    print(json.dumps(dict(passed=True,elapsedSeconds=report['elapsedSeconds'],sources=len(sources),
                         queries=len(queries),activeModes=report['activeQueryModes'],CFD=0,builds=0),ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args = parser.parse_args()
    main(args.output)
