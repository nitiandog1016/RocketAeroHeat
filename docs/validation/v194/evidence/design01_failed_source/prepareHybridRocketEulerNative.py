# -*- coding: utf-8 -*-
"""真实原全壁Euler的具名准备与审核；只准备初场，不启动发展或授予热资格。"""
from pathlib import Path
import argparse
import json
import math
import shutil
import time
import numpy as np
from readEilmerNative import sha, need, read_grid, read_field, compare_source
from runEilmerProcess import COMMIT
from eilmerAirThermo import bind_model, evaluate
from hybridEdgeSampling import map_original_wall

ROOT = Path(__file__).resolve().parents[1]
PROFILE = 'original-wall-rocket-euler-01'
SCHEMA = 'RAH.HybridRocketEulerPreparation.v1'


def local(path):
    """模块一：只转换本项目E盘在Windows与WSL间的已知路径。"""
    text = str(path).replace('\\', '/')
    return Path(text.replace('E:/', '/mnt/e/', 1) if Path('/mnt/e').exists()
                else text.replace('/mnt/e/', 'E:/', 1))


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def verify_design(source):
    """核对网格接受证据和固定工况；既有小控制/旧黏性域不能冒用此合同。"""
    source = Path(source)
    m = json.loads((source/'manifest.json').read_text())
    need(m['schema'] == SCHEMA and m['profile'] == PROFILE and m['sourceCommit'] == COMMIT,
         '非具名真实Euler准备合同')
    need(m['executionAuthorized'] is False and m['nativeConfig']['max_step'] == 0,
         '准备合同不能授权发展')
    for name, digest in m['filesSHA256'].items():
        need(sha(source/name) == digest, '设计文件改变: '+name)
    for name, digest in m['sourceSHA256'].items():
        need(sha(local(name)) == digest, '来源证据改变: '+name)
    c = m['condition']
    need((c['Mach'], c['alphaDeg'], c['betaDeg'], c['pressurePa'], c['temperatureK']) ==
         (5., 20., 0., 1197.0279007383137, 226.50908361133003), '目标工况改变')
    g = m['geometry']
    need((g['cells'], g['vertices'], g['wallTriangles'], g['farfieldTriangles']) ==
         (64982, 18027, 19752, 6572), '非已审核原全壁网格')
    need(abs(g['wallAreaM2']-1.9646292565715551) < 2e-13, '原面积改变')
    return m


def design(mesh, proof, out):
    """模块二：复用已接受网格和原大气，只生成新的Euler输入，不修改旧资产。"""
    mesh, proof, out = map(Path, (mesh, proof, out))
    need(not out.exists(), '禁止覆盖设计目录')
    r = json.loads((mesh/'report.json').read_text())
    a = json.loads(proof.read_text())
    need(r['passed'] and r['qualityScreenPassed'] and r['originalWallExactlyPreserved'] and
         r['originalFarfieldExactlyPreserved'] and a['passed'] and a['geometryQualityAccepted'],
         '网格没有同时通过几何/质量')
    need(sha(mesh/'mesh.npz') in a['inputsSHA256'].values() and
         sha(mesh/'report.json') in a['inputsSHA256'].values(), '独立审核未绑定本候选')
    for name, digest in r['filesSHA256'].items():
        need(sha(mesh/name) == digest, '已审核网格变化')
    old_path = ROOT/'outputs/v178_validation/real01/case_manifest.json'
    old = json.loads(old_path.read_text())['condition']
    c = {k: old[k] for k in ['pressurePa', 'temperatureK', 'Mach', 'alphaDeg', 'betaDeg', 'wallTemperatureK']}
    alpha, beta = map(math.radians, (c['alphaDeg'], c['betaDeg']))
    c.update(altitudeM=30000., velocityDirection=[-math.cos(alpha)*math.cos(beta), -math.sin(beta),
                                               math.sin(alpha)*math.cos(beta)])
    out.mkdir(parents=True)
    for name in ['mesh.su2', 'surface.json']:
        shutil.copyfile(mesh/name, out/name)
    (out/'gas.inp').write_text('model = "ThermallyPerfectGas"\nspecies = {"air"}\n')
    (out/'grid.lua').write_text('''-- 原19752面全部保留；只读取新的全域Euler体网格。
config.dimensions=3
config.grid_format="rawbinary"
registerFluidGrid{grid=UnstructuredGrid:new{filename="mesh.su2",fmt="su2text",
 scale=1.0,expect_gmsh_order_for_wedges=false},fsTag="initial"}
''', encoding='utf-8')
    (out/'job.lua').write_text('''-- 模块一：原生CEA空气和真实20度来流；壁温仅为后续热目标，不施加给Euler。
config.dimensions=3
config.solver_mode="time-marching"
config.field_format="rawbinary"
config.grid_format="rawbinary"
config.viscous=false
config.reacting=false
config.turbulence_model="none"
config.adjust_invalid_cell_data=false
config.max_invalid_cells=0
local nsp,nmodes,gm=setGasModel("rah-air.gas")
local gs=GasState:new{gm}
gs.p=%.17g; gs.T=%.17g; gs.massf={air=1.0}
gm:updateThermoFromPT(gs); gm:updateSoundSpeed(gs); gm:updateTransCoeffs(gs)
local speed=5.0*gs.a
local alpha,beta=math.rad(20),math.rad(0)
local ux=-speed*math.cos(alpha)*math.cos(beta)
local uy=-speed*math.sin(beta)
local uz=speed*math.sin(alpha)*math.cos(beta)
local inflow=FlowState:new{p=gs.p,T=gs.T,massf=gs.massf,velx=ux,vely=uy,velz=uz}
local receipt=assert(io.open("native-input-state.json","w"))
receipt:write(string.format('{"p":%%.17g,"T":%%.17g,"rho":%%.17g,"e":%%.17g,"a":%%.17g,"cp":%%.17g,"velx":%%.17g,"vely":%%.17g,"velz":%%.17g}',
 gs.p,gs.T,gs.rho,gs.u,gs.a,gm:Cp(gs),ux,uy,uz)); receipt:close()
-- 模块二：成熟反射滑移壁及环境外边界；无热壁/湍流/局部截断。
makeFluidBlocks({wall=WallBC_WithSlip:new{group="wall"},
 farfield=InOutFlowBC_Ambient:new{flowState=inflow}},{initial=inflow})
-- 模块三：保持已核对的Euler算子；本轮max_step=0，准备不授予推进权限。
config.flux_calculator="ausmdv"
config.interpolation_order=2
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
config.gasdynamic_update_scheme="tvd_rk3"
config.with_local_time_stepping=true
config.local_time_stepping_limit_factor=1000
config.with_super_time_stepping=false
config.residual_smoothing=false
config.stringent_cfl=true
config.cfl_value=0.3
config.cfl_count=1
config.fixed_time_step=false
config.dt_init=1.0e-10
config.dt_max=1.0e-3
config.max_time=1.0
config.max_step=0
config.max_attempts_for_step=1
config.print_count=1
config.write_transient_residuals=true
config.dt_plot_value=1.0
''' % (c['pressurePa'], c['temperatureK']), encoding='utf-8')
    paths = [mesh/'report.json', mesh/'mesh.npz', proof, old_path, Path(__file__).resolve(),
             ROOT/'docs/validation/v194/PREPARATION_PLAN.md']
    m = dict(schema=SCHEMA, profile=PROFILE, sourceCommit=COMMIT, condition=c,
        executionAuthorized=False, physicalQualification=False, newSolvedHeatAreaM2=0,
        wallHeatFluxAvailable=False, temperatureIdentity='uniform initial fluid state; no solved wall temperature',
        geometry=dict(cells=r['cells'], vertices=r['vertices'], faces=r['faces'], ghostCells=r['ghostCells'],
                      wallTriangles=19752, farfieldTriangles=6572, wallAreaM2=r['originalAreaM2'],
                      volumeM3=a['volumeM3'], cellObjectsIncludingGhost=r['cellObjectsIncludingGhost'],
                      stageArrayPayloadBytes=r['stageArrayPayloadBytes']),
        nativeConfig=dict(max_step=0, conservedVariables=5, retainedTimeLevels=4),
        filesSHA256={n:sha(out/n) for n in ['mesh.su2','surface.json','gas.inp','grid.lua','job.lua']},
        sourceSHA256={str(p.resolve()):sha(p) for p in paths})
    save(out/'manifest.json', m)
    verify_design(out)
    return dict(passed=True, scope='design only', geometry=m['geometry'])


def prepare(source, case, record):
    """模块四：三项原生准备共用90秒预算；验证固定包，无构建和求解器启动。"""
    from runEilmerProcess import run
    from eilmerWorkspaceRuntime import audit_workspace_runtime
    start=time.monotonic()
    source, case, record = map(lambda p: Path(p).resolve(), (source, case, record))
    m=verify_design(source)
    need(not case.exists() and not record.exists(), '禁止覆盖原生目录')
    runtime=ROOT/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json'
    package, locks=audit_workspace_runtime(dict(explicit=dict(profile='rah-explicit-workspace-v1',
        runtimePackage=str(runtime), runtimePackageSHA256=sha(runtime))))
    case.mkdir(parents=True); record.mkdir(parents=True)
    for name in list(m['filesSHA256'])+['manifest.json']:
        shutil.copyfile(source/name, case/name)
    shutil.copyfile(__file__, record/'executed_prepare.py')
    locks += [local(p) for p in m['sourceSHA256']]
    locks += [case/n for n in list(m['filesSHA256'])+['manifest.json']]
    locks += [Path(__file__), Path(__file__).with_name('runEilmerProcess.py')]
    result=dict(schema='RAH.HybridRocketEulerNativePreparation.v1', profile=PROFILE, passed=False,
        case=str(case), sourceCommit=COMMIT, commands=[], runtimeSHA256=sha(runtime),
        executionAuthorized=False, physicalQualification=False, newSolvedHeatAreaM2=0)
    try:
        for i, command in enumerate([['prep-gas','gas.inp','rah-air.gas'],
                                    ['lmr','prep-grids','--job=grid.lua'],
                                    ['lmr','prep-sim','--job=job.lua']]):
            left=90-(time.monotonic()-start); need(left>0, '原生准备预算耗尽')
            receipt=run(command, case, record/('process%02d'%i), left, 128, locks, 4096, 4096)
            result['commands'].append(receipt)
            need(receipt['passed'], '原生准备失败')
        for name in ['rah-air.gas','native-input-state.json','lmrsim/config','lmrsim/control',
                     'lmrsim/grid/grid-0000','lmrsim/snapshots/fluid.metadata','lmrsim/snapshots/0000/fluid-0000']:
            need((case/name).is_file() and (case/name).stat().st_size>0, '缺少原生产物: '+name)
        result['nativeInputSHA256']={str(p):sha(p) for p in case.rglob('*') if p.is_file()}
        result.update(passed=True, scope='native preparation only; no flow development')
    except Exception as error:
        result['error']=repr(error)
    result['elapsedSeconds']=time.monotonic()-start
    save(record/'report.json', result)
    return result


def audit(case, receipt, out):
    """模块五：读取真实原生网格和初场，独立核对工况、原面、体积与采样位置。"""
    start=time.monotonic(); case, receipt, out=map(Path,(case,receipt,out))
    need(not out.exists(), '禁止覆盖审核目录'); out.mkdir(parents=True)
    shutil.copyfile(__file__, out/'executed_audit.py')
    m=verify_design(case); r=json.loads(receipt.read_text())
    need(r['schema']=='RAH.HybridRocketEulerNativePreparation.v1' and r['passed'], '原生准备未通过')
    need(local(r['case']).resolve()==case.resolve(), '原生目录身份不符')
    for path,digest in r['nativeInputSHA256'].items():
        need(sha(local(path))==digest, '原生准备文件改变')
    need(len(r['commands'])==3 and all(p['passed'] and p['nativeExitCode']==0 and p['inputsStable']
         and not p['liveGroupAfterCleanup'] for p in r['commands']), '原生准备退出/来源不完整')
    cfg=json.loads((case/'lmrsim/config').read_text()); ctl=json.loads((case/'lmrsim/control').read_text())
    expected=dict(dimensions=3, axisymmetric=False, viscous=False, reacting=False, turbulence_model='none',
        solver_mode='time-marching', grid_format='rawbinary', field_format='rawbinary',
        flux_calculator='ausmdv', interpolation_order=2, gasdynamic_update_scheme='tvd_rk3',
        spatial_deriv_calc='least_squares', spatial_deriv_locn='cells', with_local_time_stepping=True,
        with_super_time_stepping=False, residual_smoothing=False, adjust_invalid_cell_data=False, max_invalid_cells=0)
    for key,value in expected.items(): need(cfg.get(key)==value, '原生算子不符: '+key)
    need(cfg['cfl_schedule_values']==[.3] and ctl['max_step']==0, '准备阶段控制改变')
    grid=read_grid(case/'lmrsim/grid/grid-0000'); mapping, error=compare_source(grid,case/'mesh.su2')
    g=m['geometry']; need((len(grid['cells']),len(grid['points']),len(grid['faces']))==
                         (g['cells'],g['vertices'],g['faces']), '原生规模改变')
    need(set(grid['boundaries'])=={'wall','farfield'}, '边界标签改变')
    for name,kind,count in [('wall','wall_with_slip',19752),('farfield','inoutflow_ambient',6572)]:
        b=grid['boundaries'][name]; bc=cfg['block_0']['boundary_%d'%b['index']]
        need(len(b['ids'])==count and bc['type']==kind and not bc['is_wall_with_viscous_effects'],
             '原生边界不符: '+name)
    f=read_field(case,0,g['cells']); native=json.loads((case/'native-input-state.json').read_text())
    model=bind_model(case/'rah-air.gas',sha(case/'rah-air.gas')); thermo=evaluate(f['p'],f['T'])
    eos={}
    for name,key in [('rho','densityKgm3'),('e','internalEnergyJkg'),('a','soundSpeedMps')]:
        eos[name]=float(np.max(abs(f[name]-thermo[key])/np.maximum(1,abs(thermo[key]))))
        need(eos[name]<2e-10, '原生EOS不相容: '+name)
    speed=np.linalg.norm([native[k] for k in ['velx','vely','velz']]); direction=np.array(m['condition']['velocityDirection'])
    need(abs(speed/native['a']-5)<2e-13 and np.max(abs(np.array([native[k] for k in ['velx','vely','velz']])/speed-direction))<2e-14,
         'Mach或攻角方向不符')
    need(native['p']==m['condition']['pressurePa'] and native['T']==m['condition']['temperatureK'], '大气输入改变')
    for field,key in [('p','p'),('T','T'),('rho','rho'),('e','e'),('a','a'),('vel.x','velx'),('vel.y','vely'),('vel.z','velz')]:
        need(np.max(abs(f[field]-native[key]))<2e-13*max(1,abs(native[key])), '初场不均匀: '+field)
    need(np.all(f['massf-air']==1) and np.all(f['vol']>0), '初场组分或体积错误')
    volume=float(f['vol'].sum()); need(abs(volume/g['volumeM3']-1)<2e-10, '原生总体积改变')
    surface=json.loads((case/'surface.json').read_text()); positions=np.c_[[f['pos.'+k] for k in 'xyz']].T
    wall=map_original_wall(surface,grid,mapping['wall'],positions)
    # 独立读原网格owner及几何中心，不以映射函数本身作真值。
    mesh_paths=[local(p) for p in m['sourceSHA256'] if p.endswith('mesh.npz')]
    need(len(mesh_paths)==1,'原网格来源歧义'); mesh=np.load(mesh_paths[0])
    owners=np.array(wall['nativeOwnerCellsZeroBased']); need(np.array_equal(owners,mesh['wallOwners']), '原面owner改变')
    np.testing.assert_allclose(positions,mesh['vertices'][mesh['tetrahedra']].mean(1),rtol=0,atol=2e-14)
    tri=mesh['vertices'][mesh['wall']]; vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    area=np.linalg.norm(vec,axis=1)/2
    distance=np.einsum('ij,ij->i',positions[owners]-tri.mean(1),vec/(2*area[:,None]))
    need(np.max(abs(distance-wall['normalDistanceM']))<2e-14, '原面法向距离改变')
    need(abs(area.sum()-g['wallAreaM2'])<2e-13,'原面积改变')
    np.savetxt(out/'original_wall_initial_sampling.csv',np.c_[np.arange(1,19753),owners,area,distance,
        positions[owners],f['p'][owners],f['T'][owners],f['vel.x'][owners],f['vel.y'][owners],f['vel.z'][owners]],
        delimiter=',',comments='',fmt='%.17g',header='face_id,owner_zero_based,area_m2,normal_distance_m,x_m,y_m,z_m,p_Pa,T_fluid_K,Ux_mps,Uy_mps,Uz_mps')
    result=dict(passed=True, profile=PROFILE, geometry=g, actualNumerics=expected, actualGasModel=model,
        sourceCoordinateError=error, originalFaceOwnerExact=True, maximumNormalDistanceErrorM=float(np.max(abs(distance-wall['normalDistanceM']))),
        nativeVolumeM3=volume, nativeMach=float(speed/native['a']), nativeVelocityMps=[native[k] for k in ['velx','vely','velz']],
        fullVolumeEOSErrors=eos, initialTotalEnthalpyJkg=float(thermo['enthalpyJkg'][0]+speed**2/2),
        executionAuthorized=False, physicalQualification=False, newSolvedHeatAreaM2=0, solverPeakMemoryBytes=None,
        note='uniform initial field, not a computed external-flow or thermal result',
        sourceSHA256={str(p.resolve()):sha(p) for p in [receipt,Path(__file__),case/'manifest.json',case/'lmrsim/config']},
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(); sub=parser.add_subparsers(dest='mode',required=True)
    for name,keys in [('design',['mesh','proof','out']),('prepare',['source','case','record']),('audit',['case','receipt','out'])]:
        command=sub.add_parser(name)
        for key in keys: command.add_argument('--'+key, required=True)
    args=vars(parser.parse_args()); mode=args.pop('mode'); result=globals()[mode](**args)
    print(json.dumps({k:result.get(k) for k in ['passed','error','elapsedSeconds','geometry']}))
    raise SystemExit(0 if result['passed'] else 1)
