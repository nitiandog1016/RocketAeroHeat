"""模块用途：原生定比热NACA空间壁温入口；保持几何，重新建立明确能量初态。"""
import argparse
import os
from pathlib import Path
import numpy as np

from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save
from naca3344Model import DATA, wall_temperature
from nacaConstantCpModel import MODEL, GAMMA, R, CP, CV, conditions, state
from prepareNACASpatialCase import fixture_geometry, expected_table, BASE, SPATIAL
from prepareSU2NACALaminar import write_mesh, config_text, capacity

KIND = 'naca-tn3344-constant-cp-spatial-wall-v1'
FLAVOR = 'rah-ideal-air-spatial-wall-v1'
PLAN = ROOT/'docs/validation/v174/PLAN.md'
RUNTIME_SHA = 'e794985c282f73913c679eb72ad1818c9d56e97371813270ec7364865bfa3b3d'
EXE_SHA = '3309e77ed0fdfdb465e6d774f2157b7ef8d2cb2252a3d4abbfab7ee49fd7a217'
STAGE = 'same-flow-update-laminar-native-ideal-gas-v1'
PREREQUISITES = [ROOT/'outputs/v174_validation'/name/'report.json' for name in
                 ('properties01','runtime_review01','independent_anchor01')]


def specification(variant):
    """模块一：仅开放已事前登记的小组合；大候选须另行定义明确的总预算。"""
    assert variant == 'small-profile', '实际发展候选尚未登记'
    return dict(variant=variant, actualValidationMesh=False, updates=4, every=2,
        maximumSeconds=30, threads=2, maximumOutputBytes=32*1024**2,
        requiredRuntimeFlavor=FLAVOR, thermalMode='axialTemperature')


def configuration(spec, c):
    """模块二：保留Riemann静PT/静压、Roe/WLS/二阶，显式切换原生IDEAL_GAS。"""
    text = config_text(c, dict(walls=[dict(marker='wall',givenTw_K=335.)]),
                       spec['updates'], spec['every'], True)
    assert text.count('FLUID_MODEL= IDEAL_AIR_VARIABLE_CP') == 1
    text = text.replace('FLUID_MODEL= IDEAL_AIR_VARIABLE_CP',
                        'FLUID_MODEL= IDEAL_GAS\nGAMMA_VALUE= %.17g' % GAMMA)
    return text + 'MARKER_PYTHON_CUSTOM= (wall)\n'


def initial_state(points, meta, labels, c):
    """模块三：解析C1物理初猜。新e=cvT，不读取任何旧OTH守恒态。"""
    s = np.asarray(meta['s'])[labels[:,1]]
    distance = np.asarray(meta['offsets'])[labels[:,0]]
    tw = wall_temperature(s)[0]
    eta = np.minimum(distance/.001,1.)
    b = eta*eta*(3-2*eta)
    t0 = c['totalTemperature_K']; ti = c['staticTemperature_K']
    temperature = (1-b)*tw+b*t0-b*b*(t0-ti)
    rho = c['pressure_Pa']/(R*temperature)
    velocity = np.column_stack((-c['speed_m_s']*b,np.zeros((len(points),2))))
    conservative = np.column_stack((rho,rho[:,None]*velocity,
        rho*(CV*temperature+.5*np.sum(velocity*velocity,axis=1))))
    assert np.isfinite(conservative).all() and np.all(rho>0)
    assert temperature.min() >= ti-1e-12 and temperature.max() <= tw.max()+1e-12
    return conservative, temperature


def checkpoint_plan(spec, nodes, wall_nodes, wall_faces):
    value = capacity(nodes,wall_nodes,wall_faces,spec['updates'],spec['every'],spec['maximumOutputBytes'])
    value['stage'] = STAGE
    return value


def source_files():
    """模块四：新模型/准备/审核与共享运行分派同时锁定，任务中不能修改。"""
    names = ('prepareNACAConstantCpCase.py','su2NACAConstantCpInterface.py',
        'nacaConstantCpModel.py','auditNACAConstantCp.py','prepareNACASpatialCase.py',
        'prepareSU2NACALaminar.py','naca3344Model.py','nacaSpatialWallFlux.py',
        'nacaWallFlux.py','su2NACALaminarInterface.py','su2IdealAirCheckpoints.py',
        'su2ThermalRuntime.py','runEmbeddedSU2.py','su2IdealAirInterface.py')
    return [ROOT/'scripts'/name for name in names]+[PLAN]+[DATA/name for name in
        ('contract.json','source_manifest.json','wall_temperature.csv','heat_parameter.csv')]+PREREQUISITES


def prepare(folder, variant):
    """模块五：独占目录，生成全新9列初场和独立ideal_gas合同，不改父运行包。"""
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    spec=specification(variant)
    assert 'V174_SMALL_NATIVE_IDEAL_GAS_REGISTERED' in PLAN.read_text(encoding='utf-8'), '小任务尚未登记'
    for path in PREREQUISITES:
        proof=read(path);assert proof['passed'] and all(proof['inputsStable'].values())
        for name,digest in proof['inputSHA256'].items():
            # 第一阶段PLAN精确快照已留存；只允许当前事前登记追加，不替换其他输入。
            if Path(name).resolve()!=PLAN.resolve():assert sha(Path(name))==digest
    assert read(PREREQUISITES[0])['developmentEntryEligible'] is True
    locks={str(p):sha(p) for p in source_files()}
    folder.mkdir(parents=True,exist_ok=False)
    (folder/'executed_prepare.py').write_bytes(Path(__file__).read_bytes())
    (folder/'plan_at_prepare.md').write_bytes(PLAN.read_bytes())
    p,cells,markers,meta,labels,ix=fixture_geometry()
    write_mesh(folder/'mesh.su2',p,cells,markers)
    np.savez_compressed(str(folder/'mesh_arrays.npz'),points=p,nodeIndices=labels,columnNodeIDs=ix,
        prisms=np.array([c[1:] for c in cells if c[0]==13]),hexes=np.array([c[1:] for c in cells if c[0]==12]))
    save(folder/'geometry.json',meta)
    c=conditions();initial,temperature=initial_state(p,meta,labels,c)
    columns=('PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy')
    np.savetxt(str(folder/'initial_flow.csv'),np.column_stack((np.arange(len(p)),p,initial)),
        delimiter=',',fmt='%.17g',comments='',header=','.join('"'+n+'"' for n in columns))
    (folder/'case.cfg').write_text(configuration(spec,c),encoding='ascii')
    np.savetxt(str(folder/'wall_temperature.csv'),expected_table('spatial-profile',meta),
        delimiter=',',fmt='%.17g',comments='',header=SPATIAL['schema'])
    wallids=np.unique([i for row in markers['wall'] for i in row[1:]])
    contract=dict(schema='RAH.NACA.NativeIdealGas.v1',initialization=KIND,
        selectedFluidModel='IDEAL_GAS',thermodynamicModel=MODEL,gamma=GAMMA,gasConstant=R,
        energyReference='e=cv*T; h=cp*T; no-OTH-offset',turbulent=False,
        developmentTemperatureRangeK=[150.,360.],physicalQualified=False,newQualifiedArea_m2=0,
        configurationSHA256=sha(folder/'case.cfg'),meshSHA256=sha(folder/'mesh.su2'),
        geometrySHA256=sha(folder/'geometry.json'),meshArraysSHA256=sha(folder/'mesh_arrays.npz'),
        checkpointContract=checkpoint_plan(spec,len(p),len(wallids),len(markers['wall'])),
        reusedRuntimeMetadataIsSelectedModel=False)
    save(folder/'ideal_gas_contract.json',contract)
    manifest=dict(caseKind=KIND,specification=spec,requiredRuntimeFlavor=FLAVOR,
        selectedFluidModel='IDEAL_GAS',idealGas=contract,idealGasContractSHA256=sha(folder/'ideal_gas_contract.json'),
        runtimeUse='native-ideal-gas-path-in-accurate-v171-spatial-wall-binary',
        runtimeIdentity=dict(runtimeSHA256=RUNTIME_SHA,executableSHA256=EXE_SHA),
        thermalBoundary=dict(mode='axialTemperature',givenTw_K=335.,spatialWall=SPATIAL),
        thermalProfileSHA256=sha(folder/'wall_temperature.csv'),
        publicBenchmark=dict(conditions=c,fixture=True,experimentalQualification=False),
        initialization=dict(kind=KIND,priorSuccessfulUpdates=0,initialStateSHA256=sha(folder/'initial_flow.csv'),
            sourceSHA256=locks,initialBlendThickness_m=.001,boundaryLayerSolved=False,
            energyReference=contract['energyReference'],reusesPreviousConservativeState=False,
            temperatureRangeK=[float(temperature.min()),float(temperature.max())]),
        nodes=len(p),elements=len(cells),expectedVelocity=[-c['speed_m_s'],0,0],
        physicalQualified=False,planSHA256=sha(PLAN))
    save(folder/'case_manifest.json',manifest)
    save(folder/'request.json',dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        maximumSeconds=spec['maximumSeconds'],threads=spec['threads'],
        sourceFiles=[p.relative_to(ROOT).as_posix() for p in source_files()]))
    from su2ThermalRuntime import resolve
    resolve(folder,(folder/'case.cfg').read_text())
    assert all(sha(Path(name))==value for name,value in locks.items())
    print(dict(prepared=True,variant=variant,nodes=len(p),cells=len(cells),selectedModel=MODEL))
    return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',required=True)
    parser.add_argument('--variant',choices=['small-profile'],required=True)
    args=parser.parse_args();prepare(args.folder,args.variant)
