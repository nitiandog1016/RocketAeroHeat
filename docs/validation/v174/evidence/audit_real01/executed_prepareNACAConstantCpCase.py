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
from prepareSU2NativeGasOutputPatch import OUTPUT_CONTRACT

KIND = 'naca-tn3344-constant-cp-spatial-wall-v1'
FLAVOR = 'rah-native-ideal-gas-spatial-wall-v1'
PLAN = ROOT/'docs/validation/v174/PLAN.md'
RUNTIME_SHA = 'faaed7f318876335f5309fb79bbbc0facb2d50a9972bfaaa94b97220f9457ebb'
EXE_SHA = 'dd092a7317568bf3ed1edefa2e4a7e09e9385867a377e94d4097fb9196b7fe8e'
RUNTIME_USE = 'native-ideal-gas-v171-parent-with-explicit-first-completed-output-extension'
STAGE = 'same-flow-update-laminar-native-ideal-gas-v1'
PREREQUISITES = [ROOT/'outputs/v174_validation'/name/'report.json' for name in
                 ('properties01','runtime_review01','independent_anchor01')]
HISTORICAL_DISPATCH = ROOT/'outputs/v174_validation/small02_source_snapshot/project_sources/scripts/su2ThermalRuntime.py'
ACTUAL_PROOF = ROOT/'docs/validation/v170/evidence_manifest.json'
ACTUAL_PROOF_SHA = '975599db91d573a26e480f82bf5a1214fcd4d0ee40a3ef46f693c38c91c9ca0e'
ACTUAL_FILES = dict(
    **{'mesh.su2':'5562d7e6a3985ec6877585db9984d273a409178e54b1fb1d336065c93cd77152',
       'mesh_arrays.npz':'4e82f33ae66d0a19bf1a755d6c0d9210ca5cadd8db977ab5377da438517e3c01',
       'geometry.json':'62cb540c4a79d214434121a9f44ad9537de2258d006bcd808cd92ce814ba966b',
       'wall_temperature.csv':'d9b4e8916609fcb64ce4e95c10badd3793975da614d8b0f696bec3004d5ee660'})
ACTUAL_MESH_REPORT_SHA = '3ee9c90182abe53d6c16f53edb748835d1ff9b8e4b317ef2f9a259f9a9801aac'
ACTUAL_PREREQUISITES = {
    ROOT/'outputs/v174_validation/audit_small03/report.json':
        'e0651feb2f426674ba6563abe5e76b1a9cba906030d86225a084e346aaeeb3f0',
    ROOT/'outputs/v174_validation/compare_small03/report.json':
        'b048d99bd15c53dff77808cdb4d9e3d3619311fbc7d256a283b910aeabeb7bdf'}
GUARD_REPORT = ROOT/'outputs/v174_validation/domain_guard03/report.json'
GUARD_PLAN = ROOT/'docs/validation/v174/DOMAIN_GUARD_PLAN.md'
ACTUAL_RESOURCE_PROFILE = 'v174-native-ideal-gas-development'
VARIANTS = ('small-profile','validation-profile')


def specification(variant):
    """模块一：仅开放已事前登记的小组合；大候选须另行定义明确的总预算。"""
    assert variant in VARIANTS, '未登记的定比热候选'
    actual = variant == 'validation-profile'
    return dict(variant=variant, actualValidationMesh=actual, updates=500 if actual else 4, every=100 if actual else 2,
        maximumSeconds=600 if actual else 30, threads=2, maximumOutputBytes=(1536 if actual else 32)*1024**2,
        requiredRuntimeFlavor=FLAVOR, thermalMode='axialTemperature')


def validate_actual_mesh_files(folder=BASE):
    """真实几何模块：只复用v170准确四文件与原报告，不重建段界、节点或Tw表。"""
    folder=Path(folder)
    assert sha(ACTUAL_PROOF)==ACTUAL_PROOF_SHA
    proof=read(ACTUAL_PROOF)
    identified={Path(row['path']).name:row['sha256'] for row in proof['localOnly']
                if Path(row['path']).parent.as_posix()=='outputs/v170_validation/mesh02'}
    identified.update({Path(row['source']).name:row['sha256'] for row in proof['files']
        if Path(row.get('source','')).parent.as_posix()=='outputs/v170_validation/mesh02'})
    assert all(identified[name]==digest for name,digest in ACTUAL_FILES.items())
    assert identified['report.json']==ACTUAL_MESH_REPORT_SHA==sha(BASE/'report.json')
    for name,digest in ACTUAL_FILES.items():
        assert sha(BASE/name)==digest and sha(folder/name)==digest, '真实网格或壁温表身份变化: '+name
    report=read(BASE/'report.json')
    assert report['geometryCandidatePrepared'] and all(report['inputsStable'].values())
    assert report['nodes']==190334 and report['cells']==179568 and report['wallFaces']==3096
    assert report['wallTemperature']['nodes']==130 and not report['physicalQualified']
    return dict(source='v170-mesh02-exact-four-files',evidenceManifestSHA256=ACTUAL_PROOF_SHA,
        fileSHA256=dict(ACTUAL_FILES),meshReportSHA256=ACTUAL_MESH_REPORT_SHA,
        nodes=190334,elements=179568,wallFaces=3096,wallTemperatureStations=130,
        previousFailedConservativeStateReused=False,physicalQualified=False)


def actual_prerequisite_summary():
    """研发前提模块：承认粗网格温域失败，只继承已验证的接口能力及守卫小检查。"""
    for path,digest in ACTUAL_PREREQUISITES.items():assert sha(path)==digest
    audit,compare=(read(path) for path in ACTUAL_PREREQUISITES)
    assert audit['interfacePassed'] and audit['normalRunCompleted'] and all(audit['inputsStable'].values())
    assert compare['comparisonPassed'] and compare['exactOutputEqualityPassed'] and all(compare['inputsStable'].values())
    assert compare['currentCompleteSavedLocalIterations']==[0,2,3]
    assert audit['allSavedFramesWithinDevelopmentModelDomain'] is False
    assert compare['allSavedFramesWithinDevelopmentModelDomain'] is False
    guard=read(GUARD_REPORT)
    assert guard['passed'] is True and guard['inputsStable'] and all(guard['inputsStable'].values())
    for path in (ROOT/'scripts/nacaConstantCpDomainGuard.py',GUARD_PLAN):
        assert guard['inputSHA256'][str(path)]==sha(path), '守卫前提执行源已变化'
    return dict(thermalInterfacePassed=True,sharedSavedFieldsExactlyEqual=True,
        coarseFixtureDevelopmentDomainPassed=False,developmentTemperatureRangeK=[150.,360.],
        basis='fresh-initial-state-on-original-validation-mesh; not-coarse-field-continuation',
        evidenceSHA256={str(p):digest for p,digest in ACTUAL_PREREQUISITES.items()},
        savedFieldGuardReportSHA256=sha(GUARD_REPORT),physicalQualified=False)


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


def source_files(variant='small-profile'):
    """模块四：新模型/准备/审核与共享运行分派同时锁定，任务中不能修改。"""
    names = ('prepareNACAConstantCpCase.py','su2NACAConstantCpInterface.py',
        'nacaConstantCpModel.py','auditNACAConstantCp.py','prepareNACASpatialCase.py',
        'prepareSU2NACALaminar.py','naca3344Model.py','nacaSpatialWallFlux.py',
        'nacaWallFlux.py','su2NACALaminarInterface.py','su2IdealAirCheckpoints.py',
        'su2ThermalRuntime.py','runEmbeddedSU2.py','su2IdealAirInterface.py',
        'prepareSU2NativeGasOutputPatch.py','nacaSavedFrameStudy.py')
    files=[ROOT/'scripts'/name for name in names]+[PLAN]+[DATA/name for name in
        ('contract.json','source_manifest.json','wall_temperature.csv','heat_parameter.csv')]+PREREQUISITES+[HISTORICAL_DISPATCH]
    if specification(variant)['actualValidationMesh']:
        files += [ACTUAL_PROOF,BASE/'report.json']+[BASE/name for name in ACTUAL_FILES]
        files += list(ACTUAL_PREREQUISITES)+[ROOT/'scripts/nacaConstantCpDomainGuard.py',GUARD_PLAN,GUARD_REPORT]
    return files


def prepare(folder, variant):
    """模块五：独占目录，生成全新9列初场和独立ideal_gas合同，不改父运行包。"""
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    spec=specification(variant)
    token='V174_ACTUAL_NATIVE_IDEAL_GAS_REGISTERED' if spec['actualValidationMesh'] else 'V174_OUTPUT_SMALL03_REGISTERED'
    assert token in PLAN.read_text(encoding='utf-8'), '本次定比热任务尚未登记'
    for path in PREREQUISITES:
        proof=read(path);assert proof['passed'] and all(proof['inputsStable'].values())
        for name,digest in proof['inputSHA256'].items():
            # 第一阶段沿用当时准确分派源；仅本次已登记输出分派换用其精确快照。
            # 旧PLAN快照保持，其余物性/基准/算法输入仍逐项检查当前原文件。
            if Path(name).resolve()==PLAN.resolve():continue
            bound=HISTORICAL_DISPATCH if Path(name).resolve()==(ROOT/'scripts/su2ThermalRuntime.py').resolve() else Path(name)
            assert sha(bound)==digest
    assert read(PREREQUISITES[0])['developmentEntryEligible'] is True
    if spec['actualValidationMesh']:
        actual_identity=validate_actual_mesh_files()
        actual_prerequisites=actual_prerequisite_summary()
    locks={str(p):sha(p) for p in source_files(variant)}
    folder.mkdir(parents=True,exist_ok=False)
    (folder/'executed_prepare.py').write_bytes(Path(__file__).read_bytes())
    (folder/'plan_at_prepare.md').write_bytes(PLAN.read_bytes())
    if spec['actualValidationMesh']:
        from su2NACALaminarInterface import mesh
        for name in ACTUAL_FILES:os.link(str(BASE/name),str(folder/name))
        p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
        with np.load(str(folder/'mesh_arrays.npz')) as arrays:
            labels=arrays['nodeIndices'];ix=arrays['columnNodeIDs']
        assert len(p)==190334 and len(cells)==179568 and len(meta['s'])==130
    else:
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
    if not spec['actualValidationMesh']:
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
        reusedRuntimeMetadataIsSelectedModel=False,
        inheritedPackageMetadataIsSelectedModel=False,selectedRuntimeMetadataIsSelectedModel=True)
    if spec['actualValidationMesh']:
        from nacaConstantCpDomainGuard import domain_guard_contract
        contract['savedFieldDomainGuard']=domain_guard_contract()
    save(folder/'ideal_gas_contract.json',contract)
    manifest=dict(caseKind=KIND,specification=spec,requiredRuntimeFlavor=FLAVOR,
        selectedFluidModel='IDEAL_GAS',idealGas=contract,idealGasContractSHA256=sha(folder/'ideal_gas_contract.json'),
        runtimeUse=RUNTIME_USE,nativeIdealGasOutputContract=dict(OUTPUT_CONTRACT),
        runtimeIdentity=dict(runtimeSHA256=RUNTIME_SHA,executableSHA256=EXE_SHA),
        thermalBoundary=dict(mode='axialTemperature',givenTw_K=335.,spatialWall=SPATIAL),
        thermalProfileSHA256=sha(folder/'wall_temperature.csv'),
        publicBenchmark=dict(conditions=c,fixture=not spec['actualValidationMesh'],experimentalQualification=False),
        initialization=dict(kind=KIND,priorSuccessfulUpdates=0,initialStateSHA256=sha(folder/'initial_flow.csv'),
            sourceSHA256=locks,initialBlendThickness_m=.001,boundaryLayerSolved=False,
            energyReference=contract['energyReference'],reusesPreviousConservativeState=False,
            temperatureRangeK=[float(temperature.min()),float(temperature.max())]),
        nodes=len(p),elements=len(cells),expectedVelocity=[-c['speed_m_s'],0,0],
        physicalQualified=False,planSHA256=sha(PLAN))
    if spec['actualValidationMesh']:
        manifest.update(actualValidationMesh=actual_identity,developmentPrerequisites=actual_prerequisites)
    save(folder/'case_manifest.json',manifest)
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',
        maximumSeconds=spec['maximumSeconds'],threads=spec['threads'],
        sourceFiles=[p.relative_to(ROOT).as_posix() for p in source_files(variant)])
    if spec['actualValidationMesh']:request['resourceProfile']=ACTUAL_RESOURCE_PROFILE
    save(folder/'request.json',request)
    from su2ThermalRuntime import resolve
    resolve(folder,(folder/'case.cfg').read_text())
    assert all(sha(Path(name))==value for name,value in locks.items())
    print(dict(prepared=True,variant=variant,nodes=len(p),cells=len(cells),selectedModel=MODEL))
    return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',required=True)
    parser.add_argument('--variant',choices=VARIANTS,required=True)
    args=parser.parse_args();prepare(args.folder,args.variant)
