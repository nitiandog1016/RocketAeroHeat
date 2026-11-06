"""模块用途：具名原生定比热入口的严格来源、物理合同和稀疏保存审核。"""
import csv
import re
from pathlib import Path
import numpy as np

from su2ThermalRuntime import ROOT,read,sha
from prepareNACAConstantCpCase import (KIND,FLAVOR,PLAN,RUNTIME_SHA,EXE_SHA,STAGE,RUNTIME_USE,OUTPUT_CONTRACT,
    specification,configuration,checkpoint_plan,initial_state,source_files,
    validate_actual_mesh_files,actual_prerequisite_summary,ACTUAL_RESOURCE_PROFILE)
from prepareNACASpatialCase import fixture_geometry,expected_table,SPATIAL
from nacaConstantCpModel import MODEL,GAMMA,R,CP,CV,conditions,state
from su2NACALaminarInterface import mesh


def validate_development_contract(manifest):
    """低成本合同模块：实际候选必须从新初态开始并启用原温区守卫，不能降级。"""
    assert manifest['caseKind']==KIND and manifest['requiredRuntimeFlavor']==FLAVOR
    spec=specification(manifest['specification']['variant']);assert manifest['specification']==spec
    contract=manifest['idealGas'];init=manifest['initialization']
    assert contract['selectedFluidModel']=='IDEAL_GAS' and contract['thermodynamicModel']==MODEL
    assert contract['developmentTemperatureRangeK']==[150.,360.]
    assert init['kind']==KIND and init['priorSuccessfulUpdates']==0
    assert init['reusesPreviousConservativeState'] is False and init['boundaryLayerSolved'] is False
    assert init['energyReference']=='e=cv*T; h=cp*T; no-OTH-offset'
    assert not any(key in init for key in ('sourceCheckpointIdentity','sourceLocalIteration','runTransition'))
    if spec['actualValidationMesh']:
        from nacaConstantCpDomainGuard import domain_guard_contract
        assert contract['savedFieldDomainGuard']==domain_guard_contract(), '真实候选缺失或改变模型温域守卫'
        assert manifest['nodes']==190334 and manifest['elements']==179568
        assert manifest['publicBenchmark']['fixture'] is False
        assert manifest['developmentPrerequisites']['coarseFixtureDevelopmentDomainPassed'] is False
    else:
        assert 'savedFieldDomainGuard' not in contract, '小例未登记运行中温域守卫'
        assert 'actualValidationMesh' not in manifest and 'developmentPrerequisites' not in manifest
    return spec


def inlet_state_check(text,c,expected_velocity):
    """模块一：独立由实际Riemann输入核对Mach向量、p/T和新能量基准。"""
    from su2IdealAirInterface import settings
    cfg=settings(text)
    assert cfg['FLUID_MODEL']=='IDEAL_GAS' and float(cfg['GAMMA_VALUE'])==GAMMA
    assert float(cfg['GAS_CONSTANT'])==R
    tokens=[v.strip() for v in cfg['MARKER_RIEMANN'].strip('() ').split(',')]
    assert len(tokens)==14 and tokens[:2]==['outer','STATIC_SUPERSONIC_INFLOW_PT']
    assert tokens[7:9]==['outlet','STATIC_PRESSURE']
    pressure,temperature=float(tokens[2]),float(tokens[3])
    assert pressure==c['pressure_Pa'] and temperature==c['staticTemperature_K']
    velocity=np.asarray(tokens[4:7],float)*np.sqrt(GAMMA*R*temperature)
    assert np.max(abs(velocity-np.asarray(expected_velocity)))<2e-10
    assert abs(np.linalg.norm(velocity)/np.sqrt(GAMMA*R*temperature)-c['Mach'])<2e-12
    assert float(tokens[9])==pressure and all(float(v)==0 for v in tokens[10:])
    total=CP*temperature+np.dot(velocity,velocity)/2
    error=abs(total/c['totalEnthalpy_Jkg']-1)
    assert error<2e-12
    return dict(velocityMS=velocity.tolist(),totalEnthalpyRelativeError=float(error),
        selectedFluidModel='IDEAL_GAS',energyReference='e=cv*T; h=cp*T; no-OTH-offset')


def validate(folder,manifest,text,runtime):
    """模块二：新模型独立拒绝OTH/旧场；已有runtime原始元数据不篡改。"""
    folder=Path(folder)
    assert manifest['caseKind']==KIND and manifest['selectedFluidModel']=='IDEAL_GAS'
    spec=validate_development_contract(manifest)
    assert manifest['requiredRuntimeFlavor']==runtime['flavor']==FLAVOR
    runtime_file=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)/'runtime.json'
    assert sha(runtime_file)==RUNTIME_SHA
    assert manifest['runtimeIdentity']==dict(runtimeSHA256=RUNTIME_SHA,executableSHA256=EXE_SHA)
    assert sha(Path(runtime['executable']))==EXE_SHA and runtime['asciiSignificantDigits']==17
    assert runtime['spatialWall']==SPATIAL and runtime['numericalParentFlavor']=='rah-ideal-air-spatial-wall-v1'
    assert runtime['nativeIdealGasOutputContract']==OUTPUT_CONTRACT==manifest['nativeIdealGasOutputContract']
    assert runtime['supportedModel']=='IDEAL_GAS' and runtime['thermodynamicModel']==MODEL
    assert runtime['gamma']==GAMMA and runtime['gasConstant']==R
    assert runtime['inheritedPackageThermodynamicModel']=='Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'
    assert runtime['inheritedPackageTemperatureRangeK']==[60,2000] and 'temperatureRangeK' not in runtime
    assert runtime['firstFrameOutputModelExtension']=='IDEAL_GAS' and runtime['samePhysicsAndNumerics'] is False
    assert runtime['numericalUpdateSourceUnchanged'] is True and runtime['sameNumericalUpdateSources'] is True
    assert runtime['outputScheduleChanged'] is True and runtime['nativeSharedFrameEqualityVerified'] is False
    assert runtime['sharedStateRestored'] is True and runtime['incrementalObjectsRebuilt']==1
    assert runtime['inheritedVerifiedObjects']==297
    assert manifest['runtimeUse']==RUNTIME_USE
    assert manifest['planSHA256']==sha(PLAN) and not manifest['physicalQualified']
    assert 'idealAir' not in manifest and 'sst' not in manifest
    assert not (folder/'ideal_air_contract.json').exists() and not (folder/'sst_contract.json').exists()
    c=conditions();assert manifest['publicBenchmark']==dict(conditions=c,fixture=not spec['actualValidationMesh'],experimentalQualification=False)
    assert text==configuration(spec,c),'未登记配置，禁止静默更改边界/物性/推进'
    inlet_state_check(text,c,manifest['expectedVelocity'])
    p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
    if spec['actualValidationMesh']:
        assert manifest['actualValidationMesh']==validate_actual_mesh_files(folder)
        assert manifest['developmentPrerequisites']==actual_prerequisite_summary()
        assert len(p)==190334 and len(cells)==179568 and len(meta['s'])==130
        assert set(markers)=={'wall','outer','outlet','sym_y','sym_z'} and len(markers['wall'])==3096
    else:
        ep,ec,em,eg,el,ei=fixture_geometry()
        assert np.array_equal(p,ep) and cells==ec and markers==em and meta==eg
    with np.load(str(folder/'mesh_arrays.npz')) as arrays:
        assert np.array_equal(arrays['points'],p)
        labels=arrays['nodeIndices']
        if spec['actualValidationMesh']:
            assert labels.shape==(len(p),3) and arrays['columnNodeIDs'].shape==(len(meta['offsets']),len(meta['s']),len(meta['azimuth']))
        else:assert np.array_equal(labels,el) and np.array_equal(arrays['columnNodeIDs'],ei)
        for key,kind in [('prisms',13),('hexes',12)]:
            assert np.array_equal(arrays[key],np.asarray([row[1:] for row in cells if row[0]==kind]))
    assert manifest['nodes']==len(p) and manifest['elements']==len(cells)
    wallids=np.unique([node for face in markers['wall'] for node in face[1:]])
    from nacaSpatialWallFlux import read_temperature_table,temperature_at
    table=read_temperature_table(folder/'wall_temperature.csv')
    assert np.array_equal(table,expected_table('spatial-profile',meta))
    if spec['actualValidationMesh']:assert table.shape==(130,2) and len(wallids)==3226
    tw=temperature_at(p[wallids,0],table);assert np.all((tw>=150)&(tw<=360))
    assert sha(folder/'wall_temperature.csv')==manifest['thermalProfileSHA256']
    assert manifest['thermalBoundary']==dict(mode='axialTemperature',givenTw_K=335.,spatialWall=SPATIAL)
    contract=read(folder/'ideal_gas_contract.json')
    assert manifest['idealGas']==contract and sha(folder/'ideal_gas_contract.json')==manifest['idealGasContractSHA256']
    assert contract['schema']=='RAH.NACA.NativeIdealGas.v1' and contract['initialization']==KIND
    assert contract['selectedFluidModel']=='IDEAL_GAS' and contract['thermodynamicModel']==MODEL
    assert contract['gamma']==GAMMA and contract['gasConstant']==R
    assert contract['energyReference']=='e=cv*T; h=cp*T; no-OTH-offset'
    assert contract['developmentTemperatureRangeK']==[150.,360.]
    assert not contract['turbulent'] and not contract['physicalQualified'] and contract['newQualifiedArea_m2']==0
    assert contract['reusedRuntimeMetadataIsSelectedModel'] is False
    assert contract['inheritedPackageMetadataIsSelectedModel'] is False
    assert contract['selectedRuntimeMetadataIsSelectedModel'] is True
    for name,key in [('case.cfg','configurationSHA256'),('mesh.su2','meshSHA256'),
                     ('geometry.json','geometrySHA256'),('mesh_arrays.npz','meshArraysSHA256')]:
        assert sha(folder/name)==contract[key]
    assert contract['checkpointContract']==checkpoint_plan(spec,len(p),len(wallids),len(markers['wall']))
    init=manifest['initialization']
    assert init['kind']==KIND and init['priorSuccessfulUpdates']==0 and init['initialBlendThickness_m']==.001
    assert not init['boundaryLayerSolved'] and init['reusesPreviousConservativeState'] is False
    assert init['energyReference']==contract['energyReference']
    assert sha(folder/'initial_flow.csv')==init['initialStateSHA256']
    with (folder/'initial_flow.csv').open() as stream:
        assert next(csv.reader(stream))==['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy']
    a=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1)
    assert np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(len(p))) and np.array_equal(a[:,1:4],p)
    wanted,t=initial_state(p,meta,labels,c)
    assert np.max(abs(a[:,4:]-wanted)/np.maximum(1,abs(wanted)))<2e-12
    assert init['temperatureRangeK']==[float(t.min()),float(t.max())]
    locks=dict(init['sourceSHA256'])
    assert {Path(k).resolve() for k in locks}=={p.resolve() for p in source_files(spec['variant'])}, '必需输入/源码锁不完整'
    assert locks and all(sha(Path(k))==v for k,v in locks.items())
    for name in ('case_manifest.json','ideal_gas_contract.json','geometry.json','mesh_arrays.npz',
                 'initial_flow.csv','wall_temperature.csv'):
        locks[str(folder/name)]=sha(folder/name)
    return locks


def resource_budget(request,manifest):
    """模块三：容量、线程、运行秒数与同一具名计划绑定。"""
    spec=validate_development_contract(manifest)
    assert manifest['planSHA256']==sha(PLAN)
    if spec['actualValidationMesh']:
        assert request.get('resourceProfile')==ACTUAL_RESOURCE_PROFILE
        assert 'V174_ACTUAL_NATIVE_IDEAL_GAS_REGISTERED' in PLAN.read_text(encoding='utf8')
    else:assert 'resourceProfile' not in request
    assert request['threads']==spec['threads'] and request['maximumSeconds']==spec['maximumSeconds']
    plan=manifest['idealGas']['checkpointContract']
    assert plan['maximumUpdates']==spec['updates'] and plan['every']==spec['every']
    assert plan['maximumOutputBytes']==spec['maximumOutputBytes']
    return float(spec['maximumSeconds']),spec['threads']


def load_state(path):
    """模块四：只检查保存态有限正值；150–360K开发域由热审另判，不删越域证据。"""
    from su2SSTInterface import phase_columns
    expected=['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy']
    expected+=phase_columns(3)[:-1]
    expected+=['Residual_Density','Residual_Momentum_x','Residual_Momentum_y','Residual_Momentum_z','Residual_Energy']
    with Path(path).open(encoding='utf-8') as stream:keys=next(csv.reader(stream))
    assert keys==expected,'未声明的原生定比热输出列：'+str(keys)
    a=np.loadtxt(str(path),delimiter=',',skiprows=1)
    assert a.ndim==2 and np.array_equal(a[:,0],np.arange(len(a))) and np.isfinite(a).all()
    fields={k:a[:,i] for i,k in enumerate(keys)}
    assert min(fields['Density'])>0 and min(fields['Pressure'])>0 and min(fields['Temperature'])>0
    return a,fields


def inventory(folder,receipt):
    """模块五：按实际历史和双文件原始SHA清点，不将异常历史行当完整场。"""
    from auditSU2SSTState import vtk_tokens
    from scipy.spatial import cKDTree
    manifest=read(folder/'case_manifest.json');spec=specification(manifest['specification']['variant'])
    p,_,markers=mesh(folder/'mesh.su2');ids=np.unique([i for f in markers['wall'] for i in f[1:]])
    history={};keys=[]
    if (folder/'history.csv').exists():
        with (folder/'history.csv').open() as stream:
            reader=csv.reader(stream,skipinitialspace=True);keys=[v.strip() for v in next(reader,[])];idx=keys.index('Inner_Iter')
            for row in reader:
                try:
                    v=np.asarray(row,float)
                    if len(v)==len(keys) and np.isfinite(v).all() and v[idx]==int(v[idx]):
                        step=int(v[idx]);assert 0<=step<spec['updates'] and step not in history;history[step]=v
                except (ValueError,OverflowError):continue
    names=receipt['outputSHA256']
    steps=sorted({int(m[1]) for name in names for m in [re.fullmatch(r'(?:restart_flow|surface_flow)_(\d+)\.(?:csv|vtk)',name)] if m})
    required={'Density','Momentum','Energy','Pressure','Temperature','Mach','Pressure_Coefficient','Velocity',
        'Laminar_Viscosity','Skin_Friction_Coefficient','Heat_Flux','Y_Plus','Residual_Density','Residual_Momentum','Residual_Energy'}
    complete=[];rejected=[]
    for step in steps:
        try:
            assert 0<=step<spec['updates'] and all(i in history for i in range(step+1))
            vp=folder/('restart_flow_%06d.csv'%step);sp=folder/('surface_flow_%06d.vtk'%step)
            assert all(path.is_file() and sha(path)==names[path.name] for path in (vp,sp))
            a,_=load_state(vp);assert len(a)==len(p) and np.allclose(a[:,1:4],p,rtol=2e-14,atol=1e-14)
            xyz,surface=vtk_tokens(sp);assert set(surface)==required and len(xyz)==len(ids)
            assert np.isfinite(xyz).all() and all(np.isfinite(v).all() for v in surface.values())
            distance,mapped=cKDTree(p[ids]).query(xyz);assert max(distance)<1e-6 and len(np.unique(mapped))==len(ids)
            complete.append(dict(localIteration=step,successfulUpdates=step+1,seriesSuccessfulUpdates=step+1,
                volume=vp.name,surface=sp.name,volumeSHA256=sha(vp),surfaceSHA256=sha(sp),stage=STAGE,physicalQualified=False))
        except (AssertionError,ValueError,KeyError,OSError,RuntimeError) as error:
            rejected.append(dict(localIteration=step,reason=type(error).__name__+': '+str(error)))
    last=complete[-1] if complete else None;last_history=max(history) if history else None
    normal=receipt['exitCode']==0 and not any(receipt.get(k,False) for k in
        ('timedOut','cancelled','outputBudgetExceeded','supervisorError','modelDomainExceeded','guardError'))
    expected=list(range(0,last_history+1,spec['every'])) if last_history is not None else []
    # 小例4步仍必须完整；实际候选只认可原10步后rho/E共同-8，不能信任自报early标志。
    early=False
    if spec['actualValidationMesh'] and normal and last_history is not None and 10<=last_history<spec['updates']-1:
        assert 'rms[Rho]' in keys and 'rms[RhoE]' in keys
        early=bool(all(history[last_history][keys.index(name)]<=-8 for name in ('rms[Rho]','rms[RhoE]')))
    full=last_history==spec['updates']-1
    accepted=bool(full or early)
    if normal:
        final_required=last_history if early else spec['updates']-1
        expected=sorted(set(list(range(0,final_required+1,spec['every']))+[final_required]))
    final=bool(normal and accepted and last and last['localIteration']==last_history)
    missing=sorted(set(expected)-{x['localIteration'] for x in complete})
    result=dict(stage=STAGE,complete=complete,rejected=rejected,lastComplete=last,
        lastFiniteHistoryLocalIteration=last_history,historyRows=len(history),historyContiguous=sorted(history)==list(range(len(history))),
        seriesUpdatesAtLastFiniteHistory=last_history+1 if history else 0,normalTerminal=normal,
        normalFinalPairPresent=final,expectedLocalIterations=expected,
        missingScheduledLocalIterations=missing,normalScheduleComplete=bool(final
            and not missing and not rejected and sorted(history)==list(range(len(history)))),
        abnormalInstantStateAvailable=False,physicsQualified=False)
    if spec['actualValidationMesh']:
        result.update(earlyConvergenceAccepted=early,completedMaximumUpdates=full,terminationRequirementPassed=accepted)
    return result
