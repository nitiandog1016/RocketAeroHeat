"""模块用途：单wall空间Tw的严格新合同和稀疏帧，不改变旧小例或SST入口。"""
import csv
import re
from pathlib import Path
import numpy as np

from prepareNACASpatialCase import KIND,FLAVOR,SPATIAL,BASE,PLAN,specification,fixture_geometry,spatial_config,expected_table
from prepareSU2NACALaminar import initial_state,capacity,STAGE
from su2NACALaminarInterface import mesh,load_state,inlet_state_check
from idealAirProperties import IdealAir
from naca3344Model import source_data
from su2ThermalRuntime import ROOT,read,sha


def read_profile(path,wall_x):
    """模块一：显式SI/有效温区/递增/全壁覆盖预检；不外推或降级到常温。"""
    with Path(path).open(encoding='ascii') as f:
        assert f.readline().strip()==SPATIAL['schema'],'空间Tw表头或单位不符'
        rows=list(csv.reader(f))
    assert 2<=len(rows)<=10000 and all(len(r)==2 for r in rows)
    a=np.array(rows,float)
    assert np.isfinite(a).all() and np.all(np.diff(a[:,0])>0)
    assert np.all((a[:,1]>=60)&(a[:,1]<=2000)),'空间Tw超出已支持物性温区'
    x=np.asarray(wall_x,float);assert np.isfinite(x).all() and x.min()>=a[0,0] and x.max()<=a[-1,0]
    return a,np.interp(x,a[:,0],a[:,1])


def validate(folder,manifest,text,runtime):
    """模块二：包、配置、网格、表、初场及固定计划全部逐项绑定后才允许启动。"""
    from su2IdealAirInterface import OBSERVATION,OBSERVED_FLAVOR,MODEL,R
    assert manifest['caseKind']==KIND
    spec=specification(manifest['specification']['variant']);assert spec==manifest['specification']
    assert manifest['requiredRuntimeFlavor']==runtime['flavor']==spec['requiredRuntimeFlavor']
    assert all(runtime[k]==v for k,v in OBSERVATION.items())
    if runtime['flavor']==FLAVOR:
        assert runtime['spatialWall']==SPATIAL and runtime['numericalParentFlavor']==OBSERVED_FLAVOR
        assert runtime['samePhysicsAndNumerics'] is False and runtime['samePhysicsAndNumericsWithoutSpatialTable'] is True
    else:assert runtime['flavor']==OBSERVED_FLAVOR and spec['variant']=='observed-scalar'
    assert runtime['thermodynamicModel']==MODEL
    assert manifest['planSHA256']==sha(PLAN) and not manifest['physicalQualified']
    c=source_data()[0]['freeStream'];assert manifest['publicBenchmark']['conditions']==c
    assert manifest['publicBenchmark']['fixture']==(not spec['actualValidationMesh'])
    assert not manifest['publicBenchmark']['experimentalQualification']
    p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json')
    assert set(markers)=={'wall','outer','outlet','sym_y','sym_z'}
    with np.load(str(folder/'mesh_arrays.npz')) as a:
        assert np.array_equal(a['points'],p);labels=a['nodeIndices'];ix=a['columnNodeIDs']
        for key,kind in [('prisms',13),('hexes',12)]:
            assert np.array_equal(a[key],np.array([r[1:] for r in cells if r[0]==kind]))
    if spec['actualValidationMesh']:
        proof=read(ROOT/'docs/validation/v170/evidence_manifest.json')
        for row in proof['localOnly']:
            assert sha(folder/Path(row['path']).name)==row['sha256']
        for name in ('geometry.json','wall_temperature.csv'):
            row=next(r for r in proof['files'] if r['source']=='outputs/v170_validation/mesh02/'+name)
            assert sha(folder/name)==row['sha256'] and sha(BASE/name)==row['sha256']
        assert meta==read(BASE/'geometry.json') and len(p)==190334 and len(cells)==179568
    else:
        expected,ec,em,eg,el,ei=fixture_geometry()
        assert np.array_equal(p,expected) and cells==ec and markers==em and meta==eg
        assert np.array_equal(labels,el) and np.array_equal(ix,ei)
    assert manifest['nodes']==len(p) and manifest['elements']==len(cells)
    assert text==spatial_config(spec['variant'],c,meta),'未登记的配置改变'
    inlet_state_check(text,c,manifest['expectedVelocity'])
    wallids=np.unique([v for face in markers['wall'] for v in face[1:]])
    thermal=manifest['thermalBoundary']
    assert thermal==dict(mode=spec['thermalMode'],givenTw_K=335.,spatialWall=SPATIAL if spec['thermalMode']=='axialTemperature' else None)
    profile=folder/'wall_temperature.csv'
    if spec['thermalMode']=='axialTemperature':
        a,_=read_profile(profile,p[wallids,0])
        assert sha(profile)==manifest['thermalProfileSHA256'] and np.array_equal(a,expected_table(spec['variant'],meta))
    else:assert not profile.exists() and manifest['thermalProfileSHA256'] is None
    contract=read(folder/'ideal_air_contract.json')
    assert contract==manifest['idealAir'] and sha(folder/'ideal_air_contract.json')==manifest['idealAirContractSHA256']
    assert contract['schemaVersion']==6 and contract['initialization']==KIND and contract['thermodynamicModel']==MODEL
    assert contract['gasConstant']==R and contract['temperatureRangeK']==[60,2000] and contract['energyEntropyReference']=='OTH'
    assert not contract['turbulent'] and not contract['physicsQualified'] and contract['newQualifiedArea']==0
    assert 'sst' not in manifest and not (folder/'sst_contract.json').exists()
    for name,key in [('case.cfg','configurationSHA256'),('mesh.su2','meshSHA256'),('geometry.json','geometrySHA256'),('mesh_arrays.npz','meshArraysSHA256')]:
        assert sha(folder/name)==contract[key]
    assert contract['checkpointContract']==capacity(len(p),len(wallids),len(markers['wall']),spec['updates'],spec['every'],spec['maximumOutputBytes'])
    init=manifest['initialization'];assert init['kind']==KIND and init['priorSuccessfulUpdates']==0
    assert init['initialBlendThickness_m']==.001 and not init['boundaryLayerSolved']
    assert sha(folder/'initial_flow.csv')==init['initialStateSHA256']
    with (folder/'initial_flow.csv').open() as f:
        assert next(csv.reader(f))==['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy']
    a=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1)
    assert np.isfinite(a).all() and np.array_equal(a[:,0],np.arange(len(p))) and np.array_equal(a[:,1:4],p)
    core=IdealAir();wanted,thermo=initial_state(p,dict(meta,nodeIndices=labels.tolist()),core,c)
    assert np.max(abs(a[:,4:]-wanted)/np.maximum(1,abs(wanted)))<2e-12
    assert min(thermo['T'])>=60 and max(thermo['T'])<=2000
    locks=dict(init['sourceSHA256']);locks.update(core.locks)
    assert all(sha(Path(p))==h for p,h in locks.items())
    for name in ('case_manifest.json','ideal_air_contract.json','geometry.json','mesh_arrays.npz','initial_flow.csv'):
        locks[str(folder/name)]=sha(folder/name)
    if profile.exists():locks[str(profile)]=sha(profile)
    return locks


def resource_budget(request,manifest):
    """模块三：只有具名2000更新候选可使用1800秒，不扩大任何旧入口限额。"""
    spec=specification(manifest['specification']['variant']);assert spec==manifest['specification']
    assert manifest['caseKind']==KIND and manifest['planSHA256']==sha(PLAN)
    assert request['threads']==spec['threads'] and request['maximumSeconds']==spec['maximumSeconds']
    if spec['actualValidationMesh']:
        assert request.get('resourceProfile')=='v171-naca-spatial-validation'
        assert manifest['requiredRuntimeFlavor']==FLAVOR and manifest['nodes']==190334
    else:assert 'resourceProfile' not in request
    plan=manifest['idealAir']['checkpointContract']
    assert plan['maximumUpdates']==spec['updates'] and plan['every']==spec['every'] and plan['maximumOutputBytes']==spec['maximumOutputBytes']
    return float(spec['maximumSeconds']),spec['threads']


def inventory(folder,receipt):
    """模块四：稀疏首/周期/正常末帧；异常最后历史不是自动应有的瞬间场。"""
    from auditSU2SSTState import vtk_tokens
    from scipy.spatial import cKDTree
    manifest=read(folder/'case_manifest.json');spec=specification(manifest['specification']['variant'])
    p,_,markers=mesh(folder/'mesh.su2');ids=np.unique([i for f in markers['wall'] for i in f[1:]])
    history={}
    if (folder/'history.csv').exists():
        with (folder/'history.csv').open() as f:
            reader=csv.reader(f,skipinitialspace=True);keys=[k.strip() for k in next(reader,[])];idx=keys.index('Inner_Iter')
            for row in reader:
                try:
                    v=np.array(row,float)
                    if len(v)==len(keys) and np.isfinite(v).all() and v[idx]==int(v[idx]):
                        step=int(v[idx]);assert 0<=step<spec['updates'] and step not in history;history[step]=v
                except (ValueError,OverflowError):continue
    names=receipt['outputSHA256'];steps=sorted({int(m[1]) for name in names for m in
        [re.fullmatch(r'(?:restart_flow|surface_flow)_(\d+)\.(?:csv|vtk)',name)] if m})
    complete=[];rejected=[]
    required={'Density','Momentum','Energy','Pressure','Temperature','Mach','Pressure_Coefficient','Velocity',
        'Laminar_Viscosity','Skin_Friction_Coefficient','Heat_Flux','Y_Plus','Residual_Density','Residual_Momentum','Residual_Energy'}
    for step in steps:
        try:
            assert 0<=step<spec['updates'] and all(i in history for i in range(step+1))
            vp=folder/('restart_flow_%06d.csv'%step);sp=folder/('surface_flow_%06d.vtk'%step)
            assert all(path.is_file() and sha(path)==names[path.name] for path in (vp,sp))
            a,fields=load_state(vp);assert len(a)==len(p) and np.allclose(a[:,1:4],p,rtol=2e-14,atol=1e-14)
            xyz,surface=vtk_tokens(sp);assert set(surface)==required and len(xyz)==len(ids)
            assert np.isfinite(xyz).all() and all(np.isfinite(v).all() for v in surface.values())
            distance,mapped=cKDTree(p[ids]).query(xyz);assert max(distance)<1e-6 and len(np.unique(mapped))==len(ids)
            complete.append(dict(localIteration=step,successfulUpdates=step+1,seriesSuccessfulUpdates=step+1,
                volume=vp.name,surface=sp.name,volumeSHA256=sha(vp),surfaceSHA256=sha(sp),stage=STAGE,physicalQualified=False))
        except (AssertionError,ValueError,KeyError,OSError,RuntimeError) as e:
            rejected.append(dict(localIteration=step,reason=type(e).__name__+': '+str(e)))
    last=complete[-1] if complete else None;last_history=max(history) if history else None
    normal=receipt['exitCode']==0 and not any(receipt.get(k,False) for k in ('timedOut','cancelled','outputBudgetExceeded','supervisorError'))
    expected=list(range(0,last_history+1,spec['every'])) if last_history is not None else []
    if normal and last_history is not None:expected=sorted(set(expected)|{last_history})
    missing=sorted(set(expected)-{r['localIteration'] for r in complete})
    return dict(stage=STAGE,complete=complete,rejected=rejected,lastComplete=last,lastFiniteHistoryLocalIteration=last_history,
        historyRows=len(history),historyContiguous=sorted(history)==list(range(len(history))),
        seriesUpdatesAtLastFiniteHistory=last_history+1 if history else 0,normalTerminal=normal,
        normalFinalPairPresent=bool(normal and last and last['localIteration']==last_history),
        expectedLocalIterations=expected,missingScheduledLocalIterations=missing,
        normalScheduleComplete=bool(normal and last and last['localIteration']==last_history and not missing and not rejected
            and sorted(history)==list(range(len(history)))),
        abnormalInstantStateAvailable=False,physicsQualified=False)
