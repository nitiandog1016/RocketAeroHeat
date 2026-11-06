"""模块用途：NACA三维层流小例的专用入口和检查点，不放宽旧SST合同。"""
from pathlib import Path
import csv
import re

import numpy as np

from idealAirProperties import IdealAir
from prepareSU2NACALaminar import KIND, STAGE, config_text, capacity, small_geometry, initial_state
from su2SSTInterface import phase_columns
from su2ThermalRuntime import ROOT, read, sha


def mesh(path):
    """模块一：读取所有具名热边界和对称面，保留共同节点的多个归属。"""
    with Path(path).open() as f:
        assert f.readline().strip()=='NDIME= 3'
        nc=int(f.readline().split('=')[1])
        cells=[list(map(int,f.readline().split())) for _ in range(nc)]
        n=int(f.readline().split('=')[1])
        rows=np.array([list(map(float,f.readline().split())) for _ in range(n)])
        assert np.array_equal(rows[:,3],np.arange(n))
        markers={}
        for _ in range(int(f.readline().split('=')[1])):
            name=f.readline().split('=')[1].strip();assert name not in markers
            size=int(f.readline().split('=')[1])
            markers[name]=[list(map(int,f.readline().split())) for _ in range(size)]
        assert not f.read().strip()
    return rows[:,:3],cells,markers


def validate(folder,manifest,text,runtime):
    """模块二：完整配置/几何/初场逐项绑定；三维/层流/壁温被改即拒绝。"""
    from su2IdealAirInterface import AXISYMMETRIC_FLAVOR,OBSERVED_FLAVOR,OBSERVATION,MODEL,R
    from naca3344Model import source_data
    assert manifest['caseKind']==KIND and runtime['flavor'] in (AXISYMMETRIC_FLAVOR,OBSERVED_FLAVOR)
    assert manifest['requiredRuntimeFlavor']==runtime['flavor']
    precise=runtime['flavor']==OBSERVED_FLAVOR
    if precise:
        assert manifest['runtimeObservation']==OBSERVATION
        assert all(runtime[k]==v for k,v in OBSERVATION.items())
    else:assert 'runtimeObservation' not in manifest
    assert manifest['runtimeUse']=='existing-general-3d-NS-with-axisymmetric-source-disabled'
    contract=read(folder/'ideal_air_contract.json')
    assert contract==manifest['idealAir'] and sha(folder/'ideal_air_contract.json')==manifest['idealAirContractSHA256']
    assert contract['initialization']==KIND and contract['schemaVersion']==5 and not contract['turbulent']
    assert contract['thermodynamicModel']==runtime['thermodynamicModel']==MODEL
    assert contract['gasConstant']==R and contract['temperatureRangeK']==[60,2000]
    assert contract['energyEntropyReference']=='OTH'
    assert not contract['physicsQualified'] and contract['newQualifiedArea']==0
    assert not (folder/'sst_contract.json').exists() and 'sst' not in manifest
    assert not manifest['physicalQualified'] and manifest['publicBenchmark']['fixture']
    assert manifest['publicBenchmark']['comparisonReady'] is False
    c=source_data()[0]['freeStream']
    assert c==manifest['publicBenchmark']['conditions']
    p,cells,markers=mesh(folder/'mesh.su2')
    expected,expected_cells,expected_markers,meta=small_geometry()
    assert np.array_equal(p,expected) and cells==expected_cells and markers==expected_markers
    assert read(folder/'geometry.json')==meta
    assert text==config_text(c,meta,3,1,precise), '专用小例配置变化，不能套用此合同'
    for name,key in [('case.cfg','configurationSHA256'),('mesh.su2','meshSHA256'),('geometry.json','geometrySHA256')]:
        assert sha(folder/name)==contract[key]
    assert manifest['nodes']==len(p) and manifest['elements']==len(cells)
    wall_ids=np.unique([i for k,v in markers.items() if k.startswith('wall_') for row in v for i in row[1:]])
    wall_faces=sum(len(v) for k,v in markers.items() if k.startswith('wall_'))
    assert contract['checkpointContract']==capacity(len(p),len(wall_ids),wall_faces,3,1,32*1024**2)
    init=manifest['initialization']
    assert init['kind']==KIND and init['priorSuccessfulUpdates']==0 and not init['boundaryLayerSolved']
    assert sha(folder/'initial_flow.csv')==init['initialStateSHA256']
    with (folder/'initial_flow.csv').open() as f:
        assert next(csv.reader(f))==['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy']
    a=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1)
    assert np.array_equal(a[:,0],np.arange(len(p))) and np.array_equal(a[:,1:4],p)
    core=IdealAir();state,thermo=initial_state(p,meta,core,c)
    assert np.max(abs(a[:,4:]-state)/np.maximum(1,abs(state)))<2e-12
    assert np.isfinite(a).all() and min(thermo['T'])>=60 and max(thermo['T'])<=2000
    locks=dict(init['sourceSHA256']);locks.update(core.locks)
    assert all(sha(Path(p))==h for p,h in locks.items())
    for file in (folder/'ideal_air_contract.json',folder/'case_manifest.json',folder/'geometry.json',
                 folder/'initial_flow.csv',Path(__file__).resolve()):
        locks[str(file)]=sha(file)
    return locks


def load_state(path):
    """模块三：没有SST也没有前后k阶段；只接受三维同次流动保存的具名列。"""
    expected=['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Energy']
    expected+=phase_columns(3)[:-1]
    expected+=['Residual_Density','Residual_Momentum_x','Residual_Momentum_y','Residual_Momentum_z','Residual_Energy']
    with path.open(encoding='utf8') as f:keys=next(csv.reader(f))
    assert keys==expected, '未声明的层流输出列：'+str(keys)
    a=np.loadtxt(str(path),delimiter=',',skiprows=1)
    assert a.ndim==2 and np.array_equal(a[:,0],np.arange(len(a))) and np.isfinite(a).all()
    fields={key:a[:,i] for i,key in enumerate(keys)}
    assert min(fields['Density'])>0 and min(fields['Pressure'])>0
    assert min(fields['Temperature'])>=60 and max(fields['Temperature'])<=2000
    return a,fields


def inventory(folder,receipt):
    """模块四：按原编号核对计划全部帧，公开基准无物理资格仍明确列出。"""
    from auditSU2SSTState import vtk_tokens
    manifest=read(folder/'case_manifest.json');plan=manifest['idealAir']['checkpointContract']
    p,_,markers=mesh(folder/'mesh.su2')
    ids=np.unique([i for k,v in markers.items() if k.startswith('wall_') for row in v for i in row[1:]])
    history={}
    if (folder/'history.csv').exists():
        with (folder/'history.csv').open() as f:
            reader=csv.reader(f,skipinitialspace=True);keys=[k.strip() for k in next(reader)]
            for row in reader:
                try:
                    values=np.array(row,float);step=values[keys.index('Inner_Iter')]
                    if len(values)==len(keys) and np.isfinite(values).all() and step==int(step):
                        assert int(step) not in history
                        history[int(step)]=values
                except ValueError:continue
    names=receipt['outputSHA256'];steps=sorted({int(m[1]) for name in names
        for m in [re.fullmatch(r'(?:restart_flow|surface_flow)_(\d+)\.(?:csv|vtk)',name)] if m})
    complete,rejected=[],[]
    expected_fields={'Density','Momentum','Energy','Pressure','Temperature','Mach','Pressure_Coefficient',
        'Velocity','Laminar_Viscosity','Skin_Friction_Coefficient','Heat_Flux','Y_Plus',
        'Residual_Density','Residual_Momentum','Residual_Energy'}
    for step in steps:
        try:
            assert 0<=step<plan['maximumUpdates'] and all(i in history for i in range(step+1))
            vp=folder/('restart_flow_%06d.csv'%step);sp=folder/('surface_flow_%06d.vtk'%step)
            assert all(path.is_file() and sha(path)==names[path.name] for path in (vp,sp))
            a,fields=load_state(vp)
            assert len(a)==len(p) and np.allclose(a[:,1:4],p,rtol=2e-14,atol=1e-14)
            xyz,surface=vtk_tokens(sp)
            assert set(surface)==expected_fields and len(xyz)==len(ids)
            assert np.isfinite(xyz).all() and all(np.isfinite(v).all() for v in surface.values())
            complete.append(dict(localIteration=step,successfulUpdates=step+1,seriesSuccessfulUpdates=step+1,
                volume=vp.name,surface=sp.name,volumeSHA256=sha(vp),surfaceSHA256=sha(sp),
                stage=STAGE,physicalQualified=False))
        except (AssertionError,ValueError,KeyError,OSError,RuntimeError) as e:
            rejected.append(dict(localIteration=step,reason=type(e).__name__+': '+str(e)))
    last=complete[-1] if complete else None
    last_history=max(history) if history else None
    normal=receipt['exitCode']==0 and not any(receipt.get(k,False) for k in
        ('timedOut','cancelled','outputBudgetExceeded','supervisorError'))
    expected=list(range(last_history+1)) if last_history is not None else []
    missing=sorted(set(expected)-{r['localIteration'] for r in complete})
    return dict(stage=STAGE,complete=complete,rejected=rejected,lastComplete=last,
        lastFiniteHistoryLocalIteration=last_history,historyRows=len(history),
        historyContiguous=sorted(history)==list(range(len(history))),
        seriesUpdatesAtLastFiniteHistory=last_history+1 if history else 0,
        normalTerminal=normal,normalFinalPairPresent=bool(normal and last and last['localIteration']==last_history),
        expectedLocalIterations=expected,missingScheduledLocalIterations=missing,
        normalScheduleComplete=bool(normal and last and not missing and not rejected),
        abnormalInstantStateAvailable=False,physicsQualified=False)
