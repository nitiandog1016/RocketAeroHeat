"""模块职责：一次真实场的原面、热诊断、能量残差及停止证据审核，不推进流场。"""
from pathlib import Path
import csv,json,re,shutil,sys,time
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from runSU2RocketPhaseV222 import ROOT,BASE,CASE,DOC,ENGINE,MESH,sha,read,save,supervise
from su2PhaseAirContract import bind
from su2PhaseAirFields import read_fields,inverse_energy,RESIDUAL
from phaseAwareViscousAir import PhaseAwareViscousAir

def error(a,b):return float(np.max(np.abs(np.asarray(a)-np.asarray(b))/np.maximum(1.,np.abs(b))))
def interval(x):return [float(np.min(x)),float(np.max(x))]

def geometry():
    """模块一：复用旧壁邻算法，只建立读取缓存；原面/体网格均不修改。"""
    from su2SSTInterface import wall_geometry,mesh
    out=BASE/'geometry01';out.mkdir(exist_ok=False)
    g=wall_geometry(MESH);points,cells,wall=mesh(MESH)
    source=ROOT/'outputs/v144_validation/actual/hybrid_volume.mat'
    original=loadmat(str(source),variable_names=['sourceVertices','sourceFaces'])
    v=original['sourceVertices'];f=original['sourceFaces'].astype(np.int64)-1
    distance,ids=cKDTree(points).query(v)
    assert np.all(distance==0) and len(set(ids.tolist()))==len(v)
    native=np.array([row[1:] for row in wall]);tri=ids[f]
    canonical=lambda a:sorted(map(tuple,np.sort(a,axis=1).tolist()))
    assert canonical(native)==canonical(tri) and len(tri)==19752
    area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)/2
    assert abs(float(sum(area))-1.9646292565715553)<1e-12
    np.savez_compressed(str(out/'cache.npz'),**g,triangles=tri,originalVertices=v,originalFaces=f,faceArea=area)
    save(out/'report.json',dict(passed=True,meshSHA256=sha(MESH),sourceSHA256=sha(source),cacheSHA256=sha(out/'cache.npz'),
        originalWallFaces=len(tri),wallNodes=len(ids),areaM2=float(sum(area)),algorithmSHA256=sha(ROOT/'scripts/su2SSTInterface.py')))

def software():
    """模块二：对既有24点保存文件验证新批量读取；无新原生查询或物性控制。"""
    import copy
    from su2SSTInterface import mesh
    from su2PhaseAirContract import read_restart
    out=BASE/'software01';out.mkdir(exist_ok=False)
    old=ROOT/'outputs/v221_validation/control01';p=old/'restart_flow.csv';metadata=read(old/'restart_flow.rah.json')
    expected=bind(p,old/'mesh.su2',old/'case.cfg',ENGINE,ROOT/'outputs/v221_validation/build/compatibility.exe','manufactured-control')
    points,_,_=mesh(old/'mesh.su2');batch=read_fields(p,metadata,expected,points)
    scalar=read_restart(p,metadata,expected,old/'mesh.su2');order=np.argsort(scalar['ids'])
    e=error(batch['current']['temperatureK'],scalar['fluidTemperatureK'][order]);assert e<1e-10
    rejected=[]
    for key in ('model','energy','meshSHA256','sourceConfigSHA256','engineSHA256','producerSHA256','payloadSHA256','role','sourceThermoSHA256'):
        bad=copy.deepcopy(metadata);bad[key]=None
        try:read_fields(p,bad,expected,points)
        except ValueError:rejected.append(key)
    for name,coords in [('changed-coordinate',points+1e-3),('missing-node',points[:-1])]:
        try:read_fields(p,metadata,expected,coords)
        except ValueError:rejected.append(name)
    assert len(rejected)==11
    save(out/'report.json',dict(passed=True,batchVsExistingScalarTemperatureScaledError=e,rejected=rejected,nativeCalls=0))

def audit():
    """模块三：只审核实际终止后完整帧，缺保存、未收敛和物理缺项均单列。"""
    from auditSU2SSTState import vtk_tokens
    out=BASE/'audit01';out.mkdir(exist_ok=False);start=time.monotonic()
    process=read(BASE/'process_engine/report.json');assert process['inputsStable'] and not process['liveGroupAfterCleanup']
    gmeta=read(BASE/'geometry01/report.json');assert gmeta['meshSHA256']==sha(MESH) and gmeta['cacheSHA256']==sha(BASE/'geometry01/cache.npz')
    g=np.load(str(BASE/'geometry01/cache.npz'));points=g['points'];wall=g['wallIDs'];tri=g['triangles'];area=g['faceArea']
    with (CASE/'history.csv').open() as stream:history=list(csv.DictReader(stream,skipinitialspace=True))
    history=[{k.strip():float(v) for k,v in row.items()} for row in history]
    steps=[int(r['Inner_Iter']) for r in history];assert steps==list(range(len(steps)))
    records=[];incomplete=[];air=PhaseAwareViscousAir();previous=None;native_files={}
    inf=air.evaluate(1197.0279007383137,226.50908361133003);speed=5*float(inf['soundSpeedMps']);rhoinf=float(inf['densityKgm3'])
    # 长度1m与速度a_inf定义研究归一化；不是修改原生-8历史门槛。
    sound=float(inf['soundSpeedMps']);scales=np.array([rhoinf*sound]+[rhoinf*sound**2]*3+[rhoinf*sound**3])
    for path in sorted(CASE.glob('restart_flow_*.csv')):
        match=re.fullmatch(r'restart_flow_(\d+)\.csv',path.name)
        if not match:continue
        step=int(match[1]);surface=CASE/('surface_flow_%06d.vtk'%step)
        if step not in steps or not surface.exists():incomplete.append(path.name);continue
        metadata=bind(path,CASE/'mesh.su2',CASE/'case.cfg',ENGINE,ENGINE,'research-state')
        frame=read_fields(path,metadata,metadata,points);v=frame['values'];native_files[path.name]=sha(path);native_files[surface.name]=sha(surface)
        save(out/(path.stem+'.rah.json'),dict(metadata=metadata,iteration=step,successfulUpdates=step+1,sourceReceiptSHA256=sha(BASE/'process_engine/report.json')))
        xyz,surf=vtk_tokens(surface);dist,sids=cKDTree(points).query(xyz)
        assert max(dist)<1e-6 and set(sids)==set(wall)
        p=v['Pressure'];t=v['Temperature'];mu=air.interpolated_transport(p,t)['dynamicViscosityPas']
        k=air.interpolated_transport(p[wall],t[wall])['conductivityWmK']
        q=k*(t[g['normalNeighbor']]-300)/g['thermalDistance']
        checks=dict(viscosity=error(v['Laminar_Viscosity'],mu),heatGradient=error(v['Heat_Flux'][wall],q),
            primitiveCurrentPressure=error(p,frame['current']['pressurePa']),primitiveCurrentTemperature=error(t,frame['current']['temperatureK']),
            primitiveCurrentVelocity=error(np.column_stack([v['Velocity_'+a] for a in ('x','y','z')]),frame['currentVelocity']))
        for name in ('Pressure','Temperature','Heat_Flux','Y_Plus'):checks['surface'+name]=error(surf[name][:,0],v[name][sids])
        cf=np.column_stack([v['Skin_Friction_Coefficient_'+a][wall] for a in ('x','y','z')])
        tau=np.linalg.norm(cf,axis=1)*.5*rhoinf*speed**2
        rho_stage=frame['primitive']['densityKgm3'][wall]
        yplus=g['nearestWallDistance']*np.sqrt(rho_stage*tau)/mu[wall]
        checks['yPlus']=error(v['Y_Plus'][wall],yplus)
        power=float(g['nodeArea']@q);checks['historyPower']=error(power,history[step]['HF'])
        checks['wallEddyFraction']=float(np.max(abs(v['Eddy_Viscosity'][wall])/mu[wall]))
        raw=np.column_stack([frame['nativeResidual'],v['Residual_Nu_Tilde']])
        rms=np.sqrt(np.mean(raw**2,axis=0));logs=np.log10(np.maximum(rms,1e-16))
        hist=np.array([history[step]['rms['+n+']'] for n in ('Rho','RhoU','RhoV','RhoW','RhoE','nu')])
        checks['historyResidualLog']=float(np.max(abs(logs-hist)))
        fp=p[tri].mean(axis=1);fq=v['Heat_Flux'][tri].mean(axis=1);ft=t[tri].mean(axis=1)
        physical=frame['physicalResidual'];normalized=np.sqrt(np.mean((physical/scales)**2,axis=0))
        record=dict(iteration=step,successfulUpdates=step+1,checks=checks,currentTemperatureRangeK=interval(frame['current']['temperatureK']),
            primitiveTemperatureRangeK=interval(t),wallFluidTemperatureRangeK=interval(t[wall]),wallPressureRangePa=interval(p[wall]),
            wallHeatRangeWm2=interval(q),givenWallTemperatureK=300,solidTemperature=None,wallYPlusRange=interval(v['Y_Plus'][wall]),
            originalFacePowerW=float(area@fq),nativeDualVectorAreaPowerW=power,originalFaceAreaM2=float(sum(area)),
            physicalRHS_RMS_SI=np.sqrt(np.mean(physical**2,axis=0)).tolist(),physicalRHS_normalizedRMS=normalized.tolist(),
            signedPhysicalRHSSumSI=physical.sum(axis=0).tolist(),nativeHistoryLogRMS=hist.tolist(),nativeMinus8Passed=bool(np.all(hist<=-8)),
            stateJointDomainAccepted=True,postUpdateFullPhysicalBalanceVerified=False,physicalQualification=False)
        required={k:v for k,v in checks.items() if not k.startswith('primitiveCurrent')}
        record['diagnosticInterfacePassed']=all(v< (2e-6 if k.startswith('surface') or k=='historyPower' else 2e-7 if k=='historyResidualLog' else 2e-10) for k,v in required.items())
        if previous is not None:
            record['previousSaveChange']=dict(pressureAreaL1Relative=float(area@abs(fp-previous[0])/max(area@abs(previous[0]),1)),
                heatAreaL1Relative=float(area@abs(fq-previous[1])/max(area@abs(previous[1]),1)))
        previous=(fp,fq)
        with (out/('original_faces_%06d.csv'%step)).open('w',newline='',encoding='utf8') as stream:
            writer=csv.writer(stream);writer.writerow(['originalFaceIDOneBased','areaM2','pressurePa','heatFluxIntoWallWm2','wallFluidTemperatureK','givenWallTemperatureK','currentConservativeWallTemperatureK'])
            writer.writerows(zip(np.arange(1,len(tri)+1),area,fp,fq,ft,np.full(len(tri),300),frame['current']['temperatureK'][tri].mean(axis=1)))
        records.append(record);save(out/'progress.json',records)
    expected=read(CASE/'case_manifest.json')['expectedSavedIterations'];actual=[r['iteration'] for r in records]
    save(out/'report.json',dict(passed=bool(records) and all(r['diagnosticInterfacePassed'] for r in records),
        engineNormal=process['passed'],historyUpdates=len(steps),expectedSaves=expected,actualCompleteSaves=actual,
        missingSaves=[n for n in expected if n not in actual],incompleteFiles=incomplete,planAcceptance=process['passed'] and actual==expected,
        records=records,outputSHA256=native_files,physicalQualification=False,newQualifiedAreaM2=0,
        fullBalanceVerified=False,spaceAndOuterDomainErrorsVerified=False,seconds=time.monotonic()-start,
        residualStage='native implicit RHS, negative pre-update spatial residual with wall equation replacement; not post-update divergence',
        energyResidualConversion='per-node RE=REstar-300000*Rmass before norm',normalization='rho_inf, a_inf, L=1m; native -8 kept separately')))
    print(json.dumps(dict(completeSaves=actual,historyUpdates=len(steps),physicalQualification=False)))

if __name__=='__main__':
    stage=sys.argv[1]
    if stage.startswith('task_'):{'geometry':geometry,'software':software,'audit':audit}[stage[5:]]()
    else:
        locks=[Path(__file__).resolve(),ROOT/'scripts/su2PhaseAirFields.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            ROOT/'scripts/su2PhaseAirContract.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
            ROOT/'scripts/su2SSTInterface.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json',MESH]
        if stage=='audit':locks += [BASE/'process_engine/report.json',CASE/'case.cfg',CASE/'history.csv',BASE/'geometry01/cache.npz']
        raise SystemExit(0 if supervise(stage,[sys.executable,str(Path(__file__).resolve()),'task_'+stage],120,256,locks,BASE,'audit' if stage=='audit' else 'offline') else 1)
