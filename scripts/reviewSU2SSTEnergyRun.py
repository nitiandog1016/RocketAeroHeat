"""模块用途。独立审核SST热内能约束、原面场与事前收益；缺失末帧保持待验。"""
import csv
import json
from pathlib import Path
import time
import numpy as np
from scipy.io import loadmat
from auditSU2IdealAirState import audit,load_state
from idealAirProperties import IdealAir
from prepareSU2WallDistanceState import boundary_faces
from reviewSU2SSTLocalDiffusion import neighbors_of
from reviewSU2AdmissibleRun import history
from reviewSU2ConservativeStepRun import volumes
from su2SSTInterface import mesh,wall_geometry
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save,run,write_rows,read_rows

BASE=ROOT/'outputs/v167_validation';CASE=BASE/'real01'
RESIDUALS=('rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]','rms[k]','rms[w]')


def trace_review(out,receipt,core):
    """模块一。按原生裁限公式重放所有实际事件，不把RCM编号当全局位置。"""
    path=CASE/'sst_energy_updates.csv'
    events=list(csv.DictReader(path.open())) if path.exists() else []
    assert len(events)<=65536
    if events:assert sha(path)==receipt['outputSHA256'][path.name]
    assert [int(x['event']) for x in events]==list(range(1,len(events)+1))
    accepted=[x for x in events if x['status']=='limited'];failed=[x for x in events if x['status']!='limited']
    bounds=core.state('PT',101325,[60,2000])['e'];error=0.;energy_error=0.;records=[];oracle_rows=[];energies=[]
    for event in events:
        d={k:float(v) for k,v in event.items() if k!='status'}
        rho=d['rho'];u=np.array([d['rhou'],d['rhov'],d['rhow']])/rho
        e0=d['rhoE']/rho-np.dot(u,u)/2
        zero_k=np.clip((d['old_k']*rho)/rho,d['lower_k'],d['upper_k'])
        old_energy=e0-zero_k
        if event['status']!='limited':
            records.append(dict(event=int(d['event']),localRCMPoint=int(d['local_point']),
                zeroStepEnergyJkg=float(old_energy),minimumEnergyJkg=float(bounds[0]),
                maximumEnergyJkg=float(bounds[1]),belowMinimumJkg=float(max(0,bounds[0]-old_energy)),
                flowConserved=[d[k] for k in ('rho','rhou','rhov','rhow','rhoE')],
                oldK=d['old_k'],actualSSTIncrement=d['increment_rhok'],
                globalPointKnown=False,status=event['status']))
            continue
        assert 0<=d['alpha']<=d['native_alpha']<=1 and d['refinements']<=64
        assert d['alpha']==0 or d['alpha']>=1e-10
        for name,inc in (('k','rhok'),('omega','rhoomega')):
            expected=np.clip((d['old_'+name]*rho+d['alpha']*d['increment_'+inc])/rho,d['lower_'+name],d['upper_'+name])
            error=max(error,abs(expected-d[name])/max(1,abs(expected)))
        energy=e0-d['k'];energy_error=max(energy_error,abs(energy-d['e']))
        assert bounds[0]<=energy<=bounds[1]
        temperature=float(core.state('rhoe',rho,energy)['T'][0])
        oracle_rows.append(('rhoT',rho,temperature));energies.append(energy)
    assert error<2e-13 and energy_error<1e-8
    oracle_error=None
    if oracle_rows:
        helper=ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
        assert sha(helper)==read(helper.parent/'runtime.json')['sha256']
        write_rows(out/'sst_oracle_input.csv',oracle_rows)
        run(out,'sst_oracle',[str(helper),'--oracle',str(out/'sst_oracle_input.csv'),str(out/'sst_oracle.csv')],30)
        oracle=read_rows(out/'sst_oracle.csv');assert np.all(oracle['accepted']==1)
        oracle_error=float(max(abs(oracle['e']-energies)/np.maximum(1,abs(np.array(energies)))))
        assert oracle_error<2e-11
    flow_path=CASE/'admissibility.csv'
    flow=list(csv.DictReader(flow_path.open())) if flow_path.exists() else []
    if flow:assert sha(flow_path)==receipt['outputSHA256'][flow_path.name]
    return dict(events=len(events),limited=len(accepted),invalidBaseOrInput=len(failed),
        zeroSteps=sum(float(x['alpha'])==0 for x in accepted),
        uniqueLocalPoints=len(set(x['local_point'] for x in events)),globalPositionsKnown=False,
        nativeUpdateRelativeError=error,energyReplayAbsoluteErrorJkg=energy_error,independentEnergyError=oracle_error,
        failedEvents=records,flowFallbackEvents=len(flow),flowInvalidOldStates=sum(x['recovered']=='0' for x in flow),
        eventLimitReached=len(events)==65536,conservationProven=False)


def main():
    out=BASE/'review01';out.mkdir(parents=True,exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes());started=time.monotonic()
    receipt=read(CASE/'engine_process.json');manifest=read(CASE/'case_manifest.json')
    assert receipt['status']=='terminal' and all(receipt['inputsStable'].values())
    normal=receipt['exitCode']==0 and receipt['checkpoints']['normalScheduleComplete']
    entries=receipt['checkpoints']['complete'];steps=[e['localIteration'] for e in entries]
    assert steps==[0,100,199] if normal else all(s in (0,100) for s in steps)
    core=IdealAir();trace=trace_review(out,receipt,core);save(out/'trace.json',trace)
    thermal={}
    for step in steps:
        thermal[str(step)]=audit(CASE,CASE/('audit%03d'%step),checkpoint=None if normal and step==199 else step)['errors']
    volume_path=ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'
    original_path=ROOT/'outputs/v56_validation/state_validity/state_validity_actual_rocket.mat'
    vol=loadmat(str(volume_path),squeeze_me=True)
    original=loadmat(str(original_path),squeeze_me=True,struct_as_record=False)['result']
    V,F,P=original.V,original.F.astype(int)-1,vol['vertices']
    assert np.array_equal(V,P[:len(V)]) and np.array_equal(F,vol['sourceFaces'].astype(int)-1)
    assert (len(F),len(P))==(19752,335726)
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    center=V[F].mean(axis=1);base=original.localGeometry.isAftTipAdjacent.astype(bool)
    nose=original.heating.noseFace.astype(bool)&~base;front=(center[:,0]>=-.6)&~base&~nose
    groups=dict(nose=nose,forward_side=front,rear_side=~base&~nose&~front,base=base)
    markers=boundary_faces(CASE/'mesh.su2');far,wall=markers['farfield'],markers['wall']
    cells=[vol['wedges'].astype(int)-1,vol['tetrahedra'].astype(int)-1]
    assert sum(len(g) for g in cells)==manifest['cells']==705490
    operators={}
    for i,adjacent in neighbors_of(cells,np.unique(far),len(P)).items():
        js=np.array(sorted(adjacent),int);delta=P[js]-P[i];weighted=delta/np.sum(delta**2,axis=1)[:,None]
        operators[i]=(js,np.linalg.solve(weighted.T@delta,weighted.T))
    bq=P[np.concatenate((far,wall))];bav=-np.cross(bq[:,1]-bq[:,0],bq[:,2]-bq[:,0])/2
    volume_error=abs(np.einsum('ij,ij->',bav,bq.mean(axis=1))/(3*manifest['fluidVolume'])-1)
    assert volume_error<2e-10
    fq=P[far];fav=-np.cross(fq[:,1]-fq[:,0],fq[:,2]-fq[:,0])/2
    geometry=wall_geometry(CASE/'mesh.su2');wallids=geometry['wallIDs']
    points,raw_cells,_=mesh(CASE/'mesh.su2');assert np.array_equal(points,P)
    dv=volumes(points,raw_cells);assert abs(sum(dv)/manifest['fluidVolume']-1)<2e-10
    bounds=core.state('PT',101325,[60,2000])['e']
    paths=[Path(__file__).resolve(),volume_path,original_path,CASE/'engine_process.json',CASE/'case_manifest.json',
           CASE/'mesh.su2',CASE/'history.csv',ROOT/'docs/validation/v167/PLAN.md']
    helpers=['scripts/reviewSU2ConservativeStepRun.py','scripts/analyzeSU2DualVolumes.py',
        'scripts/reviewSU2SSTLocalDiffusion.py','scripts/reviewSU2AdmissibleRun.py',
        'scripts/prepareSU2WallDistanceState.py','scripts/auditSU2IdealAirState.py',
        'scripts/su2SSTInterface.py','scripts/idealAirProperties.py','scripts/validateIdealAir.py']
    for name in helpers:
        p=ROOT/name;paths.append(p);target=out/'executed_sources'/name
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(p.read_bytes())
    locks={str(p):sha(p) for p in paths};locks.update(core.locks)
    frames=[];surfaces=[]
    for entry in entries:
        path=CASE/entry['volume'];assert sha(path)==entry['volumeSHA256'];locks[str(path)]=sha(path)
        _,f,nd=load_state(path,True,captured=True);assert nd==3
        rho,T,p=f['Density'],f['Temperature'],f['Pressure']
        u=np.column_stack([f['Momentum_'+a] for a in 'xyz'])/rho[:,None];thermo=core.state('PT',p,T)
        surface={a:f[b][F].mean(axis=1) for a,b in [('p','Pressure'),('q','Heat_Flux'),('T','Temperature')]};surfaces.append(surface)
        regions={}
        for name,mask in groups.items():
            regions[name]=dict(faces=int(sum(mask)),areaM2=float(sum(area[mask])),
                pressureRangePa=[float(min(surface['p'][mask])),float(max(surface['p'][mask]))],
                heatRangeWm2=[float(min(surface['q'][mask])),float(max(surface['q'][mask]))],
                meanWallFluidTemperatureRangeK=[float(min(surface['T'][mask])),float(max(surface['T'][mask]))],
                physicalHeatPowerW=float(area[mask]@surface['q'][mask]),physicsQualified=False)
        values=np.column_stack((u,T));grad=np.zeros((len(P),4,3))
        for i,(js,op) in operators.items():grad[i]=(op@(values[js]-values[i])).T
        mu=f['Laminar_Viscosity']+f['Eddy_Viscosity']
        stress=mu[:,None,None]*(grad[:,:3]+grad[:,:3].transpose(0,2,1))
        stress-=(2*mu*np.trace(grad[:,:3],axis1=1,axis2=2)/3)[:,None,None]*np.eye(3)
        conductivity=thermo['cp']*(f['Laminar_Viscosity']/.72+f['Eddy_Viscosity']/.9)
        mass_faces=np.einsum('ij,ij->i',fav,(rho[:,None]*u)[far].mean(axis=1))
        ef=(f['Energy']+p)[:,None]*u-np.einsum('nij,nj->ni',stress,u)-conductivity[:,None]*grad[:,3]
        energy_faces=np.einsum('ij,ij->i',fav,ef[far].mean(axis=1));power=float(area@surface['q'])
        mi,ei=float(sum(np.maximum(-mass_faces,0))),float(sum(np.maximum(-energy_faces,0)))
        assert mi>0 and ei>0
        balance=dict(massRelativeImbalance=float(abs(sum(mass_faces))/mi),
            energyRelativeImbalance=float(abs(sum(energy_faces)+power)/ei),
            massNetKgS=float(sum(mass_faces)),energyNetW=float(sum(energy_faces)+power),
            method='v166 physical WLS boundary diagnostic, OTH including captured k',nativeConservationProven=False)
        error=T[F]-300;square_mean=(np.sum(error,axis=1)**2+np.sum(error**2,axis=1))/12
        # 独立按保存后的k计算靠界数量/体积，不用事件次数代替体积。
        energy=f['Energy']/rho-np.sum(u*u,axis=1)/2-f['Turb_Kin_Energy']
        near=(energy-bounds[0]<1.)|(bounds[1]-energy<1.)
        frame=dict(localIteration=entry['localIteration'],seriesSuccessfulUpdates=entry['seriesSuccessfulUpdates'],
            givenWallTemperatureK=300,wallFluidTemperatureRangeK=[float(min(T[wallids])),float(max(T[wallids]))],
            maximumWallTemperatureErrorK=float(max(abs(T[wallids]-300))),
            areaRMSWallTemperatureErrorK=float(np.sqrt(area@square_mean/sum(area))),
            yPlusP99=float(np.quantile(f['Y_Plus'][wallids],.99)),physicalSurfaceHeatPowerW=power,
            boundaryBalanceDiagnostic=balance,regions=regions,postUpdateNearBoundEnergyThresholdJkg=1.,
            postUpdateNearBoundPoints=int(sum(near)),postUpdateNearBoundVolumeM3=float(sum(dv[near])),
            postUpdateOutsideBoundsPoints=int(sum((energy<bounds[0])|(energy>bounds[1]))),physicsQualified=False)
        frames.append(frame);save(out/('frame_%06d.json'%entry['localIteration']),frame)
    keys,h=history(CASE/'history.csv');ck,c=history(ROOT/'outputs/v166_validation/real02/history.csv')
    assert keys==ck and np.isfinite(h).all() and len(h)<=200
    common={n:float(np.mean(h[160:180,keys.index(n)]-c[160:180,keys.index(n)])) if len(h)>=180 else None for n in RESIDUALS}
    mean={w:{n:float(np.mean(rows[:,keys.index(n)])) for n in RESIDUALS} for w,rows in
          [('first20',h[:20]),('previous20',h[-40:-20]),('last20',h[-20:])]}
    delta={n:mean['last20'][n]-mean['first20'][n] for n in RESIDUALS}
    tail={n:mean['last20'][n]-mean['previous20'][n] for n in RESIDUALS}
    changes={}
    for field in ('p','q'):
        def change(a,b,updates):return float(area@abs(surfaces[b][field]-surfaces[a][field])/(area@abs(surfaces[a][field]))*100/updates)
        changes[field]=dict(first100=change(0,1,100) if len(frames)>1 else None,last99Per100=change(1,2,99) if normal else None)
    power=h[:,keys.index('HF')];prange=dict(first20=float(np.ptp(power[:20])/np.mean(abs(power[:20]))),last20=float(np.ptp(power[-20:])/np.mean(abs(power[-20:]))))
    gases=('rms[Rho]','rms[RhoE]');turb=('rms[k]','rms[w]');first,last=frames[0],frames[-1]
    gates=dict(normal200Updates=bool(normal and len(h)==200),strictStateAndSource=True,originalGeometry=True,
        noInvalidOrTruncatedEvents=trace['invalidBaseOrInput']==0 and trace['flowInvalidOldStates']==0 and not trace['eventLimitReached'],
        commonFlow=all(common[n]<=.05 for n in gases) if len(h)>=180 else None,
        commonTurbulence=all(common[n]<=.1 for n in turb) if len(h)>=180 else None,
        densityEnergyDevelopment=all(delta[n]<=-.3 for n in gases) if normal else None,
        turbulenceDevelopment=all(delta[n]<=.1 for n in turb) if normal else None,
        flowLateWindow=(all(tail[n]<=.05 for n in gases) and any(tail[n]<=-.05 for n in gases)) if normal else None,
        turbulenceLateWindow=all(tail[n]<=.1 for n in turb) if normal else None,
        pressureStability=changes['p']['last99Per100']<=.9*changes['p']['first100'] if normal else None,
        heatStability=changes['q']['last99Per100']<=.9*changes['q']['first100'] if normal else None,
        powerStability=prange['last20']<=.9*prange['first20'] if normal else None,
        maximumWallTemperature=last['maximumWallTemperatureErrorK']<=first['maximumWallTemperatureErrorK']+1 if normal else None,
        rmsWallTemperature=last['areaRMSWallTemperatureErrorK']<=.95*first['areaRMSWallTemperatureErrorK'] if normal else None,
        yPlusNoRegression=last['yPlusP99']<=1.05*first['yPlusP99'] if normal else None)
    for field in ('mass','energy'):
        key=field+'RelativeImbalance'
        gates[field+'BalanceNoRegression']=last['boundaryBalanceDiagnostic'][key]<=first['boundaryBalanceDiagnostic'][key]*1.1+1e-5 if normal else None
    assert all(sha(Path(p))==v for p,v in locks.items())
    report=dict(passedChecks=True,developmentBenefitPassed=all(v is True for v in gates.values()),gates=gates,
        thermalErrors=thermal,trace=trace,frames=frames,common160to179LogDifference=common,
        lastMinusFirst20Log=delta,lastMinusPrevious20Log=tail,surfaceChangePer100Updates=changes,
        nativePowerRelativeRanges=prange,finalResiduals={n:float(h[-1,keys.index(n)]) for n in RESIDUALS},
        originalMinus8Passed=all(h[-1,keys.index(n)]<=-8 for n in ('rms[Rho]','rms[RhoE]','rms[k]','rms[w]')),
        nativeNormalTerminal=bool(normal),finiteHistoryUpdates=len(h),completeCheckpoints=steps,priorUpdates=200,
        failureInstantFullStateAvailable=False,originalFaces=len(F),originalAreaM2=float(sum(area)),
        boundaryVolumeRelativeError=volume_error,inputSHA256=locks,inputsStable=True,elapsedSeconds=time.monotonic()-started,
        automaticAdditionalTasks=0,physicsQualified=False,newQualifiedArea=0)
    save(out/'report.json',report)
    print(json.dumps({k:report[k] for k in ('developmentBenefitPassed','gates','trace','finiteHistoryUpdates','completeCheckpoints','elapsedSeconds')},indent=2))


if __name__=='__main__':main()
