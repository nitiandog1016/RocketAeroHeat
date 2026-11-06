"""模块用途。独立复核实际限步/几何/积分修正及同源收益，保留所有未过条件。"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from analyzeSU2DualVolumes import FACES
from auditSU2IdealAirState import audit,load_state
from idealAirProperties import IdealAir
from reviewSU2AdmissibleRun import history
from su2SSTInterface import mesh
from su2ThermalRuntime import ROOT,read,sha
from validateIdealAir import save,run,write_rows,read_rows

BASE=ROOT/'outputs/v165_validation'
FIELDS=('rho','rhou','rhov','rhow','rhoE')
RESIDUALS=('rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]','rms[k]','rms[w]')


def volumes(points,cells):
    """从固定SU2原连接独立构造中位对偶体积，不使用事件输出的体积作为输入。"""
    result=np.zeros(len(points))
    for kind,faces in FACES.items():
        count=max(max(face) for face in faces)+1
        ids=np.array([row[1:1+count] for row in cells if row[0]==kind],int)
        if not len(ids):continue
        q=points[ids];center=q.mean(axis=1)
        for face in faces:
            fc=q[:,face].mean(axis=1)
            for a,b in zip(face,face[1:]+face[:1]):
                midpoint=(q[:,a]+q[:,b])/2
                for j in (a,b):
                    p=q[:,j]
                    volume=abs(np.einsum('ij,ij->i',np.cross(midpoint-p,fc-p),center-p))/6
                    assert min(volume)>0
                    np.add.at(result,ids[:,j],volume)
    assert min(result)>0
    return result


def main():
    out=BASE/'review01';out.mkdir(parents=True,exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    case=BASE/'real01';control=ROOT/'outputs/v164_validation/real01'
    receipt=read(case/'engine_process.json');manifest=read(case/'case_manifest.json')
    assert receipt['exitCode']==0 and all(receipt['inputsStable'].values())
    assert receipt['checkpoints']['normalScheduleComplete']
    assert [x['localIteration'] for x in receipt['checkpoints']['complete']]==[0,49]
    thermal={}
    thermal['0']=audit(case,case/'audit00',checkpoint=0)
    thermal['49']=audit(case,case/'audit49')
    events=list(csv.DictReader((case/'conservative_steps.csv').open()))
    assert events and len(events)<=4096
    assert sha(case/'conservative_steps.csv')==receipt['outputSHA256']['conservative_steps.csv']
    assert [int(e['event']) for e in events]==list(range(1,len(events)+1))
    def array(key):return np.array([float(e[key]) for e in events])
    def states(prefix):return np.column_stack([array(prefix+'_'+v) for v in FIELDS])
    old,trial,accepted,linear=map(states,('old','trial','accepted','linear'))
    theta,k,alpha=array('theta'),array('k'),array('alpha')
    assert all(e['context_known']=='1' and e['old_valid']=='1' and e['dimension']=='3' and e['domain']=='1' for e in events)
    assert np.all((theta>0)&(theta<=.5))
    assert np.array_equal(theta,2.**(-array('attempts')))
    assert min(array('inner_iter'))>=0 and max(array('inner_iter'))<50 and np.all(array('mesh')==0)
    segment=float(np.max(abs(accepted-(old+theta[:,None]*(trial-old)))/np.maximum(1,abs(accepted))))
    assert segment<2e-13
    # 线性增量是实际读取量，单独检查；不先假定trial-old就是该增量。
    increment_error=float(np.max(abs((trial-old)-alpha[:,None]*linear)/np.maximum(1,abs(trial))))
    def e(u):return u[:,4]/u[:,0]-.5*np.sum((u[:,1:4]/u[:,0,None])**2,axis=1)-k
    old_e,trial_e,accepted_e=map(e,(old,trial,accepted))
    replay=max(float(np.max(abs(actual-array(name)))) for actual,name in
               ((old_e,'old_e'),(trial_e,'trial_e'),(accepted_e,'accepted_e')))
    assert replay<1e-8
    core=IdealAir();new_state=core.state('rhoe',accepted[:,0],accepted_e)
    old_state=core.state('rhoe',old[:,0],old_e)
    bounds=core.state('PT',101325,[60,2000])['e']
    invalid_trials=(trial[:,0]<=0)|(trial_e<bounds[0])|(trial_e>bounds[1])
    assert np.all(invalid_trials)
    helper=ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    assert sha(helper)==read(ROOT/'outputs/v153_validation/native03/runtime.json')['sha256']
    write_rows(out/'oracle_input.csv',[('rhoT',r,t) for r,t in zip(accepted[:,0],new_state['T'])])
    run(out,'oracle',[str(helper),'--oracle',str(out/'oracle_input.csv'),str(out/'oracle.csv')],30)
    oracle=read_rows(out/'oracle.csv');assert np.all(oracle['accepted']==1)
    oracle_error=float(np.max(abs(oracle['e']-accepted_e)/np.maximum(1,abs(accepted_e))))
    assert oracle_error<2e-11
    points,cells,wall=mesh(case/'mesh.su2');v=volumes(points,cells)
    assert len(points)==335726 and len(cells)==705490 and len(wall)==19752
    assert abs(v.sum()/manifest['fluidVolume']-1)<2e-10
    ids=array('global_point').astype(int);xyz=np.column_stack([array(a) for a in 'xyz'])
    location_error=float(np.max(abs(points[ids]-xyz)))
    volume_error=float(np.max(abs(v[ids]/array('volume')-1)))
    assert location_error<1e-13 and volume_error<2e-10
    hk,h=history(case/'history.csv');ck,c=history(control/'history.csv')
    assert hk==ck and len(h)==len(c)==50 and np.isfinite(h).all()
    changes={name:float(np.mean(h[-20:,hk.index(name)])-np.mean(c[-20:,hk.index(name)])) for name in RESIDUALS}
    source=ROOT/'outputs/v118_validation/actual01/surface.mat'
    s=loadmat(str(source),squeeze_me=True,struct_as_record=False)['s']
    V,F=s.V,s.F.astype(int)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    distance,wallids=cKDTree(points).query(V)
    assert distance.max()<1e-13 and len(set(wallids))==len(V)
    frames={step:load_state(case/('restart_flow_%06d.csv'%step),True,captured=True) for step in (0,49)}
    surface={}
    for step,(_,f,nd) in frames.items():
        assert nd==3
        surface[step]={a:f[b][wallids][F].mean(axis=1) for a,b in [('p','Pressure'),('q','Heat_Flux')]}
    stability={name+'L1Change':float(area@abs(surface[49][field]-surface[0][field])/(area@abs(surface[0][field])))
               for name,field in [('pressure','p'),('heat','q')]}
    power=h[:,hk.index('HF')]
    stability.update(last20NativePowerRelativeRange=float(np.ptp(power[-20:])/np.mean(abs(power[-20:]))),
        physicalSurfacePowerW=[float(area@surface[i]['q']) for i in (0,49)],
        nativeDiscretePowerW=[float(power[0]),float(power[-1])],
        finalResiduals={name:float(h[-1,hk.index(name)]) for name in RESIDUALS})
    reference=read(ROOT/'outputs/v164_validation/review01/report.json')['stability']
    heat_ratio=stability['heatL1Change']/reference['heatL1Change']
    range_ratio=stability['last20NativePowerRelativeRange']/reference['last20NativePowerRelativeRange']
    final=frames[49][1]
    integrals=np.column_stack([final[n] for n in ('Density','Momentum_x','Momentum_y','Momentum_z','Energy')])
    inventory=v@integrals
    correction=v[ids,None]*(accepted-trial)
    gates=dict(normalRun=True,completeLocatedPositiveSteps=True,strictThermal=True,
        originalGeometry=True,
        densityEnergyNoRegression=all(changes[n]<=.05 for n in ('rms[Rho]','rms[RhoE]')),
        densityOrEnergyBenefit=any(changes[n]<=-.10 for n in ('rms[Rho]','rms[RhoE]')),
        turbulenceNoRegression=all(changes[n]<=.10 for n in ('rms[k]','rms[w]')),
        heatStabilityNoRegression=heat_ratio<=1 and range_ratio<=1,
        heatStabilityBenefit=min(heat_ratio,range_ratio)<=.95)
    trace=dict(events=len(events),uniqueGlobalPoints=len(set(ids.tolist())),
        globalPointIDs=sorted(set(ids.tolist())),localPointIDs=sorted(set(array('local_point').astype(int).tolist())),
        eventIterationRange=[int(min(array('inner_iter'))),int(max(array('inner_iter')))],
        outputPhaseCounts={str(i):int(sum(array('output')==i)) for i in (0,1)},
        thetaRange=[float(min(theta)),float(max(theta))],alphaRange=[float(min(alpha)),float(max(alpha))],
        rejectedTrialEnergyRangeJkg=[float(min(trial_e)),float(max(trial_e))],
        acceptedTemperatureRangeK=[float(min(new_state['T'])),float(max(new_state['T']))],
        oldTemperatureRangeK=[float(min(old_state['T'])),float(max(old_state['T']))],
        coordinatesBoundsM=[xyz.min(axis=0).tolist(),xyz.max(axis=0).tolist()],
        maximumSegmentError=segment,maximumEnergyReplayJkg=replay,oracleEnergyError=oracle_error,
        maximumLocationErrorM=location_error,maximumVolumeRelativeError=volume_error,
        maximumObservedLinearIncrementError=increment_error,
        rawIncrementAgreement=increment_error<2e-13,
        integratedCorrection=correction.sum(axis=0).tolist(),absoluteIntegratedCorrection=abs(correction).sum(axis=0).tolist(),
        finalDomainInventory=inventory.tolist(),
        normalizedMassEnergyCorrection=[float(sum(correction[:,j])/inventory[j]) for j in (0,4)],
        nativeBoundaryFluxAvailable=False,globalConservationProven=False,
        correctionSemantics='sum over pseudo-time events V*(accepted-trial); kg, kg*m/s, J; not physical-time flux or steady conservation')
    paths=[case/'engine_process.json',case/'history.csv',case/'mesh.su2',case/'conservative_steps.csv',
           control/'history.csv',ROOT/'outputs/v164_validation/review01/report.json',source,helper,Path(__file__).resolve()]
    save(out/'report.json',dict(passed=True,numericalCandidateAccepted=all(gates.values()),gates=gates,
        trace=trace,stability=stability,controlStability=reference,meanLast20LogChanges=changes,
        heatStabilityRatios=dict(areaL1=heat_ratio,powerRange=range_ratio),
        originalResidualThresholdPassed=thermal['49']['originalResidualThresholdPassed'],
        yPlusP99=thermal['49']['yPlusP99'],thermalErrors={s:r['errors'] for s,r in thermal.items()},
        inputSHA256={str(p):sha(p) for p in paths},physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(dict(accepted=all(gates.values()),gates=gates,trace=trace,stability=stability,
                         meanLogChanges=changes,heatRatios=[heat_ratio,range_ratio]),indent=2))


if __name__=='__main__':main()
