"""模块职责：独立核对恢复五量、实际延拓方向、后续完整阶段及分层状态扰动。"""
from pathlib import Path
import csv
import numpy as np
from scipy.optimize import brentq
from phaseAwareDryAir import pressure_ceiling
from su2PhaseAirFields import inverse_energy
from fasRecoveryFieldsV248 import binary,state
from runSU2RocketPhaseV222 import read,save


def rows(path):
    with Path(path).open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit_events(case,out,air,check,completed,success):
    events=rows(case/'rejections.csv');observed=rows(case/'recovery_checks.csv');count=len(events)
    assert [int(e['eventID']) for e in events]==list(range(1,count+1))
    fine=state(case/'state_000000.bin',initial=True)
    coarse=binary(case/'coarse_initial.bin',b'RAH246LE',44484,33)
    low=air.H[0]-air.R*air.T[0];high=air.H[-1]-air.R*air.T[-1]
    def ucols(data,prefix):return np.array([[float(e[prefix+str(j)]) for j in range(5)] for e in data])
    def decode(u):
        vel=u[:,1:4]/u[:,0,None];energy=u[:,4]/u[:,0]-.5*np.sum(vel*vel,axis=1)-300000
        temp=inverse_energy(air,energy);p=u[:,0]*air.R*temp;phys=air.evaluate(p,temp)
        return energy,temp,p,phys,vel
    if not count:
        result=dict(passed=True,eventCount=0,firstReturnChecks=0,completeCycleChecks=0,physicalQualification=False)
        save(out/'recovery_report.json',result);return result
    trial=ucols(events,'trial');old=ucols(events,'old');returned=ucols(events,'returned')
    assert np.array_equal(old.view(np.uint64),returned.view(np.uint64))
    assert all(int(e['recovered'])==1 for e in events)
    eo,to,po,_,_=decode(old)
    velocity=trial[:,1:4]/trial[:,0,None];et=trial[:,4]/trial[:,0]-.5*np.sum(velocity*velocity,axis=1)-300000
    illegal=(trial[:,0]<=0)|(et<low)|(et>high)|~np.isfinite(trial).all(axis=1)
    available=~illegal
    if np.any(available):
        t=inverse_energy(air,et[available]);p=trial[available,0]*air.R*t;illegal[available]=(p<=0)|(p>pressure_ceiling(t))
    assert illegal.all(),'NATIVE_REJECTION_NOT_INDEPENDENTLY_REPRODUCED'
    by_event={i:[] for i in range(1,count+1)}
    for r in observed:by_event[int(r['eventID'])].append(r)
    volume=np.array([float(e['volume']) for e in events]);layer=np.array([int(e['level']) for e in events]);updates=np.array([int(e['update']) for e in events])
    matched=np.array([int(e['prolongationMatched'])==1 for e in events]);before=ucols(events,'prolongBefore');direction=ucols(events,'prolongDirection');after=ucols(events,'prolongAfter')
    if matched.any():
        factor=np.array([float(e['prolongationFactor']) for e in events])[matched];assert np.all(factor==.5)
        check('actualProlongationDirection',after[matched],before[matched]+factor[:,None]*direction[matched])
        assert np.array_equal(trial[matched].view(np.uint64),after[matched].view(np.uint64))
    assert np.isnan(before[~matched]).all() and np.isnan(direction[~matched]).all() and np.isnan(after[~matched]).all()
    delta=returned-trial;delta[:,4]-=300000*delta[:,0];delta*=volume[:,None]
    check('perEventVolumePhysicalDelta',ucols(events,'volumeDeltaPhysical'),delta)
    first=[];finished=[];details=[]
    for k,e in enumerate(events):
        lvl=layer[k];i=int(e['localID']);gid=int(e['globalID']);geometry=fine if lvl==0 else coarse
        assert lvl in (0,1) and int(geometry[i,0])==gid;check('eventVolume',volume[k],geometry[i,4])
        checks=by_event[k+1];f=[r for r in checks if r['checkpoint']!='complete_flow_SA'];c=[r for r in checks if r['checkpoint']=='complete_flow_SA']
        assert len(f)<=1 and len(c)<=1
        if f:
            a=np.array([float(f[0]['U%d'%j]) for j in range(5)]);assert np.array_equal(a.view(np.uint64),returned[k].view(np.uint64));first.append(f[0])
        if c:assert updates[k]<=completed;finished.append(c[0])
        if success:assert len(f)==len(c)==1
        details.append(dict(eventID=k+1,update=int(updates[k]),level=int(lvl),globalID=gid,xyzM=geometry[i,1:4].tolist(),
          wallDistanceM=float(geometry[i,24] if lvl==0 else geometry[i,32]),trialPhysicalEnergyJkg=float(et[k]),oldTemperatureK=float(to[k]),
          oldPressurePa=float(po[k]),prolongationMatched=bool(matched[k]),returnStageChecked=bool(f),completeCycleChecked=bool(c),
          volumeMassDeltaKg=float(delta[k,0]),volumePhysicalEnergyDeltaJ=float(delta[k,4])))
    # 第一个后续原生阶段的五量/物性是已恢复状态；粗末缓存不冒充更新后primitive。
    firstU=ucols(first,'U');ef,tf,pf,phys,vf=decode(firstU)
    check('restoredTemperature',[float(r['T']) for r in first],tf);check('restoredPressure',[float(r['p']) for r in first],pf)
    check('restoredH',[float(r['Hinternal'])-300000 for r in first],np.asarray(phys['enthalpyJkg'])+.5*np.sum(vf*vf,axis=1))
    transport=air.interpolated_transport(pf,tf);check('restoredK',[float(r['k']) for r in first],transport['conductivityWmK']);check('restoredCp',[float(r['Cp']) for r in first],phys['cpJkgK'])
    fineFinished=[r for r in finished if int(r['level'])==0]
    if fineFinished:
        _,tc,pc,phys,vc=decode(ucols(fineFinished,'U'));check('completeFineT',[float(r['T']) for r in fineFinished],tc);check('completeFineP',[float(r['p']) for r in fineFinished],pc)
        trans=air.interpolated_transport(pc,tc);rho=ucols(fineFinished,'U')[:,0];nu=np.array([float(r['SAcurrent']) for r in fineFinished]);chi=nu/(np.asarray(trans['dynamicViscosityPas'])/rho)
        mut=np.where(nu<0,0,rho*nu*chi**3/(chi**3+7.1**3));check('completeFineSAeddy',[float(r['mutSA']) for r in fineFinished],mut)
    chosen=np.unique(np.r_[np.linspace(0,count-1,min(16,count)).astype(int),np.argmin(eo),np.argmax(eo)])
    brent=[brentq(lambda t:float(air.thermal(t)[1]-air.R*t)-eo[i],60,2000,xtol=1e-10) for i in chosen];check('recoveryOldBrent',to[chosen],brent)
    totals={str(l):dict(events=int(sum(layer==l)),signedVolumeDelta=delta[layer==l].sum(axis=0).tolist(),absoluteVolumeDelta=abs(delta[layer==l]).sum(axis=0).tolist()) for l in (0,1)}
    if success:
        native=read(case/'recovery_summary.json');assert native['events']==count and native['firstReturnChecks']==len(first) and native['completeCycleChecks']==len(finished) and native['pendingEvents']==0
        check('nativeCombinedDiagnosticSum',native['physicalVolumeStateDisturbance'],delta.sum(axis=0));assert native['maximumProlongationScaledDifference']<=2e-14
    result=dict(passed=True,eventCount=count,firstReturnChecks=len(first),completeCycleChecks=len(finished),matchedProlongationEvents=int(sum(matched)),
      uniquePointsByLevel={str(l):len(set(int(e['globalID']) for e in events if int(e['level'])==l)) for l in (0,1)},
      initialEventUpdate=int(min(updates)),lastEventUpdate=int(max(updates)),trialBelowEnergyFloor=int(sum(et<low)),
      oldTemperatureRangeK=[float(min(to)),float(max(to))],stateDisturbanceByLevel=totals,
      combinedAcrossLevelsIsNotPhysicalBudget=True,allTrialRejectedIndependently=True,oldReturnedBitEqual=True,
      completeFineSAEddyChecked=len(fineFinished),brentSamples=len(chosen),physicalQualification=False)
    save(out/'recovery_details.json',details);save(out/'recovery_report.json',result);return result
