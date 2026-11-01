"""复用一次完整域原生更新记录，定位触底/冻结及鼻部能量推进；不增加求解。"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.io import loadmat
from prepareSU2UpdateTraceCase import ROOT,BASE,PLAN,read,sha,config,write
from reviewSU2UpdateTrace import replay

def main():
    folder=BASE/'actual';target=BASE/'native_update_analysis.json';assert not target.exists()
    receipt=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
    assert audit['interfacePassed'] and audit['checkpointLocalIteration']==0
    result,data,a,keys=replay(folder);n=len(a)
    geometry=ROOT/'outputs/v145_validation/mesh/hybrid_volume.mat'
    mesh=loadmat(str(geometry),variable_names=['vertices','sourceVertices','sourcePreparedNodeIDsOneBased'])
    points=mesh['vertices'];nw=len(mesh['sourceVertices']);assert n==335726 and nw==9878
    assert np.all(abs(a[:,1:4]-points)<=8*np.finfo(float).eps*np.maximum(1,abs(points)))
    off=np.arange(n)>=nw
    s=data['sst'];f=data['flow'];floor=1e-10*(1+8*np.finfo(float).eps)
    oldfloor=(s['old'][:,0]<=floor)&off;newfloor=(s['after'][:,0]<=floor)&off
    zero=(s['alpha']==0)&off;newonly=newfloor&~oldfloor;recovered=oldfloor&~newfloor
    limiting=np.argmax(s['ratios'],axis=1)
    with (folder/'history.csv').open() as stream:rows=list(csv.reader(stream,skipinitialspace=True))
    h=dict(zip([v.strip() for v in rows[0]],map(float,rows[1])));assert len(rows)==2 and h['Inner_Iter']==0
    names=['Density','Momentum_x','Momentum_y','Momentum_z','Energy','TKE','Omega']
    histories=['rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]','rms[k]','rms[w]']
    rhs=np.column_stack([a[:,keys.index('Residual_'+q)] for q in names]);rms=np.log10(np.sqrt(np.mean(rhs**2,axis=0)))
    error=float(max(abs(rms-[h[k] for k in histories])));assert error<2e-8
    dt=a[:,keys.index('Delta_Time')];tdt=a[:,keys.index('Turb_Delta_Time')];assert min(dt)>0 and np.array_equal(dt,tdt)
    cfg=config(folder/'case.cfg');R=float(cfg['GAS_CONSTANT']);gamma=float(cfg['GAMMA_VALUE']);Ti=float(cfg['FREESTREAM_TEMPERATURE'])
    cp=gamma*R/(gamma-1);speed=float(cfg['MACH_NUMBER'])*np.sqrt(gamma*R*Ti);hi=cp*Ti+.5*speed**2
    velocity=np.column_stack([a[:,keys.index('Velocity_'+q)] for q in 'xyz']);T=a[:,keys.index('Temperature')]
    enthalpy=cp*T+.5*np.sum(velocity**2,axis=1)
    source=read(BASE.parent/'v146_validation/startup_comparison.json')['variants']['control']
    prepared=mesh['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1
    layers=np.where((prepared>=0)&(prepared<33*nw),prepared//nw,-1)
    def summary(mask):
        ids=np.flatnonzero(mask);values,counts=np.unique(layers[mask],return_counts=True)
        return dict(nodes=len(ids),zeroRelaxation=int(sum(mask&zero)),kControlsRelaxation=int(sum(mask&(limiting==0))),
            omegaControlsRelaxation=int(sum(mask&(limiting==1))),positiveRawKIncrement=int(sum(mask&(s['delta'][:,0]>0))),
            negativeRawKIncrement=int(sum(mask&(s['delta'][:,0]<0))),zeroRawKIncrement=int(sum(mask&(s['delta'][:,0]==0))),
            positiveIntegratedKRHS=int(sum(mask&(rhs[:,5]>0))),negativeIntegratedKRHS=int(sum(mask&(rhs[:,5]<0))),
            proposedKBelowFloor=int(sum(mask&(s['proposed'][:,0]<1e-10))),
            xRangeM=[float(min(points[mask,0])),float(max(points[mask,0]))] if len(ids) else None,
            referenceLayerCounts={str(int(k)):int(v) for k,v in zip(values,counts)})
    def sample(i):
        return dict(nodeZeroBased=int(i),xyzM=points[i].tolist(),referenceLayer=int(layers[i]),
            nativeDtS=float(dt[i]),nextCFL=float(a[i,keys.index('CFL')]),temperatureK=float(T[i]),primitiveEnthalpyRatio=float(enthalpy[i]/hi),
            flowCurrent=f['current'][i].tolist(),flowOld=f['old'][i].tolist(),flowRawIncrement=f['delta'][i].tolist(),
            flowRelaxation=float(f['alpha'][i]),flowAfter=f['after'][i].tolist(),
            turbulenceDensityUsed=float(s['rho'][i]),turbulenceCurrent=s['current'][i].tolist(),turbulenceOld=s['old'][i].tolist(),
            turbulenceRawIncrement=s['delta'][i].tolist(),turbulenceRelaxation=float(s['alpha'][i]),
            turbulenceUnclippedProposal=s['proposed'][i].tolist(),turbulenceAfter=s['after'][i].tolist(),
            controllingVariable='k' if limiting[i]==0 else 'omega',integratedRHS=dict(zip(names,rhs[i].tolist())))
    groups={'sourceKFloor':oldfloor,'remainingKFloor':oldfloor&newfloor,'recoveredKFloor':recovered,'newKFloor':newonly,
        'frozenOffWall':zero,'positiveKIncrementButFrozen':zero&(s['delta'][:,0]>0),
        'floorPositiveKIncrementButFrozen':oldfloor&zero&(s['delta'][:,0]>0),
        'floorNegativeKIncrement':oldfloor&(s['delta'][:,0]<0)}
    examples={name:[sample(i) for i in np.flatnonzero(mask)[:12]] for name,mask in groups.items() if any(mask)}
    hotspots={}
    for j,key in enumerate(names):
        w=rhs[:,j]**2;hotspots[key]=dict(largest=[sample(i) for i in np.argsort(w)[-3:][::-1]],
            sourceKFloorSquareFraction=float(sum(w[oldfloor])/sum(w)),newKFloorSquareFraction=float(sum(w[newonly])/sum(w)))
    locks=[Path(__file__).resolve(),ROOT/'scripts/reviewSU2UpdateTrace.py',PLAN,geometry,
        folder/'engine_process.json',folder/'sst_state_audit.json',folder/'restart_flow.csv',folder/'history.csv']
    locks+=list(folder.glob('*_update_*.csv'))
    out=dict(passed=True,sourceLocalIteration=99,newLocalIteration=0,nodes=n,originalFaces=19752,
        replay=result,sevenLogRMS=dict(zip(names,rms.tolist())),maximumLogRMSDifference=error,
        groups={name:summary(mask) for name,mask in groups.items()},examples=examples,hotspots=hotspots,
        sourceEnergyExtreme=sample(source['enthalpyExtremeNodeZeroBased']),currentEnergyExtreme=sample(int(np.argmax(enthalpy))),
        maximumTemperatureK=float(max(T)),maximumPrimitiveEnthalpyRatio=float(max(enthalpy)/hi),
        currentVsOldTurbulenceDifferentNodes=int(sum(np.any(s['current']!=s['old'],axis=1))),
        omegaFloorOffWall=int(sum(off&(s['after'][:,1]<=1e-4*(1+8*np.finfo(float).eps)))),
        semantics='Native linear increments and relaxation at each subsystem update; source terms not individually recorded. RHS is pre-update integrated right-hand side; saved CFL is adapted for next update.',
        referenceLayerMinus1MeansFanOrExterior=True,physicsQualified=False,newQualifiedArea=0,solverRuns=0,
        sourceSHA256={str(p):sha(p) for p in locks})
    write(target,out);print(json.dumps({k:out[k] for k in ['passed','replay','maximumLogRMSDifference','groups','maximumTemperatureK','maximumPrimitiveEnthalpyRatio']}))
if __name__=='__main__':main()
