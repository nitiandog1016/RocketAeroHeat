"""模块用途。独立核对回退事件、原生共同历史、物面稳定性；不把局部改变量称为全域守恒。"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree

from auditSU2IdealAirState import load_state
from idealAirProperties import IdealAir
from su2ThermalRuntime import ROOT, read, sha
from validateIdealAir import save, run, write_rows, read_rows

BASE = ROOT/'outputs/v164_validation'


def history(path):
    with path.open() as stream:
        rows = list(csv.reader(stream, skipinitialspace=True))
    return [key.strip() for key in rows[0]], np.array(rows[1:], float)


def main():
    out = BASE/'review01'
    out.mkdir(parents=True, exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    folder, old = BASE/'real01', ROOT/'outputs/v163_validation/real02'
    receipt = read(folder/'engine_process.json')
    assert receipt['exitCode'] == 0 and all(receipt['inputsStable'].values())
    assert receipt['checkpoints']['normalScheduleComplete']
    assert [r['localIteration'] for r in receipt['checkpoints']['complete']] == [0,49]
    paths = [folder/'engine_process.json', folder/'case.cfg', folder/'admissibility.csv',
             folder/'history.csv', old/'history.csv', folder/'audit49/report.json']
    assert sha(folder/'admissibility.csv') == receipt['outputSHA256']['admissibility.csv']
    with (folder/'admissibility.csv').open() as stream:
        events = list(csv.DictReader(stream))
    assert len(events) == 14 and [int(e['event']) for e in events] == list(range(1,15))
    assert all(e['dimension']=='3' and e['recovered']=='1' for e in events)
    def array(key):
        return np.array([float(e[key]) for e in events])
    k=array('k');trial=np.column_stack([array('trial_'+x) for x in ('rho','rhou','rhov','rhow','rhoE')])
    restored=np.column_stack([array('old_'+x) for x in ('rho','rhou','rhov','rhow','rhoE')])
    def thermal(a):
        return a[:,4]/a[:,0]-.5*np.sum((a[:,1:4]/a[:,0,None])**2,axis=1)-k
    trial_e, old_e = thermal(trial), thermal(restored)
    replay=max(float(np.max(abs(trial_e-array('trial_e')))),float(np.max(abs(old_e-array('old_e')))))
    assert replay < 1e-8 and min(trial[:,0])>0 and min(restored[:,0])>0
    core=IdealAir()
    endpoints=core.state('PT',101325,np.array([60.,2000.]))
    recovered=core.state('rhoe',restored[:,0],old_e)
    assert np.all(trial_e < endpoints['e'][0]) and np.all(old_e>=endpoints['e'][0])
    helper=ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    assert sha(helper)==read(ROOT/'outputs/v153_validation/native03/runtime.json')['sha256']
    write_rows(out/'oracle_input.csv',[('rhoT',r,t) for r,t in zip(restored[:,0],recovered['T'])])
    run(out,'oracle',[str(helper),'--oracle',str(out/'oracle_input.csv'),str(out/'oracle.csv')],30)
    oracle=read_rows(out/'oracle.csv');assert np.all(oracle['accepted']==1)
    oracle_error=float(np.max(abs(oracle['e']-old_e)/np.maximum(1,abs(old_e))))
    assert oracle_error<2e-11
    delta=restored-trial
    # 没有RCM全局映射、节点控制体积及调用迭代，不能给出位置或kg/J积分。
    trace=dict(events=len(events),uniqueLocalPoints=len(set(int(e['local_point']) for e in events)),
        recoveredEvents=14,failedOldStates=0,allTrialsBelow60KDomain=True,
        localPointIDs=sorted(set(int(e['local_point']) for e in events)),
        globalPositionsKnown=False,iterationPerEventKnown=False,
        minimumEnergyAt60K=float(endpoints['e'][0]),trialEnergyRangeJkg=[float(trial_e.min()),float(trial_e.max())],
        oldEnergyRangeJkg=[float(old_e.min()),float(old_e.max())],
        recoveredTemperatureRangeK=[float(recovered['T'].min()),float(recovered['T'].max())],
        currentKRange=[float(k.min()),float(k.max())],energyReplayMaximumJkg=replay,oracleEnergyError=oracle_error,
        massDensityChangeRangeKgm3=[float(delta[:,0].min()),float(delta[:,0].max())],
        totalEnergyDensityChangeRangeJm3=[float(delta[:,4].min()),float(delta[:,4].max())],
        maxRelativeDensityChange=float(np.max(abs(delta[:,0]/trial[:,0]))),
        maxRelativeTotalEnergyDensityChange=float(np.max(abs(delta[:,4]/trial[:,4]))),
        globallyIntegratedConservationChange=None,
        conservationNote='eventwise density changes only; no volumes/global mapping/iteration, not integrated mass or energy balance')
    hk,h=history(folder/'history.csv');ok,o=history(old/'history.csv')
    assert hk==ok and len(h)==50 and len(o)==17
    residual_names=['rms[Rho]','rms[RhoU]','rms[RhoV]','rms[RhoW]','rms[RhoE]','rms[k]','rms[w]']
    residual_ix=[hk.index(key) for key in residual_names]
    common_error=float(np.max(abs(h[:17,residual_ix]-o[:,residual_ix])))
    assert common_error<2e-8
    locks={str(p):sha(p) for p in paths}
    frames={}
    for step in (0,49):
        path=folder/('restart_flow_%06d.csv'%step)
        assert sha(path)==receipt['outputSHA256'][path.name]
        frames[step]=load_state(path,True,captured=True)
        locks[str(path)]=sha(path)
    old_path=old/'restart_flow_000000.csv'
    old_data,_,_=load_state(old_path,True,captured=True)
    first_difference=float(np.max(abs(frames[0][0]-old_data)/np.maximum(1,abs(old_data))))
    assert first_difference<2e-10
    locks[str(old_path)]=sha(old_path)
    surface_path=ROOT/'outputs/v118_validation/actual01/surface.mat'
    s=loadmat(str(surface_path),squeeze_me=True,struct_as_record=False)['s']
    vertices,faces=s.V,s.F.astype(int)-1
    areas=np.linalg.norm(np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]]),axis=1)/2
    distance,ids=cKDTree(frames[49][0][:,1:4]).query(vertices)
    assert distance.max()<1e-13 and len(set(ids))==len(vertices)
    p0=frames[0][1]['Pressure'][ids][faces].mean(axis=1)
    p1=frames[49][1]['Pressure'][ids][faces].mean(axis=1)
    q0=frames[0][1]['Heat_Flux'][ids][faces].mean(axis=1)
    q1=frames[49][1]['Heat_Flux'][ids][faces].mean(axis=1)
    native_power=h[:,hk.index('HF')]
    stability=dict(faces=len(faces),areaM2=float(areas.sum()),
        pressureL1Change=float(areas@abs(p1-p0)/(areas@abs(p0))),
        heatL1Change=float(areas@abs(q1-q0)/(areas@abs(q0))),
        physicalSurfacePowerW=[float(areas@q0),float(areas@q1)],
        nativeDiscretePowerW=[float(native_power[0]),float(native_power[-1])],
        last20NativePowerRelativeRange=float(np.ptp(native_power[-20:])/np.mean(abs(native_power[-20:]))),
        originalResidualThresholdPassed=False,
        finalResiduals={name:float(h[-1,hk.index(name)]) for name in residual_names})
    assert len(faces)==19752 and abs(areas.sum()-1.9646292566)<1e-9
    locks[str(surface_path)]=sha(surface_path)
    locks.update(core.locks)
    locks[str(helper)]=sha(helper);locks[str(Path(__file__).resolve())]=sha(Path(__file__).resolve())
    save(out/'report.json',dict(passed=True,trace=trace,common17ResidualLogDifference=common_error,
        firstFrameRelativeDifference=first_difference,stability=stability,thermal=read(folder/'audit49/report.json')['errors'],
        inputSHA256=locks,physicsQualified=False,newQualifiedArea=0))
    print(json.dumps(dict(trace=trace,common17ResidualLogDifference=common_error,stability=stability),indent=2))


if __name__=='__main__':
    main()
