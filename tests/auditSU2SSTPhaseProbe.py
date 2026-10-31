"""独立读取一次SST更新的输入/输出k，逐壁点比较两种阶段的热梯度。"""
from pathlib import Path
import argparse,json,sys
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import sha,read,config
from su2SSTInterface import state,wall_geometry
from auditSU2SSTState import vtk_tokens


def audit(folder):
    folder=Path(folder).resolve();r=read(folder/'engine_process.json');c=config(folder/'case.cfg')
    assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and not r['cancelled']
    assert all(r['inputsStable'].values()) and r['physicsModel']['diagnosticOnly']
    for name in ['mesh.su2','initial_flow.csv','case.cfg','sst_initialization.json','sst_contract.json']:
        assert sha(folder/name)==r['inputSHA256'][str(folder/name)]
    for name in ['restart_flow.csv','restart_flow_000000.csv','surface_flow.vtk','surface_flow_000000.vtk','history.csv']:
        assert sha(folder/name)==r['outputSHA256'][name]
    old,keys,nd,_,_,_=state(folder/'initial_flow.csv',True)
    new,keys2,nd2,rho,u,ecurrent=state(folder/'restart_flow.csv',True)
    assert keys==keys2 and nd==nd2 and np.array_equal(old[:,:nd+1],new[:,:nd+1])
    kcol=keys.index('Turb_Kin_Energy');stage=ecurrent+new[:,kcol]-old[:,kcol]
    gamma=float(c['GAMMA_VALUE']);R=float(c['GAS_CONSTANT']);factor=(gamma-1)/R
    T=stage*factor;Tcurrent=ecurrent*factor
    g=wall_geometry(folder/'mesh.su2');xyz,f=vtk_tokens(folder/'surface_flow.vtk')
    d,ids=cKDTree(new[:,1:4]).query(xyz);assert max(d)<1e-6
    order={v:i for i,v in enumerate(g['wallIDs'])};ix=np.array([order[j] for j in ids]);near=g['normalNeighbor'][ix]
    def mu(t):return float(c['MU_REF'])*(t/float(c['MU_T_REF']))**1.5*(float(c['MU_T_REF'])+float(c['SUTHERLAND_CONSTANT']))/(t+float(c['SUTHERLAND_CONSTANT']))
    cp=gamma*R/(gamma-1);conduct=mu(T[ids])*cp/float(c['PRANDTL_LAM'])
    q=conduct*(T[near]-300)/g['thermalDistance'][ix]
    qcurrent=mu(Tcurrent[ids])*cp/float(c['PRANDTL_LAM'])*(Tcurrent[near]-300)/g['thermalDistance'][ix]
    error=abs(q-f['Heat_Flux'][:,0])/np.maximum(abs(q),1)
    currentError=abs(qcurrent-f['Heat_Flux'][:,0])/np.maximum(abs(qcurrent),1)
    eos=max(float(max(abs(T[ids]/f['Temperature'][:,0]-1))),float(max(abs(rho[ids]*R*T[ids]/f['Pressure'][:,0]-1))))
    hist=np.atleast_1d(np.genfromtxt(str(folder/'history.csv'),names=True,delimiter=','));assert len(hist)==1 and hist['Inner_Iter'][0]==0
    native=float(g['nodeArea'][ix]@f['Heat_Flux'][:,0]);powerError=abs(native-float(hist['HF'][-1]))/max(abs(native),1)
    worst=np.argsort(currentError)[-10:]
    report=dict(nativeExitCode=r['exitCode'],nativeElapsedSeconds=r['elapsedSeconds'],nodes=len(new),wallNodes=len(ids),
        sourceLocalIteration=99,diagnosticLocalIteration=0,uninterruptedTrajectoryEquivalent=False,
        threshold=2e-6,inputKHeatGradientError=float(max(error)),updatedKHeatGradientError=float(max(currentError)),
        inputKMaximumAbsoluteHeatErrorWm2=float(max(abs(q-f['Heat_Flux'][:,0]))),
        updatedKMaximumAbsoluteHeatErrorWm2=float(max(abs(qcurrent-f['Heat_Flux'][:,0]))),
        nativeHeatPowerW=native,nativeHistoryPowerError=powerError,wallEOSError=eos,
        minimumDensity=float(min(rho)),minimumStageThermalEnergy=float(min(stage)),
        maximumKThermalTemperatureChangeK=float(max(abs(T-Tcurrent))),
        nativeHighPrecisionTemperatureWritten=(folder/'volume_flow.csv').is_file(),
        phaseProvedForThisUpdate=bool(max(error)<2e-6 and max(currentError)>2e-6 and eos<2e-6 and powerError<2e-6),
        finalResiduals=[float(hist[k][-1]) for k in ['rmsRho','rmsRhoE','rmsk','rmsw']],
        originalResidualPassed=False,physicsQualified=False,newQualifiedArea=0,
        worstUpdatedK=[dict(nodeID=int(ids[j]),neighborID=int(near[j]),savedHeat=float(f['Heat_Flux'][j,0]),
            inputKHeat=float(q[j]),updatedKHeat=float(qcurrent[j]),oldK=float(old[near[j],kcol]),newK=float(new[near[j],kcol]),
            inputKError=float(error[j]),updatedKError=float(currentError[j])) for j in worst],
        sourceSHA256={str(p):sha(p) for p in [Path(__file__),ROOT/'scripts/su2SSTInterface.py',folder/'engine_process.json']})
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    target=Path(a.output);assert not target.exists();report=audit(a.folder)
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','worstUpdatedK']},indent=2))
