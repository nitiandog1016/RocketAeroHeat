"""按事前门槛比较完整箭体下限口袋初猜；短程改善不自动续跑或授予资格。"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha
from su2SSTInterface import state
from prepareSU2SSTPocketStudy import floor_nodes
from auditSU2SSTState import vtk_tokens


def main():
    base=ROOT/'outputs/v128_validation';target=base/'pocket_comparison.json';assert not target.exists()
    source=ROOT/'outputs/v127_validation/control'
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s']
    V=s.V;F=s.F.astype(int)-1;area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
    init=read(base/'pocket_repaired_roundoff/pocket_initialization.json');pocket=np.array(init['nodeIDs']);old111=np.array(init['exactFloorNodeIDs'])
    variants={};heat={};locks={};res=['rms[Rho]','rms[RhoE]','rms[k]','rms[w]'];profiles={}
    for name in ['control','repaired']:
        folder=base/('pocket_control' if name=='control' else 'pocket_repaired_roundoff')
        r=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert audit['interfacePassed'] and audit['checkpointLocalIteration']==9
        a,keys,_,_,_,_=state(folder/'restart_flow.csv',True,True)
        with (folder/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
        h=dict(zip([v.strip() for v in rows[0]],np.array(rows[1:],float).T))
        assert np.array_equal(h['Inner_Iter'],np.arange(10))
        omega=a[:,keys.index('Omega')];initial=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1,usecols=10)
        xy,field=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xy);assert max(d)<1e-6 and len(set(ids))==len(V)
        q=np.empty(len(V));q[ids]=field['Heat_Flux'][:,0];heat[name]=q[F].mean(axis=1)
        variants[name]=dict(elapsedSeconds=r['elapsedSeconds'],last5MeanLogResidual={k:float(np.mean(h[k][-5:])) for k in res},
            finalLogResidual={k:float(h[k][-1]) for k in res},history={k:h[k].tolist() for k in res+['HF','Min CFL','Max CFL']},
            floorEquivalentNodeCount=len(floor_nodes(omega)),repairedNodeCount=len(pocket),
            recoveredAbove100FloorCount=int(sum(omega[pocket]>.01)),old111RecoveredCount=int(sum(omega[old111]>.01)),
            changedOmegaAtPocketCount=int(sum(omega[pocket]!=initial[pocket])),pocketFinalOmegaRange=[float(min(omega[pocket])),float(max(omega[pocket]))],
            nativeHeatPowerW=audit['nativeHeatPowerW'],originalFaceHeatPowerW=float(area@heat[name]),
            sourceToFinalHeatL1=audit['diagnosticSourceToOutputHeatRelativeL1Change'],yPlusP99=audit['yPlusPercentiles']['p99'],
            originalResidualPassed=audit['originalResidualThresholdPassed'])
        # 固定原壁节点2的一列，仅作尾缘停滞/恢复剖面，不解释为一维真实解。
        ii=np.arange(12)*9878+2;profiles[name]=dict(nodeIDs=ii.tolist(),omega=omega[ii].tolist(),initialOmega=initial[ii].tolist(),
            distanceAlongColumnM=np.linalg.norm(a[ii,1:4]-a[2,1:4],axis=1).tolist())
        for p in [folder/'engine_process.json',folder/'sst_state_audit.json',folder/'history.csv',folder/'restart_flow.csv',
            folder/'sst_initialization.json',folder/'surface_flow.vtk']:locks[str(p)]=sha(p)
    c=variants['control'];n=variants['repaired'];gain={k:c['last5MeanLogResidual'][k]-n['last5MeanLogResidual'][k] for k in res}
    heatL1=float(area@abs(heat['repaired']-heat['control'])/(area@abs(heat['control'])))
    gates=dict(omegaMeanImprovedOneLog=gain['rms[w]']>=1,
        densityEnergyKNotWorse03=min(gain[k] for k in res[:3])>=-.3,
        atLeast95PercentPocketRecovered=n['recoveredAbove100FloorCount']/len(pocket)>=.95,
        noAdditionalFloorNodes=n['floorEquivalentNodeCount']<=c['floorEquivalentNodeCount'],
        pairedHeatL1Below5Percent=heatL1<=.05,yPlusP99Below1=n['yPlusP99']<=1)
    report=dict(variants=variants,meanResidualGainLog=gain,screeningGates=gates,allScreeningGatesPassed=all(gates.values()),
        pairedFaceHeatRelativeL1Difference=heatL1,originalFaces=len(F),originalAreaM2=float(sum(area)),columnProfiles=profiles,
        originalExactFloorCount=111,roundoffEquivalentFloorCount=146,physicsQualified=False,newQualifiedArea=0,
        defaultChanged=False,automaticContinuationAllowed=False,sourceSHA256=dict(locks,**{str(Path(__file__).resolve()):sha(Path(__file__).resolve())}))
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['sourceSHA256','columnProfiles']},indent=2))


if __name__=='__main__':main()
