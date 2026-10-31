"""按事前规则比较完整箭体近壁初猜；原面积、两次保存及原门槛保持。"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2CheckpointContinuation import ROOT,read,sha
from prepareSU2SSTColumnStudy import SOURCE,STUDY
from prepareSU2SSTPocketStudy import floor_nodes
from su2SSTInterface import state
from auditSU2SSTState import vtk_tokens


def main():
    base=ROOT/'outputs/v130_validation';target=base/'column_comparison.json';assert not target.exists()
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s']
    V=s.V;F=s.F.astype(int)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
    a,keys,_,_,_,_=state(SOURCE/'restart_flow.csv',True,True)
    pocket=floor_nodes(a[:,keys.index('Omega')]);assert len(pocket)==146
    profiles={};ii=np.arange(13)*9878+2
    def profile(data):
        return dict(nodeIDs=ii.tolist(),k=data[ii,keys.index('Turb_Kin_Energy')].tolist(),
                    omega=data[ii,keys.index('Omega')].tolist(),
                    distanceAlongColumnM=np.linalg.norm(data[ii,1:4]-data[2,1:4],axis=1).tolist())
    profiles['source']=profile(a);del a
    variants={};heat={};locks={str(SOURCE/'restart_flow.csv'):sha(SOURCE/'restart_flow.csv')}
    residuals=['rms[Rho]','rms[RhoE]','rms[k]','rms[w]']
    for name in ['control','column']:
        folder=base/('control' if name=='control' else 'column_final')
        r=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert audit['interfacePassed'] and not audit['diagnosticOnly'] and audit['checkpointLocalIteration']==99
        assert r['physicsModel']['initializationStudy']['studyID']==STUDY
        assert r['physicsModel']['initializationStudy']['variant']==name
        with (folder/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
        h=dict(zip([v.strip() for v in rows[0]],np.array(rows[1:],float).T))
        assert np.array_equal(h['Inner_Iter'],np.arange(100))
        snapshots={}
        for step in [50,99]:
            path=folder/('restart_flow_%06d.csv'%step)
            assert sha(path)==r['outputSHA256'][path.name]
            a,newkeys,_,_,_,_=state(path,True,True);assert keys==newkeys
            omega=a[:,keys.index('Omega')];floor=floor_nodes(omega)
            new=np.setdiff1d(floor,pocket);remaining=np.intersect1d(floor,pocket)
            snapshots[str(step)]=dict(floorEquivalentNodeCount=len(floor),
                recoveredOriginalNodesAbove100Floor=int(sum(omega[pocket]>.01)),
                originalFloorNodesRemaining=remaining.tolist(),newFloorNodeIDs=new.tolist(),
                originalNodeOmegaRange=[float(min(omega[pocket])),float(max(omega[pocket]))])
            profiles[name+'_'+str(step)]=profile(a);del a
            locks[str(path)]=sha(path)
        xy,field=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xy)
        assert max(d)<1e-6 and len(set(ids))==len(V)
        q=np.empty(len(V));q[ids]=field['Heat_Flux'][:,0];heat[name]=q[F].mean(axis=1)
        variants[name]=dict(elapsedSeconds=r['elapsedSeconds'],snapshots=snapshots,
            last20MeanLogResidual={k:float(np.mean(h[k][-20:])) for k in residuals},
            finalLogResidual={k:float(h[k][-1]) for k in residuals},
            history={k:h[k].tolist() for k in residuals+['HF','Min CFL','Max CFL']},
            nativeHeatPowerW=audit['nativeHeatPowerW'],originalFaceHeatPowerW=float(area@heat[name]),
            withinRunHeatRelativeL1Change=audit['withinRunHeatRelativeL1Change'],
            yPlusP99=audit['yPlusPercentiles']['p99'],originalResidualPassed=audit['originalResidualThresholdPassed'])
        for fn in ['engine_process.json','sst_state_audit.json','history.csv','sst_initialization.json','surface_flow.vtk']:
            locks[str(folder/fn)]=sha(folder/fn)
    c=variants['control'];n=variants['column']
    gain={k:c['last20MeanLogResidual'][k]-n['last20MeanLogResidual'][k] for k in residuals}
    heatL1=float(area@abs(heat['column']-heat['control'])/(area@abs(heat['control'])))
    initial=read(base/'column_final/column_initialization.json')
    gates=dict(omegaMeanImprovedOneLog=gain['rms[w]']>=1,
        densityEnergyKNotWorse03=min(gain[k] for k in residuals[:3])>=-.3,
        all146RecoveredAtBothSnapshots=all(v['recoveredOriginalNodesAbove100Floor']==146 for v in n['snapshots'].values()),
        noNewFloorPositionsAtEitherSnapshot=all(not v['newFloorNodeIDs'] for v in n['snapshots'].values()),
        pairedHeatL1Below5Percent=heatL1<=.05,yPlusP99Below1=n['yPlusP99']<=1,
        restartThermalEnergyPreserved=initial['thermalInternalEnergyMaxRelativeDifference']<=2e-12)
    for p in [base/'column_final/column_initialization.json',Path(__file__).resolve(),
              ROOT/'docs/validation/v130/PLAN.md',ROOT/'docs/validation/v130/NATIVE_FLOOR.md']:locks[str(p)]=sha(p)
    report=dict(variants=variants,meanResidualGainLog=gain,screeningGates=gates,allScreeningGatesPassed=all(gates.values()),
        pairedFaceHeatRelativeL1Difference=heatL1,originalFaces=len(F),originalAreaM2=float(sum(area)),
        originalFloorNodeIDs=pocket.tolist(),columnProfiles=profiles,initialization=initial,
        sourceLocalIteration=149,newLocalIteration=99,uninterruptedTrajectoryEquivalent=False,
        physicsQualified=False,newQualifiedArea=0,defaultChanged=False,automaticContinuationAllowed=False,sourceSHA256=locks)
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(gates=gates,gain=gain,heatL1=heatL1,
                         snapshots={k:v['snapshots'] for k,v in variants.items()}),indent=2))


if __name__=='__main__':main()
