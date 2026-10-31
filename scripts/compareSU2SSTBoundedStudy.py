"""v131整箭受控有界输运的事前筛选；复用历史控制，不重跑、不改面积分母。"""
import csv,json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2SSTBoundedStudy import ROOT,BASE,REFERENCE,STUDY,sha,read,write,validate_source,config
from prepareSU2SSTPocketStudy import floor_nodes
from su2SSTInterface import state
from auditSU2SSTState import vtk_tokens
def main():
    target=BASE/'bounded_comparison.json';assert not target.exists()
    validate_source(BASE/'actual',config(BASE/'actual/case.cfg'),read(BASE/'actual/sst_initialization.json'))
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s']
    V=s.V;F=s.F.astype(int)-1;area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
    source=ROOT/'outputs/v127_validation/control/restart_flow.csv'
    a,cols,_,_,_,_=state(source,True,True);pocket=floor_nodes(a[:,cols.index('Omega')]);assert len(pocket)==146;del a
    variants={};heat={};pressure={};fluidT={};locks={str(source):sha(source)}
    residuals=['rms[Rho]','rms[RhoE]','rms[k]','rms[w]']
    for name,folder in [('reference',REFERENCE),('bounded',BASE/'actual')]:
        r=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert audit['interfacePassed'] and audit['checkpointLocalIteration']==99
        if name=='bounded':assert r['physicsModel']['numericalStudy']['studyID']==STUDY
        with (folder/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
        h=dict(zip([v.strip() for v in rows[0]],np.array(rows[1:],float).T));assert np.array_equal(h['Inner_Iter'],np.arange(100))
        snapshots={}
        for step in [50,99]:
            p=folder/('restart_flow_%06d.csv'%step);assert sha(p)==r['outputSHA256'][p.name]
            a,newcols,_,_,_,_=state(p,True,True);assert cols==newcols
            w=a[:,cols.index('Omega')];floor=floor_nodes(w)
            snapshots[str(step)]=dict(recoveredOriginalNodesAbove100Floor=int(sum(w[pocket]>.01)),
                originalFloorNodesRemaining=np.intersect1d(floor,pocket).tolist(),newFloorNodeIDs=np.setdiff1d(floor,pocket).tolist(),
                floorEquivalentNodeCount=len(floor),originalOmegaRange=[float(min(w[pocket])),float(max(w[pocket]))])
            locks[str(p)]=sha(p);del a
        xy,field=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xy);assert max(d)<1e-6 and len(set(ids))==len(V)
        q=np.empty(len(V));q[ids]=field['Heat_Flux'][:,0];heat[name]=q[F].mean(axis=1)
        p=np.empty(len(V));p[ids]=field['Pressure'][:,0];pressure[name]=p[F].mean(axis=1)
        t=np.empty(len(V));t[ids]=field['Temperature'][:,0];fluidT[name]=t[F].mean(axis=1)
        variants[name]=dict(elapsedSeconds=r['elapsedSeconds'],snapshots=snapshots,
            last20MeanLogResidual={k:float(np.mean(h[k][-20:])) for k in residuals},
            finalLogResidual={k:float(h[k][-1]) for k in residuals},
            history={k:h[k].tolist() for k in residuals+['HF','Min CFL','Max CFL']},
            nativeHeatPowerW=audit['nativeHeatPowerW'],originalFaceHeatPowerW=float(area@heat[name]),
            withinRunHeatRelativeL1Change=audit['withinRunHeatRelativeL1Change'],
            yPlusP99=audit['yPlusPercentiles']['p99'],originalResidualPassed=audit['originalResidualThresholdPassed'])
        for n in ['case.cfg','sst_initialization.json','engine_process.json','sst_state_audit.json','history.csv','surface_flow.vtk']:locks[str(folder/n)]=sha(folder/n)
    a=variants['reference'];b=variants['bounded'];gain={k:a['last20MeanLogResidual'][k]-b['last20MeanLogResidual'][k] for k in residuals}
    heatL1=float(area@abs(heat['bounded']-heat['reference'])/(area@abs(heat['reference'])))
    gates=dict(omegaMeanImprovedOneLog=gain['rms[w]']>=1,densityEnergyKNotWorse03=min(gain[k] for k in residuals[:3])>=-.3,
        all146RecoveredAtBothSnapshots=all(v['recoveredOriginalNodesAbove100Floor']==146 for v in b['snapshots'].values()),
        noNewFloorPositionsAtEitherSnapshot=all(not v['newFloorNodeIDs'] for v in b['snapshots'].values()),
        pairedHeatL1Below5Percent=heatL1<=.05,yPlusP99Below1=b['yPlusP99']<=1)
    report=dict(variants=variants,meanResidualGainLog=gain,screeningGates=gates,allScreeningGatesPassed=all(gates.values()),
        pairedFaceHeatRelativeL1Difference=heatL1,originalFaces=len(F),originalAreaM2=float(sum(area)),
        pairedFacePressureRelativeL1Difference=float(area@abs(pressure['bounded']-pressure['reference'])/(area@abs(pressure['reference']))),
        pairedFaceFluidNodeTemperatureRelativeL1Difference=float(area@abs(fluidT['bounded']-fluidT['reference'])/(area@abs(fluidT['reference']))),
        sourceLocalIteration=149,referenceRunLocalIteration=99,newLocalIteration=99,
        physicsQualified=False,newQualifiedArea=0,defaultChanged=False,automaticContinuationAllowed=False,
        sourceSHA256=dict(locks,**{str(Path(__file__).resolve()):sha(Path(__file__))}))
    write(target,report);print(json.dumps(dict(gates=gates,gain=gain,heatL1=heatL1,snapshots=b['snapshots']),indent=2))
if __name__=='__main__':main()
