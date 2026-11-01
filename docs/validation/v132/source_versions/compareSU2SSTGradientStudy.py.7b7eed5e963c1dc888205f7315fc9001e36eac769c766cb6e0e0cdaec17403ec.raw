"""按v132事前门槛比较原生梯度：保留全部原面、下限和源/本段身份。"""
import csv,json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from prepareSU2SSTGradientStudy import ROOT,BASE,SOURCE,PLAN,config,read,sha,write,validate_source
from prepareSU2SSTPocketStudy import floor_nodes
from su2SSTInterface import state
from auditSU2SSTState import vtk_tokens

def main():
    target=BASE/'gradient_comparison.json';assert not target.exists()
    a,cols,_,_,_,_=state(SOURCE/'restart_flow.csv',True,True)
    original=floor_nodes(a[:,cols.index('Omega')]);del a
    s=loadmat(str(ROOT/'outputs/v118_validation/actual01/surface.mat'),squeeze_me=True,struct_as_record=False)['s']
    V=s.V;F=s.F.astype(int)-1
    area=np.linalg.norm(np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]]),axis=1)/2
    assert len(F)==19752 and abs(sum(area)-1.9646292565715553)<1e-12
    residuals=['rms[Rho]','rms[RhoE]','rms[k]','rms[w]'];variants={};faces={};locks={}
    for name in ('wls','gg'):
        folder=BASE/name;validate_source(folder,config(folder/'case.cfg'),read(folder/'sst_initialization.json'))
        r=read(folder/'engine_process.json');audit=read(folder/'sst_state_audit.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert audit['interfacePassed'] and audit['checkpointLocalIteration']==9
        with (folder/'history.csv').open() as f:rows=list(csv.reader(f,skipinitialspace=True))
        h=dict(zip([v.strip() for v in rows[0]],np.array(rows[1:],float).T));assert np.array_equal(h['Inner_Iter'],np.arange(10))
        a,keys,_,_,_,_=state(folder/'restart_flow.csv',True,True);assert keys==cols
        w=a[:,cols.index('Omega')];k=a[:,cols.index('Turb_Kin_Energy')];floors=floor_nodes(w)
        snapshot=dict(recoveredSourceFloorNodesAbove100Floor=int(sum(w[original]>.01)),
            sourceFloorNodesRemaining=np.intersect1d(floors,original).tolist(),
            newFloorNodeIDs=np.setdiff1d(floors,original).tolist(),floorEquivalentNodeCount=len(floors),
            sourceFloorOmegaValues=w[original].tolist(),minimumK=float(min(k)),minimumOmega=float(min(w)))
        xyz,native=vtk_tokens(folder/'surface_flow.vtk');d,ids=cKDTree(V).query(xyz)
        assert max(d)<1e-6 and len(set(ids))==len(V)
        fields={}
        for field in ('Pressure','Heat_Flux','Temperature'):
            n=np.empty(len(V));n[ids]=native[field][:,0];fields[field]=n[F].mean(axis=1)
        faces[name]=fields
        variants[name]=dict(elapsedSeconds=r['elapsedSeconds'],snapshot=snapshot,
            last5MeanLogResidual={key:float(np.mean(h[key][-5:])) for key in residuals},
            finalLogResidual={key:float(h[key][-1]) for key in residuals},
            history={key:h[key].tolist() for key in residuals+['HF']},
            originalFaceHeatPowerW=float(area@fields['Heat_Flux']),nativeHeatPowerW=audit['nativeHeatPowerW'],
            originalResidualPassed=audit['originalResidualThresholdPassed'],interfacePassed=audit['interfacePassed'])
        for file in ['engine_process.json','sst_state_audit.json','history.csv','restart_flow.csv','surface_flow.vtk']:
            locks[str(folder/file)]=sha(folder/file)
    gain={key:variants['wls']['last5MeanLogResidual'][key]-variants['gg']['last5MeanLogResidual'][key] for key in residuals}
    differences={field:float(area@abs(faces['gg'][field]-faces['wls'][field])/(area@abs(faces['wls'][field]))) for field in faces['wls']}
    b=variants['gg']['snapshot']
    gates=dict(allSourceFloorRecovered=b['recoveredSourceFloorNodesAbove100Floor']==len(original),
        noNewFloorPositions=not b['newFloorNodeIDs'],omegaMeanImprovedOneLog=gain['rms[w]']>=1,
        densityEnergyKNotWorse02=min(gain[key] for key in residuals[:3])>=-.2,
        pairedHeatL1Below5Percent=differences['Heat_Flux']<=.05,
        pairedPressureL1Below5Percent=differences['Pressure']<=.05,
        legalStateAndInterface=all(v['interfacePassed'] and v['snapshot']['minimumK']>=1e-10 for v in variants.values()))
    result=dict(variants=variants,sourceFloorNodeIDs=original.tolist(),sourceLocalIteration=99,newLocalIteration=9,
        meanResidualGainLog=gain,screeningGates=gates,allScreeningGatesPassed=all(gates.values()),
        pairedFaceRelativeL1Differences=differences,originalFaces=len(F),originalAreaM2=float(sum(area)),
        physicsQualified=False,newQualifiedArea=0,defaultChanged=False,automaticContinuationAllowed=False,
        sourceSHA256=dict(locks,**{str(p):sha(p) for p in [Path(__file__).resolve(),PLAN,SOURCE/'restart_flow.csv']}))
    write(target,result);print(json.dumps(dict(gates=gates,gain=gain,differences=differences,snapshot=b),indent=2))

if __name__=='__main__':main()
