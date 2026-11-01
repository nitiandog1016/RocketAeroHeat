"""复用v145两个具名状态及独立初猜定位k下限；不追加求解或改下限。"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.io import loadmat
from prepareSU2CheckpointContinuation import ROOT,read,sha

def main():
    base=ROOT/'outputs/v145_validation';target=base/'startup_diagnostic.json';assert not target.exists()
    case=Path(read(base/'main01/engine_folder.json')['folder']);receipt=read(case/'engine_process.json')
    d=loadmat(str(base/'mesh/hybrid_volume.mat'));P=d['vertices'];nv=len(d['sourceVertices'])
    prepared=d['sourcePreparedNodeIDsOneBased'].ravel().astype(np.int64)-1
    original=(prepared>=0)&(prepared<33*nv);layers=np.full(len(P),-1,np.int64);layers[original]=prepared[original]//nv
    locked={str(case/'engine_process.json'):sha(case/'engine_process.json')};cases={}
    init=case/'initial_flow.csv';initial=np.loadtxt(str(init),delimiter=',',skiprows=1,usecols=[9,10])
    cases['initial_guess']=dict(kFloorOffWallNodes=int(sum(initial[nv:,0]<=1e-10*(1+8*np.finfo(float).eps))),
        omegaFloorCount=int(sum(initial[:,1]<=1e-4*(1+8*np.finfo(float).eps))))
    locked[str(init)]=sha(init)
    for step in [75,149]:
        path=case/('restart_flow_%06d.csv'%step);assert sha(path)==receipt['outputSHA256'][path.name]
        with path.open() as f:cols=next(csv.reader(f))
        keys=['Turb_Kin_Energy','Omega','Temperature','Pressure']
        a=np.loadtxt(str(path),delimiter=',',skiprows=1,usecols=[cols.index(k) for k in keys]);assert len(a)==len(P)
        floor=np.flatnonzero(a[:,0]<=1e-10*(1+8*np.finfo(float).eps));off=floor[floor>=nv]
        indices,count=np.unique(layers[off],return_counts=True);j=int(np.argmax(a[:,2]))
        cases[str(step)]=dict(kFloorAllNodes=len(floor),kFloorOffWallNodes=len(off),
            offWallReferenceLayerCounts={str(int(i)):int(n) for i,n in zip(indices,count)},
            offWallKFloorXRangeM=[float(min(P[off,0])),float(max(P[off,0]))] if len(off) else None,
            omegaFloorCount=int(sum(a[:,1]<=1e-4*(1+8*np.finfo(float).eps))),
            maxTemperatureK=float(a[j,2]),maxTemperatureNodeZeroBased=j,maxTemperatureXYZM=P[j].tolist(),
            pressureRangePa=[float(min(a[:,3])),float(max(a[:,3]))])
        locked[str(path)]=sha(path)
    report=dict(cases=cases,referenceLayerMinus1MeansFanOrNewExterior=True,
        limitations='reference layer identities only; no native local RHS/update saved; does not isolate turbulence production, diffusion or convection causality',
        solverRuns=0,physicsQualified=False,sourceSHA256=dict(locked,**{str(Path(__file__).resolve()):sha(Path(__file__).resolve())}))
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(cases,indent=2))
if __name__=='__main__':main()
