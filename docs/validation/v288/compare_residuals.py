"""将两份局部方程右端与同次 SU2 历史配对，区分原生残差、节点相关及全周首步成本。"""
import argparse
import csv
from pathlib import Path

import numpy as np
from matched_step import FULL, PREVIOUS, RMS, frame, read, save, sha, terminal


def one_step(folder):
    rec=terminal(folder)
    endpoint=folder/"restart_00001.csv"
    log=(folder/"debugger.log").read_text(encoding="utf8",errors="replace")
    assert "Exit Success (SU2_CFD)" in log and endpoint.name in log
    with (folder/"history_00001.csv").open(encoding="ascii") as f:
        keys=[k.strip() for k in next(csv.reader(f,skipinitialspace=True))]
        history=np.loadtxt(f,delimiter=",",ndmin=2)
    assert history.shape==(1,len(keys)) and history[0,keys.index("Inner_Iter")]==0 and history[0,keys.index("Time_Iter")]==1
    native={k:float(history[0,keys.index(k)]) for k in RMS}
    with np.load(folder/"nodes.npz") as z:
        xyz,ids=z["xyz"],z["ids"]
    seed=frame(folder/"seed_00000.csv",xyz)
    final=frame(endpoint,xyz)
    assert np.any(seed[:,4:]!=final[:,4:])
    assert (101325+1.18*36*final[:,4]).min()>0 and final[:,9].min()>=-1e-12 and final[:,9].max()<=1+1e-12
    assert np.all(final[ids[:-1,:65],9]==1)
    columns=["Residual_"+k for k in ["Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"]]
    with endpoint.open(encoding="ascii") as f:
        keys=[v.strip() for v in next(csv.reader(f,skipinitialspace=True))]
        data=np.loadtxt(f,delimiter=",",usecols=[keys.index(k) for k in columns])
    assert len(data)==len(xyz) and np.isfinite(data).all()
    computed=np.log10(np.sqrt(np.mean(data**2,axis=0)))
    pairs=dict(zip(RMS,computed))
    geom=read(folder/"geometry.json")
    left,right=data[:,0],data[:,4]
    peak=int(np.argmax(np.abs(right)))
    # 采用该次输入态的局部比焓，检查 R_h 与 h_i R_p 数值接近程度；只对账，不替换原残差。
    h=seed[:,8]
    projected=h*left
    fraction=float(np.linalg.norm(right-projected)/np.linalg.norm(right))
    levels=np.empty(len(xyz),int)
    levels[ids[:-1]]=np.arange(ids.shape[1])//65
    regions=dict(wall=levels==0,nextTwoLayers=(levels>=1)&(levels<=2),
                 outerBoundary=levels==48,axis=np.hypot(xyz[:,1],xyz[:,2])==0,
                 downstream=xyz[:,0]<-0.5*1.179248562795375)
    norm2=np.dot(right,right)
    return dict(seedSHA256=sha(folder/"seed_00000.csv"),frameSHA256=sha(endpoint),
                fullMeshSHA256=sha(folder/"mesh.su2"),initialMode=float(
                    (2*np.abs(np.fft.rfft(seed[ids[:-1],5],axis=0)[1])/geom["angularCells"]).max()),
                nativeSeconds=rec["elapsedSeconds"],nativeLogRMS=native,nodalLogRMS=pairs,
                pressureEnthalpyLogMatch=all(abs(native[k]-pairs[k])<1e-7 for k in (RMS[0],RMS[4])),
                relativeNormAfterSubtractingLocalHtimesMass=fraction,
                cosineOfEnthalpyAndLocalHtimesMass=float(np.dot(right,projected)/(np.linalg.norm(right)*np.linalg.norm(projected))),
                leastSquaresEnthalpyPerMass=float(np.dot(left,right)/np.dot(left,left)),
                seedNativeEnthalpyRange=[float(h.min()),float(h.max())],
                maxAbsEnthalpyRhsNode=peak,maxNodeXYZoverD=xyz[peak].tolist(),
                fractionsOfSquaredEnthalpyRhs={k:float(np.dot(right[m],right[m])/norm2) for k,m in regions.items()})


def compare(root):
    ray=one_step(root/"probe_ray01")
    mirror=one_step(root/"probe_mirror01")
    new=read(root/"audit01/result.json")
    old=read(PREVIOUS/"evidence/first_step.json")
    assert ray["nativeLogRMS"]==new["firstLogRMS"] and ray["seedSHA256"]==new["sha256"]["seed_00000.csv"]
    assert ray["fullMeshSHA256"]==mirror["fullMeshSHA256"]==sha(FULL/"mesh.su2")
    assert abs(ray["initialMode"]-mirror["initialMode"])<1e-12
    change={k:10**(ray["nativeLogRMS"][k]-mirror["nativeLogRMS"][k]) for k in RMS}
    win={k:10**(old["last20MaxLogRMS"][k]-new["last20MaxLogRMS"][k]) for k in RMS}
    return dict(ray=ray,mirror=mirror,rayOverMirrorFirstRMS=change,requiredReduction=10,
                runFull=bool(change[RMS[0]]>=10 and change[RMS[4]]>=10),
                oldOverMatchedWindowRMS=win,originalStepAccepted=old["acceptedFirstStep"],
                matchedStepAccepted=new["acceptedFirstStep"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("result",type=Path)
    args=parser.parse_args()
    save(args.result,compare(args.root.resolve()))
