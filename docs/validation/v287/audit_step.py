"""读取全周首步的原生进程、逐次残差和保存帧，独立判定数值门槛。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha

RMS=["rms["+name+"]" for name in ("P","U","V","W","h","rho*Y_0")]


def frame(file, xyz):
    """核对 SU2 的列、原始节点号和坐标，不能把原初场作为新保存态。"""
    with file.open(encoding="ascii") as stream:
        header=[s.strip() for s in next(csv.reader(stream,skipinitialspace=True))]
    required=["PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"]
    cols=[header.index(key) for key in required]
    a=np.loadtxt(file,delimiter=",",skiprows=1,usecols=cols)
    assert a.shape==(len(xyz),len(cols)) and np.isfinite(a).all()
    np.testing.assert_array_equal(a[:,0],np.arange(len(a)))
    np.testing.assert_allclose(a[:,1:4],xyz,atol=2e-14,rtol=1e-14)
    return a


def audit(folder):
    rec=json.loads((folder/"debugger_process.json").read_bytes())
    assert rec["status"]=="terminal" and not rec["survivingProcessIDs"]
    assert all(rec["inputsStable"].values()) and rec["supervisorError"] is None
    for path,digest in rec["inputSHA256"].items():
        assert sha(Path(path))==digest,path
    log=(folder/"debugger.log").read_text(encoding="utf8",errors="replace")
    normal=rec["debuggerExitCode"]==0 and rec["stopReason"] is None and "Exit Success (SU2_CFD)" in log
    geom=json.loads((folder/"geometry.json").read_bytes())
    init=folder/"seed_00000.csv"
    out=dict(nativeExitCode=rec["debuggerExitCode"],normalTerminal=normal,stopReason=rec["stopReason"],
             nativeSeconds=rec["elapsedSeconds"],peakProcessTreeRSSBytes=rec["peakFamilyRSSBytes"],
             newOutputBytes=rec["peakOutputBytes"],originalFaces=5216,newComputedFaces=0,newComputedAreaM2=0,
             acceptedFirstStep=False,statisticalWakeEstablished=False,sha256={init.name:sha(init)})
    hist=folder/"history_00001.csv"
    if hist.exists():
        with hist.open(encoding="ascii") as f:
            keys=[s.strip() for s in next(csv.reader(f,skipinitialspace=True))]
        h=np.loadtxt(hist,delimiter=",",skiprows=1,ndmin=2)
        assert np.isfinite(h).all() and np.all(h[:,keys.index("Time_Iter")]==1)
        np.testing.assert_array_equal(h[:,keys.index("Inner_Iter")],np.arange(len(h)))
        cols=[keys.index(name) for name in RMS]
        win=h[-20:,cols].max(axis=0)
        out.update(innerRows=len(h),firstLogRMS=dict(zip(RMS,h[0,cols])),
                   finalLogRMS=dict(zip(RMS,h[-1,cols])),last20MaxLogRMS=dict(zip(RMS,win)),
                   residualWindowPassed=bool(len(h)>=20 and np.max(win)<-8))
        out["sha256"][hist.name]=sha(hist)
    endpoint=folder/"restart_00001.csv"
    out["savedStepOneFrame"]=endpoint.is_file()
    if endpoint.is_file():
        with np.load(folder/"nodes.npz") as grid:
            ids,xyz=grid["ids"],grid["xyz"]
        start=frame(init,xyz);end=frame(endpoint,xyz)
        assert np.array_equal(ids[0],ids[-1]) and len(ids)==geom["angularCells"]+1
        assert np.any(start[:,4:]!=end[:,4:])
        wall=ids[:-1,:65].ravel()
        def mode(a):
            # 对每个 (x,r) 圆环取轴向速度的一阶谐波幅值，全周有独立自由度。
            harm=np.fft.rfft(a[ids[:-1],5],axis=0)
            return float(2*np.max(np.abs(harm[1]))/geom["angularCells"])
        init_mode=mode(start)
        assert abs(init_mode-geom["peakAxialPerturbationOverU"])<1e-12
        p=geom["airPa"]+geom["rho"]*geom["speedMps"]**2*end[:,4]
        t=end[:,9];bounded=bool(t.min()>=-1e-12 and t.max()<=1+1e-12)
        wall_error=float(np.max(np.abs(t[wall]-1)))
        no_slip=float(np.max(np.abs(end[wall,5:8])))
        out.update(thetaRange=[float(t.min()),float(t.max())],positiveAbsolutePressure=bool(np.all(p>0)),
                   pressureRangePa=[float(p.min()),float(p.max())],strongWallError=wall_error,
                   wallVelocityMaxOverU=no_slip,scalarBoundsPassed=bounded,initialM1=init_mode,
                   finalM1=mode(end),maxDimensionlessAxialVelocityChange=float(np.max(np.abs(end[:,5]-start[:,5]))))
        out["acceptedFirstStep"]=bool(normal and out.get("residualWindowPassed") and bounded and
                                      p.min()>0 and wall_error<=1e-12 and no_slip<=1e-12)
        out["sha256"][endpoint.name]=sha(endpoint)
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    parser.add_argument("result",type=Path)
    args=parser.parse_args()
    result=audit(args.folder.resolve())
    with args.result.open("x",encoding="utf8") as f:
        json.dump(result,f,indent=2)
        f.write("\n")
    print(json.dumps(result))
