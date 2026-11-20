"""对同一原生时间层核对六个残差、空间位置和完整场，保存可追溯的二次启动种子。"""
import argparse
import csv
import os
from pathlib import Path
import sys

import numpy as np
from run_relaxation import (ASSETS, ENGINE, FULL, HERE, OLD, PLAN, ROOT, RUN, frame, frozen,
                            mode, observe, read, save, sha, steady_config, terminal, write_seed)
from matched_step import RMS


def columns(path,names):
    with path.open(encoding="ascii") as stream:
        header=[key.strip() for key in next(csv.reader(stream,skipinitialspace=True))]
        return np.loadtxt(stream,delimiter=",",usecols=[header.index(name) for name in names],ndmin=2)


def rhs_at_endpoint(path,xyz,control):
    names=["Residual_"+v for v in ("Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0")]
    values=columns(path,names)
    assert values.shape==(len(xyz),6) and np.isfinite(values).all()
    logs=np.log10(np.sqrt(np.mean(values**2,axis=0)))
    h=values[:,4];squared=h*h;total=squared.sum()
    axis=np.hypot(xyz[:,1],xyz[:,2])==0
    downstream=xyz[:,0]<-0.5*1.179248562795375
    peak=int(np.argmax(squared))
    return dict(logRMS=dict(zip(RMS,logs)),maxLogDifferenceFromHistory=float(np.max(np.abs(logs-control))),
                peakHNode=peak,peakHPositionOverD=xyz[peak].tolist(),axisShare=float(squared[axis].sum()/total),
                downstreamShare=float(squared[downstream].sum()/total),
                downstreamAxisShare=float(squared[axis&downstream].sum()/total))


def inspect(folder,steady):
    rec=terminal(folder)
    count=2 if steady else 1
    assert rec["command"]==[str(ENGINE),"case.cfg","--threads",str(count)]
    hist=folder/("history.csv" if steady else "history_00001.csv")
    h=columns(hist,["Time_Iter","Inner_Iter"]+RMS)
    assert len(h)>=20 and np.isfinite(h).all()
    assert np.all(h[:,0]==(0 if steady else 1))
    np.testing.assert_array_equal(h[:,1],np.arange(len(h)))
    last_idx=len(h)-1
    end=folder/(f"restart_{last_idx:06d}.csv" if steady else "restart_00001.csv")
    log=(folder/"debugger.log").read_text(encoding="utf8",errors="replace")
    assert end.name in log and "Exit Success (SU2_CFD)" in log
    if steady:
        assert (folder/"case.cfg").read_text("ascii")==steady_config() and sha(folder/"restart.csv")==sha(end)
    else:
        expected=(FULL/"case.cfg").read_text("ascii").replace(
            "VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)","VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL)")
        assert (folder/"case.cfg").read_text("ascii")==expected
    with np.load(folder/"nodes.npz") as z:
        xyz,ids=z["xyz"],z["ids"]
    before=frame(folder/("initial_flow.csv" if steady else "seed_00000.csv"),xyz)
    state=frame(end,xyz)
    assert np.any(before[:,4:]!=state[:,4:])
    geom=read(FULL/"geometry.json")
    p=geom["airPa"]+geom["rho"]*geom["speedMps"]**2*state[:,4]
    t=state[:,9];wall=ids[:-1,:65]
    wall_err=float(np.max(abs(t[wall]-1)))
    slip=float(np.max(abs(state[wall,5:8])))
    win=h[-20:,2:].max(axis=0)
    snapshot=rhs_at_endpoint(end,xyz,h[-1,2:])
    passed=bool(win.max()<-8 and p.min()>0 and t.min()>=-1e-12 and t.max()<=1+1e-12 and
                wall_err<=1e-12 and slip<=1e-12 and snapshot["maxLogDifferenceFromHistory"]<1e-7)
    result=dict(normalNativeTermination=True,lastSavedIndex=last_idx if steady else 1,
                innerRows=len(h),firstLogRMS=dict(zip(RMS,h[0,2:])),lastLogRMS=dict(zip(RMS,h[-1,2:])),
                last20MaxLogRMS=dict(zip(RMS,win)),nodeRHS=snapshot,
                pressureRangePa=[float(p.min()),float(p.max())],thetaRange=[float(t.min()),float(t.max())],
                maxWallScalarError=wall_err,maxWallVelocityOverU=slip,passesWindowAndState=passed,
                initialM1=mode(before[:,5],ids),finalM1=mode(state[:,5],ids),
                maxAxialChangeOverU=float(np.max(abs(state[:,5]-before[:,5]))),frameName=end.name,
                frameSHA256=sha(end),historySHA256=sha(hist),elapsedSeconds=rec["elapsedSeconds"],
                peakRSSBytes=rec["peakFamilyRSSBytes"],outputBytes=rec["peakOutputBytes"])
    return state,ids,result


def prepare_step(folder,out):
    state,ids,result=inspect(folder,True)
    ctrl=read(RUN/"control01/control.json")
    indexes=("rms[P]","rms[h]")
    gain={k:10**(ctrl["nativeLogRMS"][k]-result["last20MaxLogRMS"][k]) for k in indexes}
    accepted=bool(result["passesWindowAndState"] and min(gain.values())>=10)
    # 各周期帧分别与自己的迭代行配对，保留下游轴线误差迁移的可核对记录。
    hist=columns(folder/"history.csv",["Inner_Iter"]+RMS)
    xyz=state[:,1:4]
    checkpoints={}
    for path in sorted(folder.glob("restart_*.csv")):
        idx=int(path.stem.split("_")[-1])
        assert hist[idx,0]==idx
        rhs=result["nodeRHS"] if path.name==result["frameName"] else rhs_at_endpoint(path,xyz,hist[idx,1:])
        assert rhs["maxLogDifferenceFromHistory"]<1e-7
        checkpoints[path.name]=dict(sha256=sha(path),nodeRHS=rhs)
    result.update(acceptedForInitialization=accepted,improvementVsCalm=gain,periodicFrames=checkpoints)
    save(out/"relaxation.json",result)
    if not accepted:
        return
    meta=read(RUN/"calm01/calm.json")
    delta=np.load(frozen(RUN/"calm01/impulse.npy",meta["impulseSHA256"]))
    old=state[:,5].copy()
    state[:,5]+=delta
    np.testing.assert_allclose(state[:,5]-old,delta,atol=2e-16,rtol=0)
    assert not np.any(delta[ids[:,:65]]) and mode(delta,ids)==meta["impulseM1"]
    write_seed(out/"seed_00000.csv",state)
    save(out/"relaxed_seed.json",dict(acceptedForInitialization=True,seedSHA256=sha(out/"seed_00000.csv"),
         relaxationFrameSHA256=result["frameSHA256"],impulseSHA256=meta["impulseSHA256"],
         baseM1=result["finalM1"],addedM1=mode(delta,ids),totalM1=mode(state[:,5],ids)))


def check_step(folder,out):
    _,_,result=inspect(folder,False)
    result["acceptedFirstStep"]=result["passesWindowAndState"]
    result["statisticalWakeEstablished"]=False
    result["newComputedFaces"]=0
    result["newComputedAreaM2"]=0
    save(out/"first_step.json",result)


def supervised(stage,folder,out):
    out.mkdir(parents=True,exist_ok=False)
    files=[PLAN,HERE/"run_relaxation.py",Path(__file__).resolve(),OLD/"matched_step.py",
           folder/"debugger_process.json",folder/"case.cfg",folder/"nodes.npz",RUN/"control01/control.json"]
    files.extend(folder.glob("restart_*.csv"))
    files.extend(folder.glob("history*.csv"))
    files.extend([RUN/"calm01/calm.json",RUN/"calm01/impulse.npy"])
    for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[key]="1"
    result=observe([sys.executable,str(Path(__file__).resolve()),"_"+stage,str(folder),str(out)],
                   out,{str(p):sha(p) for p in files},90,128*1024**2)
    return result["debuggerExitCode"]==0 and result["stopReason"] is None


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=("ready","finish","_ready","_finish"))
    parser.add_argument("source",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    folder,out=args.source.resolve(),args.output.resolve()
    folder.relative_to(RUN);out.relative_to(RUN)
    if args.stage=="_ready":
        prepare_step(folder,out)
    elif args.stage=="_finish":
        check_step(folder,out)
    else:
        raise SystemExit(0 if supervised(args.stage,folder,out) else 1)
