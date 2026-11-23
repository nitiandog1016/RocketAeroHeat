"""固定真实物理端点读取银球壁反力、七个测站以及含 BDF1 储能项的完整收支。"""
import argparse
import csv
from pathlib import Path

import numpy as np

from shortSphereWake import DT, SAMPLES, CASE, ROOT, build_observer, native, observer_dir, read, sha, write_json
from shortSphere65 import conductivity
from thermistorStartupSurface import columns


def saved_evaluation(folder,step,source,limit,mode):
    seal=read(folder/"native.json");receipt=read(folder/"debugger_process.json")
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
    assert not(receipt["stopReason"] or receipt["survivingProcessIDs"] or receipt["supervisorError"])
    path=folder/f"restart_{step:05d}.csv";history=folder/f"history_{step:05d}.csv"
    assert seal["step"]==step and seal["mode"]==mode and seal["old"]==source==sha(folder/f"seed_{step-1:05d}.csv")
    assert seal["frame"]==sha(path) and seal["history"]==sha(history) and seal["inner"]==limit
    return path,history


def readout(work,branch,step):
    dt=DT[branch];j=read(work/"job.json");s=j["condition"];prep=Path(j["prepared"])
    physical=work/branch/f"step{step:03d}";r=read(physical/"native.json")
    assert max(r["window"].values())<-8 and r["inner"]>=20
    start=physical/f"seed_{step-1:05d}.csv";frame=physical/f"restart_{step:05d}.csv"
    assert sha(start)==r["old"] and sha(frame)==r["frame"]
    area_dir=physical/"readout"
    if (area_dir/"accepted.json").exists():
        seal=read(area_dir/"accepted.json")
        assert seal["frame"]==r["frame"] and seal["old"]==r["old"]
        for name,h in seal["files"].items():assert sha(area_dir/name)==h
        saved_evaluation(area_dir/seal.get("observerDirectory","observer"),step,r["old"],r["inner"]+1,"observe")
        saved_evaluation(area_dir/"release_native",step,r["frame"],1,"release")
        return seal
    area_dir.mkdir(exist_ok=True);release=area_dir/"release_native";n=r["inner"]
    raw=columns(physical/f"history_{step:05d}.csv",list(r["window"]))
    for workers in (1,4):
        build=observer_dir(work,workers);observer=area_dir/build.name
        if observer.exists():
            path,history=saved_evaluation(observer,step,r["old"],n+1,"observe")
            observed=np.loadtxt(path,delimiter=",",skiprows=1,usecols=range(10));hist=columns(history,list(r["window"]))
        else:
            if not (build/"build.json").exists():build_observer(work,workers)
            observed,hist,_=native(work,branch,step,start,observer,n+1,"observe",workers)
        history_error=float(abs(raw-hist[:-1]).max())
        if history_error<.01:break
        reject=observer/"prefix_rejected.json"
        if not reject.exists():
            write_json(reject,dict(maxLogDifference=history_error,limit=.01,workers=workers,frame=r["frame"]))
        if workers==4:raise ValueError(f"同线程回放历史差 {history_error:.8g} 超过门槛")
    new=np.loadtxt(frame,skiprows=1,delimiter=",",usecols=range(10))
    prior=np.loadtxt(start,skiprows=1,delimiter=",",usecols=range(10))
    events=list(csv.reader((observer/"readout_updates.csv").read_text("ascii").splitlines()))
    assert events==[[str(step),str(n),"flow"],[str(step),str(n),"scalar"]]
    rows=np.loadtxt(observer/"before_wall.csv",delimiter=",",skiprows=1);rows=rows[rows[:,0].argsort()]
    after=np.loadtxt(observer/"after_wall.csv",delimiter=",",skiprows=1);after=after[after[:,0].argsort()]
    np.testing.assert_array_equal(rows[:,0],new[:,0]);np.testing.assert_array_equal(after[:,0],new[:,0])
    differences=dict(flow=float(abs(rows[:,5:9]-new[:,4:8]).max()),scalar=float(abs(rows[:,2]-new[:,9]).max()),
                     previous=float(abs(rows[:,3]-prior[:,9]).max()),frozenFrame=float(abs(observed-new).max()))
    assert differences["flow"]<1e-8 and differences["scalar"]<2e-10 and differences["previous"]<1e-13 and differences["frozenFrame"]<1e-8,differences
    assert np.all(rows[:,4]>0) and np.all(rows[:,9]==1) and np.all(rows[:,10]==dt)
    forces=rows[:,1]-after[:,1]
    with np.load(prep/"mapping.npz") as z:
        wet,first,f,centers,area,lump,line,angle=(z[k] for k in ("wet","first","F","centers","area","lump","line","angles"))
    assert np.count_nonzero(forces)==np.count_nonzero(forces[wet]) and np.all(after[wet,1]==0)
    np.testing.assert_allclose(rows[wet,2],rows[wet,3],rtol=0,atol=2e-14)
    # 解除壁的强标量约束，用同一端点重新装配。壁温不随时间变化，壁节点的 BDF1 项为零。
    if release.exists():saved_evaluation(release,step,r["frame"],1,"release")
    else:native(work,branch,step,frame,release,1,"release")
    check=columns(release/f"restart_{step:05d}.csv",["PointID","Residual_Species_0"])
    np.testing.assert_array_equal(check[:,0],new[:,0])
    release_diff=float(np.max(abs(check[wet,1]+forces[wet]))/np.max(abs(forces[wet])))
    assert release_diff<1e-5,release_diff
    scale=s["rho"]*s["cp"]*s["speedMps"]*s["diameterM"]**2*s["scalarUnitK"]
    sums={i:{} for i in (n-1,n)}
    for i,k,val in csv.reader((observer/"stages.csv").read_text("ascii").splitlines()):
        assert k not in sums[int(i)];sums[int(i)][k]=float(val)*scale
    keys=["start","beforeTime","afterTime","beforeFar","afterFar","beforeWall","afterWall"]
    assert all(list(d)==keys for d in sums.values());ledger=sums[n]
    power=ledger["beforeWall"]-ledger["afterWall"];far=ledger["afterFar"]-ledger["beforeFar"]
    storage=ledger["afterTime"]-ledger["beforeTime"];independent=rows[:,4]@(new[:,9]-prior[:,9])*scale/dt
    previous=sums[n-1];change=abs((previous["beforeWall"]-previous["afterWall"])/power-1)
    balance=(power-far-storage)/power
    assert power>0 and abs(independent-storage)<1e-9 and abs(balance)<.005 and change<.002
    assert max(abs(ledger["start"]),abs(ledger["beforeTime"]),abs(ledger["beforeFar"]-ledger["afterTime"]),
               abs(ledger["beforeWall"]-ledger["afterFar"]))<1e-10
    node_q=-forces[wet]*scale/lump;assert abs(node_q@lump+power)<1e-11
    # 三角面与 v321 逐号一致。首层、强温界与压力均取同一接受帧。
    def faces(x):return np.r_[x,x[centers].mean(axis=1)][f].mean(axis=1)
    gauge=s["rho"]*s["speedMps"]**2*new[wet,4];wall=s["airK"]+s["scalarUnitK"]*new[wet,9]
    q=faces(node_q);assert abs(q@area+power)<1e-11
    near=faces(s["airK"]+s["scalarUnitK"]*new[first,9]);p=faces(gauge)
    np.savez_compressed(area_dir/"fields.npz",q=q,gauge=p,first=near,wall=faces(wall),volume=rows[:,4])
    theta=np.asarray(read(CASE)["angleDeg"]);axis=angle[line]
    sample_wall=np.interp(theta,axis,wall[line]);sample_q=np.interp(theta,axis,node_q[line])
    nu=-sample_q*s["diameterM"]/(conductivity(sample_wall)*(sample_wall-s["airK"]))
    seal=dict(frame=r["frame"],old=r["old"],physicalStep=step,dt=dt,timeDU=step*dt,seconds=step*dt*s["diameterM"]/s["speedMps"],
              wallW=power,farW=far,storageW=storage,storageFromFramesW=independent,relativeBalance=balance,
              relativeWallChange=change,releaseError=release_diff,historyError=history_error,stateErrors=differences,
              localNu=nu.tolist(),stationFluxIn=sample_q.tolist(),stagesW=ledger,previousStagesW=previous,window=r["window"],
              observerDirectory=observer.name,
              files={p.relative_to(area_dir).as_posix():sha(p) for p in [area_dir/"fields.npz",release/"native.json",
                     *(observer/k for k in ("before_wall.csv","after_wall.csv","stages.csv","readout_updates.csv","native.json"))]})
    write_json(area_dir/"accepted.json",seal);print(branch,step,nu[-2:],power,balance,flush=True);return seal


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("work",type=Path);p.add_argument("branch",choices=DT)
    p.add_argument("--through",type=float,required=True);a=p.parse_args();w=a.work.resolve();w.relative_to(ROOT/"outputs")
    from shortSphereWake import advance
    assert a.through in SAMPLES
    for time in SAMPLES:
        if time>a.through:break
        advance(w,a.branch,time);readout(w,a.branch,int(round(time/DT[a.branch])))
