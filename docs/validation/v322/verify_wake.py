"""独立由物理时间层、自由壁残差、原生面几何与测站坐标复核银球尾迹交付。"""
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(HERE.parent/"v321"))
from verify import close, digest, read, wall_faces


def table(file):
    with file.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))


def state(file):
    return np.loadtxt(file,skiprows=1,delimiter=",",usecols=range(10))


def ring_rms(velocity,xyz,volume):
    radius=np.linalg.norm(xyz[:,1:],axis=1)
    use=(xyz[:,0]>=-3)&(xyz[:,0]<=-.6)&(radius<=1)&(radius>1e-10)
    pairs=np.round(np.c_[xyz[use,0],radius[use]],11)
    groups,ids=np.unique(pairs,axis=0,return_inverse=True)
    u=velocity[use];unit=xyz[use,1:]/radius[use,None]
    local=np.c_[u[:,0],(u[:,1:]*unit).sum(1),u[:,2]*unit[:,0]-u[:,1]*unit[:,1]]
    weights=volume[use];variance=0.
    for i in range(len(groups)):
        ix=ids==i;v=local[ix];w=weights[ix]
        variance+=np.sum(w[:,None]*(v-np.average(v,axis=0,weights=w))**2)
    return np.sqrt(variance/weights.sum())


def verify(work,result,ui,branch,mode):
    job=read(work/"job.json");prep=Path(job["prepared"]);s=job["condition"];p=read(prep/"prepared.json")
    summary=read(result/branch/"summary.json");comparison=read(result/"comparison.json")
    dt=summary["dtDU"];end=summary["physicalStep"];dt_s=dt*s["diameterM"]/s["speedMps"]
    proof=read(ui.parent/"ui.json")[mode];assert proof["error"]==0 and proof["pixels"]>0
    start=state(work/"seed.csv");base=Path(job["base"]);endpoint=read(base/"endpoints.json");steady=state(base/endpoint["stage"]/"restart.csv")
    diff=start-steady;close(diff[:,[0,1,2,3,4,6,8,9]],0,0)
    xyz=start[:,1:4];dx=xyz[:,0]+1.15;inside=(dx/.55)**2+(xyz[:,1]/.45)**2+(xyz[:,2]/.45)**2<1
    assert np.all(diff[~inside]==0);close(np.linalg.norm(diff[:,5:8],axis=1).max(),.01,2e-15)
    assert digest(work/"seed.csv")==comparison["seed"]==proof.get("seed",job["seedSHA256"])
    # 从 SU2 的原边界面重建分摊，与计算侧保存的三角权重独立。
    walls=wall_faces(prep/"mesh.su2");wet=np.unique([i for f in walls for i in f]);pts=xyz*s["diameterM"]
    dual=np.zeros(len(pts));centroids=[];face_nodes=[]
    for face in walls:
        tri=pts[face]
        if len(face)==3:
            a=np.linalg.norm(np.cross(tri[1]-tri[0],tri[2]-tri[0]))/2
            dual[face]+=a/3;centroids.append(tri.mean(0));face_nodes.append((face,[]))
        else:
            center=tri.mean(0);v=tri-center;a=np.linalg.norm(np.cross(v,np.roll(v,-1,axis=0)),axis=1)/2
            dual[face]+=(a+np.roll(a,1))/3+a.sum()/12
            for j in range(4):
                ids=[face[j],face[(j+1)%4]];centroids.append((pts[ids].sum(0)+center)/3);face_nodes.append((ids,face))
    assert np.all(diff[wet]==0)
    raw=(ui/"workbench/case-data.js").read_text("utf8").strip()
    page=json.loads(raw[len("window.ROCKET_CASE="):-1]);geometry=np.array(page["vertices"])[page["faces"]]
    area=np.linalg.norm(np.cross(geometry[:,1]-geometry[:,0],geometry[:,2]-geometry[:,0]),axis=1)/2
    dist,order=cKDTree(centroids).query(geometry.mean(1));assert max(dist)<1e-14 and len(np.unique(order))==4464
    dist,first=cKDTree(pts).query(pts[wet]*1.004);assert max(dist)<1e-14
    angle=np.degrees(np.arctan2(np.linalg.norm(pts[:,1:],axis=1),pts[:,0]))
    line=wet[(abs(pts[wet,1])<1e-15)&(pts[wet,2]>=-1e-15)];line=line[np.argsort(angle[line])]
    assert len(line)==49
    levels=read(ROOT/"cases/short_sphere/test65.json")["angleDeg"]
    old=start;last=job["seedSHA256"];max_window=-np.inf;checks=[];times=[r["timeDU"] for r in summary["history"]]
    total=0.;peak=0;trace=np.array([r["nonAxisymmetricRmsOverU"] for r in summary["wake"]]);assert len(trace)==end+1
    for step in range(1,end+1):
        folder=work/branch/f"step{step:03d}";rec=read(folder/"native.json");proc=read(folder/"debugger_process.json")
        assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
        assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
        assert last==rec["old"]==digest(folder/f"seed_{step-1:05d}.csv")
        frame=folder/f"restart_{step:05d}.csv";last=digest(frame);assert last==rec["frame"]
        hist=np.loadtxt(folder/f"history_{step:05d}.csv",delimiter=",",skiprows=1,ndmin=2)
        assert len(hist)==rec["inner"] and np.all(hist[:,0]==step) and np.all(hist[:,2]==np.arange(len(hist)))
        assert rec["step"]==step and rec["dt"]==dt and rec["history"]==digest(folder/f"history_{step:05d}.csv")
        close(hist[-20:,3:9].max(axis=0),list(rec["window"].values()),0)
        m=hist[-20:,3:9].max();assert m<-8;max_window=max(max_window,m);total+=rec["seconds"];peak=max(peak,rec["peakRSS"])
        new=state(frame);np.testing.assert_array_equal(new[:,0:4],old[:,0:4]);assert np.isfinite(new).all()
        assert min(s["airPa"]+s["rho"]*s["speedMps"]**2*new[:,4])>0
        close(new[wet,9],start[wet,9],2e-14);assert new[:,9].min()>=-1e-12 and new[:,9].max()<=start[wet,9].max()+1e-12
        if step*dt in times:
            d=folder/"readout";a=read(d/"accepted.json");assert a["frame"]==last and a["old"]==rec["old"]
            observer=a.get("observerDirectory","observer")
            for tag in (observer,"release_native"):
                proc=read(d/tag/"debugger_process.json");sub=read(d/tag/"native.json")
                assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
                assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
                assert sub["frame"]==digest(d/tag/frame.name)
            replay=np.loadtxt(d/observer/f"history_{step:05d}.csv",delimiter=",",skiprows=1,ndmin=2)
            err=abs(replay[:-1,3:9]-hist[:,3:9]).max();assert len(replay)==len(hist)+1 and err<.01 and err==a["historyError"]
            obs=np.loadtxt(d/observer/"before_wall.csv",delimiter=",",skiprows=1);obs=obs[obs[:,0].argsort()]
            np.testing.assert_array_equal(obs[:,0],new[:,0]);close(obs[:,3],old[:,9],1e-13);close(obs[:,2],new[:,9],2e-10);close(obs[:,5:9],new[:,4:8],1e-8)
            volume=obs[:,4];close(volume.sum(),p["geometry"]["volumeD3"],1e-7)
            scale=s["rho"]*s["cp"]*s["scalarUnitK"]*s["speedMps"]*s["diameterM"]**2
            storage=(volume@(new[:,9]-old[:,9]))*scale/dt
            with (d/"release_native"/frame.name).open() as f:
                idx=next(csv.reader(f)).index("Residual_Species_0");rhs=np.loadtxt(f,skiprows=0,delimiter=",",usecols=idx)
            reaction=-rhs[wet]*scale;close(storage,a["storageW"],1e-9)
            assert abs(a["wallW"]-a["farW"]-storage)<.005*a["wallW"]
            close(reaction.sum(),a["wallW"],1e-5*a["wallW"])
            # 完整环上三分量旋转到局部坐标，直接积分各环速度方差。
            close(trace[step],ring_rms(new[:,5:8],xyz,volume),1e-13)
            q=np.zeros(len(pts));q[wet]=-reaction/dual[wet];qobs=np.zeros(len(pts));qobs[wet]=-obs[wet,1]*scale/dual[wet]
            assert np.max(abs(q[wet]-qobs[wet])*dual[wet])/np.max(reaction)<1e-5
            temp=s["airK"]+s["scalarUnitK"]*new[:,9]
            tw=np.interp(levels,angle[line],temp[line]);nu=-np.interp(levels,angle[line],qobs[line])*s["diameterM"]/((tw-s["airK"])*(.0241*(tw/273.15)**1.5*467.15/(tw+194)))
            close(nu,a["localNu"],2e-9)
            checks.append(dict(timeDU=step*dt,storageW=float(storage),releaseCheckW=float(reaction.sum()),historyRows=rec["inner"],
                               prefixError=float(err),observer=observer,nu1198=nu[-2],nu1667=nu[-1]))
        old=new
    # 导出的每个原面和截线都在同一真实时刻；稳态标签不得出现在瞬时表中。
    assert last==proof["frame"] and proof["previous"]==rec["old"]
    val=np.zeros((len(pts),4));val[wet,0]=qobs[wet];val[:,1]=new[:,4]*s["rho"]*s["speedMps"]**2;val[:,2]=temp;val[wet,3]=temp[first]
    expected=np.array([(val[n].sum(0)+(val[quad].mean(0) if quad else 0))/3 for n,quad in face_nodes])[order]
    fields={f["key"]:f for f in page["fields"]}
    for j,key in enumerate(("heatFlux","pressureGauge","wallInput","nearWallCellTemperature")):close(fields[key]["values"],expected[:,j]/(1000 if j==0 else 1),1e-8)
    rows=table(ui/"wall.csv");assert len(rows)==4464
    for i,r in enumerate(rows):
        assert int(r["face_id"])==i+1 and r["transient_frame_sha256"]==last and "steady_frame_sha256" not in r
        assert r["joint_verified"] in ("0","false") and abs(float(r["physical_time_s"])-end*dt_s)<1e-13
        assert abs(float(r["heat_flux"])*1000-expected[i,0])<1e-7 and abs(float(r["pressure"])*1000-s["airPa"]-expected[i,1])<1e-8
    station=table(ui/"stations.csv");close([float(r["computedNu"]) for r in station],nu,2e-9)
    assert all(r["frameSHA256"]==last and abs(float(r["physical_time_s"])-end*dt_s)<1e-13 for r in station)
    cut=table(ui/"section.csv");ids=np.array([int(r["face_id"])-1 for r in cut]);x=.00137
    np.testing.assert_array_equal(np.sort(ids),np.flatnonzero((geometry[:,:,0].min(axis=1)<x)&(geometry[:,:,0].max(axis=1)>x)))
    close([float(r["value"]) for r in cut],expected[ids,0]/1000,1e-9)
    return dict(step=end,seconds=end*dt_s,frame=last,previous=rec["old"],faces=len(rows),sectionRows=len(cut),peakRSS=peak,
                nativeAdvanceSeconds=total,worstAcceptedWindow=max_window,wholeOutW=float(-area@expected[:,0]),readouts=checks,
                steadyWakeRmsOverU=ring_rms(steady[:,5:8],xyz,volume),initialWakeRmsOverU=trace[0],finalWakeRmsOverU=trace[-1],
                finalMinusSteadyNonAxisymmetricRmsOverU=ring_rms(new[:,5:8]-steady[:,5:8],xyz,volume),csvSHA256=digest(ui/"wall.csv"))


def compare_exports(ui,result):
    pair=[table(ui/m/"wall.csv") for m in ("wakeCoarse","wakeFine")]
    area=np.array([float(r["area_m2"]) for r in pair[0]])
    assert all(a["face_id"]==b["face_id"] and a["area_m2"]==b["area_m2"] for a,b in zip(*pair))
    def difference(key,scale):
        d=scale*np.array([float(b[key])-float(a[key]) for a,b in zip(*pair)])
        return float(np.sqrt(np.average(d*d,weights=area))),float(max(abs(d)))
    rms_q,max_q=difference("heat_flux",1000);rms_p,max_p=difference("pressure",1000);rms_t,max_t=difference("first_fluid_K",1)
    actual=dict(heatRmsWm2=rms_q,heatMaxWm2=max_q,pressureRmsPa=rms_p,pressureMaxPa=max_p,firstRmsK=rms_t,firstMaxK=max_t)
    target=read(result/"comparison.json")["comparison"][-1]
    for key,value in actual.items():close(value,target[key],2e-9)
    return actual


if __name__=="__main__":
    work,result,ui,dest=map(Path,sys.argv[1:5]);proof=read(ui/"ui.json");rec=read(ui/"debugger_process.json")
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values()) and not(rec["stopReason"] or rec["survivingProcessIDs"] or rec["supervisorError"])
    out={}
    for mode in proof:
        b="coarse" if mode=="wakeCoarse" else "fine"
        out[mode]=verify(work,result,ui/mode,b,mode)
    if "wakeFine" in proof and "wakeCoarse" in proof:out["exportedPair"]=compare_exports(ui,result)
    with dest.open("x",encoding="utf8") as f:json.dump(out,f,indent=2);f.write("\n")
