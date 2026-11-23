"""在同一 BDF1 端点重算储热、独立壁反力、逐边限幅及真实远场能量。"""
from pathlib import Path
import numpy as np
from build_time import PREP, read, sha, write_json
from thermistorWireSurface import residual


def endpoint(work,branch,step):
    folder=work/branch/f"step{step}";release=folder/"release";job=read(work/"job.json")
    cfg=job["condition"];meta=read(folder/"native.json");dt=meta["dt"]
    frame=folder/f"restart_{step:05d}.csv";seed=folder/f"seed_{step-1:05d}.csv"
    for place in (folder,release):
        rec=read(place/"debugger_process.json")
        assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
        assert not(rec["supervisorError"] or rec["survivingProcessIDs"] or rec["stopReason"])
    assert meta["old"]==sha(seed) and meta["frame"]==sha(frame)==read(release/"native.json")["old"]
    assert max(meta["window"]) < -8 and meta["rawMax"]<1e-8
    state=np.loadtxt(frame,delimiter=",",skiprows=1,usecols=range(10))
    previous=np.loadtxt(seed,delimiter=",",skiprows=1,usecols=9)
    raw_nodes=np.fromfile(folder/"nodes.bin","<f8").reshape(-1,18);order=np.argsort(raw_nodes[:,0])
    n=raw_nodes[order];N=len(n);np.testing.assert_array_equal(n[:,0],state[:,0])
    assert N==job["nodes"] and np.isfinite(n).all() and (n[:,5:7]>0).all()
    with np.load(PREP/"mapping.npz") as m:wet=m["wet"];theta=m["theta"]
    wall=n[:,4].astype(bool);np.testing.assert_array_equal(np.flatnonzero(wall),wet)
    np.testing.assert_allclose(state[wet,9],theta,atol=2e-14,rtol=0)
    def native_row(name):return np.fromfile(folder/(name+".bin"),"<f8")[order]
    np.testing.assert_array_equal(n[:,1],n[:,3])
    np.testing.assert_allclose(n[:,1],state[:,9],rtol=0,atol=2e-14)
    np.testing.assert_allclose(n[:,16],previous,rtol=0,atol=2e-14)
    for place in (folder,release):
        flow=np.fromfile(place/"read_flow.bin","<f8").reshape(-1,5)
        flow=flow[np.argsort(flow[:,0])];np.testing.assert_array_equal(flow[:,0],state[:,0])
        np.testing.assert_allclose(flow[:,1:],state[:,4:8],atol=1e-12,rtol=1e-13)
    scale=cfg["rho"]*cfg["cp"]*cfg["speedMps"]*cfg["diameterM"]**2*cfg["scalarUnitK"]
    storage=n[:,5]*(state[:,9]-previous)/dt
    native_storage=native_row("after_time")-native_row("before_time")
    np.testing.assert_allclose(n[:,17],storage,atol=2e-13,rtol=1e-12)
    np.testing.assert_allclose(native_storage,storage,atol=2e-13,rtol=1e-12)
    actual_storage=native_storage.sum()*scale;storage_w=storage.sum()*scale
    assert abs(actual_storage-storage_w)<1e-9
    before=native_row("before_wall");after=native_row("after_wall")
    np.testing.assert_array_equal(before[~wall],after[~wall]);assert np.all(after[wall]==0)
    pre=before[wet]-after[wet]
    independent=-residual(release/f"restart_{step:05d}.csv")[wet]
    release_error=max(abs(independent-pre))/max(abs(pre));assert release_error<1e-5,release_error
    e=np.fromfile(folder/"edges.bin","<f8").reshape(-1,6)
    assert len(e)==921290 and np.isfinite(e).all()
    i,j=e[:,:2].astype(int).T;mass,high,corr,alpha=e[:,2:].T
    def acc(a,b):return np.bincount(i,a,minlength=N)+np.bincount(j,b,minlength=N)
    def net(a):return acc(a,-a)
    div=-n[:,1]*(net(mass)+n[:,11]);spatial=net(high)+div+n[:,12]
    np.testing.assert_allclose(spatial+storage,n[:,8],atol=2e-13,rtol=1e-12)
    np.testing.assert_allclose(spatial+storage,before,atol=2e-13,rtol=1e-12)
    low=n[:,2];upper=max(theta);assert low.min()>=-1e-12 and low.max()<=upper+1e-12
    shares=[]
    for outgoing,room in ((acc(np.maximum(corr,0),np.maximum(-corr,0)),np.minimum(n[:,1],low)),
                          (acc(np.maximum(-corr,0),np.maximum(corr,0)),upper-np.maximum(n[:,1],low))):
        active=(outgoing>0)&~wall;a=np.ones(N)
        a[active]=np.minimum(1,n[active,5]/n[active,6]*np.maximum(room[active],0)/outgoing[active]);shares.append(a)
    give,take=shares;check=np.where(corr>=0,np.minimum(give[i],take[j]),np.minimum(take[i],give[j]))
    np.testing.assert_allclose(check,alpha,rtol=1e-12,atol=1e-12)
    correction=net(corr*(check-1))
    np.testing.assert_allclose(correction,n[:,10]-n[:,9],atol=2e-14,rtol=1e-12)
    node_in=-(pre+correction[wet])*scale;wall_w=-node_in.sum()
    far=np.loadtxt(folder/"flow_far.csv",delimiter=",");ids=far[:,0].astype(int)
    assert len(np.unique(ids))==len(ids)==2690 and not wall[ids].any()
    np.testing.assert_array_equal(far[:,1],n[ids,11]);far_w=(np.maximum(far[:,1],0)@state[ids,9])*scale
    assert abs(far_w-n[:,12].sum()*scale)<1e-10
    total=spatial+correction+storage;free=total[~wall].sum()*scale;source=div.sum()*scale
    closure=free+wall_w-far_w-actual_storage-source
    gap=(wall_w-far_w-storage_w)/wall_w
    assert wall_w>0 and abs(closure)<1e-9 and abs(gap)<.005,(closure,gap)
    assert abs(total[wet].sum()*scale-wall_w)<1e-10
    fields=dict(wet=wet,nodeInW=node_in,rawWall=pre,released=independent,edgeCorrection=correction[wet],
                farID=ids,farMass=far[:,1],farTheta=state[ids,9])
    budget=folder/"wall_budget.npz"
    if budget.exists():
        with np.load(budget) as archived:
            for key,value in fields.items():np.testing.assert_array_equal(archived[key],value)
    else:np.savez_compressed(budget,**fields)
    files=["nodes.bin","edges.bin","before_time.bin","after_time.bin","before_wall.bin","after_wall.bin","read_flow.bin",
           "flow_far.csv","stages.csv","far_pairs.csv","wall_budget.npz","native.json","release/native.json"]
    result=dict(**meta,wallW=wall_w,farW=far_w,storageW=storage_w,storageFromNativeW=actual_storage,
                storageDifferenceW=actual_storage-storage_w,relativeBalance=gap,freeEquationW=free,divergenceSourceW=source,
                assemblyGapW=closure,releasePeakRelativeError=release_error,scale=scale,
                edgeCount=len(e),farNodes=len(ids),files={name:sha(folder/name) for name in files})
    file=folder/"energy.json"
    if file.exists():assert read(file)==result
    else:write_json(file,result)
    return result
