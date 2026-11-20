"""将 6 m/s 热敏珠已接受的启动终场和同阶段原生壁反力保守传递到完整原面。"""
import argparse
import csv
import json
from pathlib import Path
import shutil
import zipfile

import numpy as np
from beadSurfaceResult import projection
from thermistorBead import ROOT, B, sha


def read(path):
    return json.loads(path.read_bytes())


def write(path, data):
    with path.open("x",encoding="utf8") as f:
        json.dump(data,f,ensure_ascii=False,allow_nan=False,indent=2)
        f.write("\n")


def columns(path,names):
    with path.open(encoding="ascii") as f:
        header=[v.strip() for v in next(csv.reader(f,skipinitialspace=True))]
        return np.loadtxt(f,delimiter=",",usecols=[header.index(k) for k in names],ndmin=2)


def checked(folder):
    receipt=read(folder/"debugger_process.json")
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and not receipt["stopReason"]
    assert not receipt["survivingProcessIDs"] and not receipt["supervisorError"] and all(receipt["inputsStable"].values())
    assert "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")
    for p,h in receipt["inputSHA256"].items():
        assert sha(Path(p))==h,p
    return receipt


def energy(source,observer,release,geom,ids,xyz):
    obs=checked(observer);ctrl=checked(release)
    accepted=read(ROOT/"docs/validation/v289/evidence/first_step_200.json")
    assert accepted["acceptedFirstStep"] and accepted["referencePrefixByteIdentical"]
    for name in ("seed_00000.csv","restart_00001.csv","history_00001.csv"):
        assert sha(source/name)==accepted["sha256"][name]
    assert sha(observer/"seed_00000.csv")==accepted["sha256"]["seed_00000.csv"]
    assert sha(release/"seed_00000.csv")==accepted["sha256"]["restart_00001.csv"]
    assert all(sha(d/"mesh.su2")=="673151f2681de2a150f7de689bea3096e3b6fa6d9b91ddb80dc635e563fe951c" for d in (source,observer,release))
    assert ctrl["inputSHA256"][ctrl["command"][0]]=="01f22d0c9de2e9a5feea666f67adf3957f503d2e14ca6fd998a9a28a11417bf7"
    base_cfg=(source/"case.cfg").read_text("ascii")
    assert (observer/"case.cfg").read_text("ascii")==base_cfg.replace("INNER_ITER= 200","INNER_ITER= 201")
    assert (release/"case.cfg").read_text("ascii")==base_cfg.replace("INNER_ITER= 200","INNER_ITER= 1").replace("(wall, VALUE, 1.0)","(wall, FLUX, 0.0)")
    return step_energy(source,observer,release,geom,ids,xyz,obs,step=1,inner=200,dt_nd=.05)


def step_energy(source,observer,release,geom,ids,xyz,obs,step,inner,dt_nd):
    """在已绑定的同一热敏珠时间层比对原生装配、完整旧新态与独立壁读数。"""
    end_name=f"restart_{step:05}.csv"
    old_name=f"seed_{step-1:05}.csv"
    history_name=f"history_{step:05}.csv"
    fields=["PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"]
    new=columns(source/end_name,fields)
    old=columns(source/old_name,fields)
    n=len(xyz);node_ids=np.arange(n)
    for state in (new,old):
        assert np.isfinite(state).all() and np.array_equal(state[:,0],node_ids)
        np.testing.assert_allclose(state[:,1:4],xyz,rtol=1e-14,atol=2e-14)
    wids=np.unique(ids[:,:65])
    sample=columns(observer/"before_wall.csv",["id","rhs","theta","theta_old","volume","p","u","v","w","rho","dt"])
    after=columns(observer/"after_wall.csv",["id","rhs"])
    sample=sample[np.argsort(sample[:,0])];after=after[np.argsort(after[:,0])]
    assert np.array_equal(sample[:,0],node_ids) and np.array_equal(after[:,0],node_ids)
    assert np.all(sample[:,4]>0) and np.all(sample[:,9]==1) and np.all(sample[:,10]==dt_nd)
    matches=dict(oldTheta=float(np.max(abs(sample[:,3]-old[:,9]))),newTheta=float(np.max(abs(sample[:,2]-new[:,9]))),
                 pressure=float(np.max(abs(sample[:,5]-new[:,4]))),velocity=float(np.max(abs(sample[:,6:9]-new[:,5:8]))))
    assert matches["oldTheta"]<=1e-13 and matches["newTheta"]<=2e-10 and max(matches["pressure"],matches["velocity"])<=1e-8,matches
    assert np.all(sample[wids,2:4]==1)
    assert np.all(after[wids,1]==0)
    delta=sample[:,1]-after[:,1]
    assert np.count_nonzero(delta)==np.count_nonzero(delta[wids]) and np.isfinite(delta).all()
    stages={i:{} for i in (inner-1,inner)}
    for i,key,value in csv.reader((observer/"stages.csv").read_text("ascii").splitlines()):
        assert key not in stages[int(i)]
        stages[int(i)][key]=float(value)
    assert all(list(t)==["start","beforeTime","afterTime","beforeFar","afterFar","beforeWall","afterWall"] for t in stages.values())
    speed=geom["speedMps"];diam=geom["transverseDiameterM"];dt=geom["dtSeconds"]
    delta_t=geom["nominalWallK"]-geom["airK"]
    scale=geom["rho"]*geom["cp"]*speed*diam**2*delta_t
    def wall(t):
        return (t["beforeWall"]-t["afterWall"])*scale
    previous=stages[inner-1]
    now=stages[inner];wall_w=wall(now);assert wall_w>0
    far=(now["afterFar"]-now["beforeFar"])*scale
    storage=(now["afterTime"]-now["beforeTime"])*scale
    from_frames=float(np.dot(sample[:,4],new[:,9]-old[:,9])*geom["rho"]*geom["cp"]*delta_t*diam**3/dt)
    assert abs(from_frames-storage)<1e-9 and abs(delta.sum()*scale-wall_w)<1e-12
    assert max(abs(now["beforeTime"]),abs(now["start"]),abs(now["beforeWall"]-now["afterFar"]),
               abs(now["beforeFar"]-now["afterTime"]))*scale<1e-10
    closure=(wall_w-far-storage)/wall_w
    assert abs(closure)<.005 and abs(wall(previous)/wall_w-1)<.002
    released=columns(release/end_name,["PointID","Residual_Species_0"])
    assert np.array_equal(released[:,0],node_ids)
    max_error=float(np.max(abs(released[wids,1]+delta[wids]))/np.max(abs(delta[wids])))
    assert max_error<1e-5,max_error
    keys=["Time_Iter","Inner_Iter"]+["rms["+s+"]" for s in ("P","U","V","W","h","rho*Y_0")]
    official=columns(source/history_name,keys)
    observed=columns(observer/history_name,keys)
    assert len(official)==inner and len(observed)==inner+1 and np.array_equal(official[:,:2],observed[:-1,:2])
    history_diff=float(np.max(abs(official[:,2:]-observed[:-1,2:])))
    assert history_diff<.01,history_diff
    receipt=dict(wallReleaseW=wall_w,farOutW=far,airStorageW=storage,independentStorageW=from_frames,
        imbalanceW=wall_w-far-storage,relativeImbalance=closure,previousAssemblyReleaseW=wall(previous),
        relativeWallChange=abs(wall(previous)/wall_w-1),releasePeakRelativeError=max_error,
        observedHistoryMaxLogDifference=history_diff,assemblyStateErrors=matches,stages=stages,
        observerSHA256=obs["inputSHA256"][obs["command"][0]],observedWallSHA256=sha(observer/"before_wall.csv"),
        reactionCheckFrameSHA256=sha(release/end_name),acceptedFrameSHA256=sha(source/end_name),
        oldFrameSHA256=sha(source/old_name),dtSeconds=dt)
    return new,-scale*delta,receipt


def collect(source,observer,release,out):
    out.mkdir(parents=True,exist_ok=False)
    geom=read(source/"geometry.json")
    with np.load(source/"nodes.npz") as grid:
        ids,xyz=grid["ids"],grid["xyz"]
    assert ids.shape==(49,3185) and np.array_equal(ids[0],ids[-1])
    state,reaction,balance=energy(source,observer,release,geom,ids,xyz)
    summary=map_surface(source,geom,ids,xyz,state,reaction,balance,out/"wall.csv")
    write(out/"summary.json",summary)
    write(out/"case.json",{k:geom[k] for k in ("speedMps","airK","airPa","nominalWallK","dtSeconds","transverseDiameterM")})
    shutil.copyfile(source/"body.json",out/"body.json")
    with zipfile.ZipFile(out/"startup_6ms.zip","x") as bundle:
        for key in ("body.json","case.json","summary.json","wall.csv"):
            info=zipfile.ZipInfo(key,date_time=(2026,1,1,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            bundle.writestr(info,(out/key).read_bytes())
    return summary


def map_surface(source,geom,ids,xyz,state,reaction,balance,wall):
    """把同一时刻的热敏珠壁反力和节点状态按真实方位映回原面，并核对功率。"""
    shape=read(source/"body.json");v=np.array(shape["V"]);f=np.array(shape["F"])-1;tris=v[f]
    assert sha(source/"body.json")=="5e66bb2f70f0774b8ebb87a77b76061fc84b89b9f379d54e90879e81628bf62f"
    area=np.linalg.norm(np.cross(tris[:,1]-tris[:,0],tris[:,2]-tris[:,0]),axis=1)/2
    assert len(f)==5216 and np.isclose(area.sum(),3.0898302613936234e-6,rtol=1e-12,atol=0)
    weights,wall_ids=projection(v,f,ids[:,:65],np.linspace(0,2*np.pi,49))
    # 线性 y/B、z/B 参考场可独立识别周向误折叠、符号和接缝错误。
    err=np.max(abs(weights@(xyz[wall_ids,1:3]*geom["transverseDiameterM"]/B)-tris.mean(axis=1)[:,1:3]/B),axis=0)
    assert err.max()<.005,err
    dual=weights.T@area;assert np.all(dual>0) and len(wall_ids)==48*63+2
    q=weights@(reaction[wall_ids]/dual)
    near=np.zeros(len(xyz),int);near[ids[:,:65]]=ids[:,65:130]
    dt=geom["nominalWallK"]-geom["airK"]
    first=geom["airK"]+dt*(weights@state[near[wall_ids],9])
    gauge=geom["rho"]*geom["speedMps"]**2*(weights@state[wall_ids,4])
    p=geom["airPa"]+gauge;power=float(area@q)
    assert abs(power+balance["wallReleaseW"])<1e-11*abs(power)
    assert np.isfinite(q).all() and (p>0).all() and first.min()>=geom["airK"]-1e-10 and first.max()<=geom["nominalWallK"]+1e-10
    np.savetxt(wall,np.column_stack((np.arange(len(f)),area,p,gauge,q,first)),
        fmt=["%d"]+["%.17g"]*5,delimiter=",",header="faceID,areaM2,pressurePa,pressureGaugePa,heatFluxWm2,nearWallCellTemperatureK",comments="")
    summary=dict(originalFaces=len(f),areaM2=float(area.sum()),powerIntoWallW=power,energy=balance,
        bodySHA256=sha(source/"body.json"),wallTableSHA256=sha(wall),meshSHA256=sha(source/"mesh.su2"),
        projectionMaxError=err.tolist(),wallNodes=len(wall_ids),wallPressurePa=[float(p.min()),float(p.max())],
        heatFluxWm2=[float(q.min()),float(q.max())],nearWallK=[float(first.min()),float(first.max())],
        peakOutwardFaceID=int(np.argmin(q)),peakFaceCentroidM=tris[np.argmin(q)].mean(axis=0).tolist(),
        jointExperimentVerifiedAreaM2=0,fullSurfaceComputed=True)
    return summary


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    for key in ("source","observer","release","output"):cli.add_argument(key,type=Path)
    a=cli.parse_args();out=a.output.resolve();out.relative_to(ROOT/"outputs")
    collect(a.source.resolve(),a.observer.resolve(),a.release.resolve(),out)
