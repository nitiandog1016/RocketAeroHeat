"""从名义双支管银球的全周流体网格启动 Test 65，不接受外部场代替初态。"""
import argparse
import os
from pathlib import Path

import numpy as np

from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import exclusive
from shortSphere65 import CASE, DRIVER, ROOT, advance, configuration, read, sha
from shortSupportMesh import BODY_HASH, prepare, reference
from thermistorOuterDomain import load_mesh


def qualify(prep):
    """从磁盘核对原生边界、逐节点热条件及原面/求解面载荷分配后才能推进。"""
    spec=read(prep/"prepared.json");s=spec["condition"]
    data,info=reference()
    assert spec["supportedSHA256"]==BODY_HASH and spec["caseSHA256"]==sha(CASE)
    assert all(sha(prep/n)==h for n,h in spec["files"].items())
    cells,coords,marks=load_mesh(prep/"mesh.su2",("wall","far"))
    with np.load(prep/"mapping.npz") as m,np.load(prep/"transfer.npz") as z:
        xyz,wet,corners,area,group,first=(m[k] for k in ("xyz","wet","corners","nativeArea","group","first"))
        np.testing.assert_array_equal(coords,xyz)
        assert len(cells)==spec["geometry"]["cells"]
        np.testing.assert_array_equal(wet,np.unique([a for face in marks["wall"] for a in face]))
        assert not np.intersect1d(wet,first).size
        for i,(face,ordered) in enumerate(zip(marks["wall"],corners)):
            assert sorted(face)==sorted(set(ordered)) and len(set(face))==len(face)
        assert len(marks["wall"])==len(corners)
        distance=np.linalg.norm(xyz[first]-xyz[wet],axis=1)
        assert distance.min()>=.00195 and distance.max()<.004
        t,j,wt,ws,old,ids=(z[k] for k in ("target","source","targetFraction","sourceFraction","targetArea","region"))
        np.testing.assert_array_equal(old,data["areaM2"]);np.testing.assert_array_equal(ids,data["regionID"])
        np.testing.assert_array_equal(ids[t],group[j]);assert min(wt)>0 and min(ws)>0
        np.testing.assert_allclose(np.bincount(t,wt),1,atol=1e-7,rtol=0)
        np.testing.assert_allclose(np.bincount(j,ws),1,atol=1e-7,rtol=0)
        # 带正负的制造通量逐区守恒；端盖还用解析正多边形面积检查。
        manufactured=(xyz[corners,0].mean(1)+xyz[corners,1].mean(1)**2-.25)*area
        mapped=np.bincount(t,ws*manufactured[j])
        err=np.bincount(group,manufactured)[1:]-np.bincount(ids,mapped)[1:]
        assert max(abs(err))<1e-15
        cap_area=32*(info["tubeDiameterM"]/2)**2*np.sin(2*np.pi/64)
        np.testing.assert_allclose([area[group==i].sum() for i in (4,5)],cap_area,rtol=1e-12,atol=0)
        profile=np.loadtxt(prep/"profile.dat");released=np.loadtxt(prep/"release.dat")
        np.testing.assert_array_equal(profile[:,0],wet);np.testing.assert_array_equal(profile[:,:2],released[:,:2])
        assert np.all(released[:,2]==0)
        source=read(CASE);temps=(np.array(source["surfaceF"])-32)*5/9+273.15
        target=np.full(len(wet),np.interp(90,source["angleDeg"],temps))
        sphere=np.unique(corners[group==1])
        angles=np.degrees(np.arccos(np.clip(xyz[sphere,0]/np.linalg.norm(xyz[sphere],axis=1),-1,1)))
        target[np.searchsorted(wet,sphere)]=np.interp(angles,source["angleDeg"],temps)
        np.testing.assert_allclose(s["airK"]+s["scalarUnitK"]*profile[:,2],target,rtol=0,atol=1e-12)
        np.testing.assert_allclose(profile[:,1]-1,profile[:,2],rtol=0,atol=1e-15)
        state=np.loadtxt(prep/"initial.csv",delimiter=",",skiprows=1)
        np.testing.assert_array_equal(state[:,0],np.arange(len(xyz)));np.testing.assert_array_equal(state[:,1:4],xyz)
        assert np.isfinite(state).all() and np.all(state[:,[4,6,7]]==0)
        fluid=np.ones(len(xyz),bool);fluid[wet]=False
        assert np.all(state[fluid,5]==-1) and np.all(state[wet,5]==0) and np.all(state[fluid,9]==0)
        np.testing.assert_array_equal(state[wet,9],profile[:,2]);np.testing.assert_array_equal(state[:,8],(1-298.15)+state[:,9])
    return dict(prepared=sha(prep/"prepared.json"),manufacturedPowerErrorsW=err.tolist(),faces=len(old),
                physicalAreaM2=old.sum(),firstNodeDistanceD=[distance.min(),distance.max()],caseSHA256=sha(CASE))


def checked_stage(prep,dest,start,build,steps,mode,species_convection="SCALAR_UPWIND"):
    if not dest.exists():
        try:
            return advance(prep,dest,build,start,steps,mode,species_convection)
        except AssertionError as error:
            if mode=="solve" and (dest/"debugger_process.json").is_file():
                with np.load(prep/"mapping.npz") as grid:
                    upper=float(grid["theta"].max())
                for path in sorted(dest.glob("restart_*.csv")):
                    scalar=np.loadtxt(path,delimiter=",",skiprows=1,usecols=9)
                    bad=(scalar < -1e-12)|(scalar > upper+1e-12)
                    if bad.any():
                        raise ValueError(f"{dest.name}/{path.name} 有 {bad.sum()} 个主温差标量超界，"
                                         f"范围 [{scalar.min():.9g}, {scalar.max():.9g}]。") from error
            raise
    if not (dest/"native.json").is_file():
        raise RuntimeError(f"{dest.name} 没有接受回执，不能从该目录接续。")
    record=read(dest/"native.json");proc=read(dest/"debugger_process.json");spec=read(prep/"prepared.json")
    exe=build/("wire_observer.exe" if mode=="observer" else "wire.exe")
    assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
    assert not(proc["stopReason"] or proc["survivingProcessIDs"] or proc["supervisorError"])
    assert proc["command"]==[str(exe),"case.cfg","--threads",str(4 if steps>20 else 1)]
    assert sha(exe)==DRIVER[exe.stem] and sha(dest/"mesh.su2")==spec["files"]["mesh.su2"]
    assert sha(dest/"wire_boundary.dat")==spec["files"]["release.dat" if mode=="release" else "profile.dat"]
    assert (dest/"case.cfg").read_text("ascii")==configuration(spec["condition"],steps,mode=="release",species_convection)
    assert record["input"]==sha(start)==sha(dest/"initial_flow.csv") and record["frame"]==sha(dest/"restart.csv")
    assert record["history"]==sha(dest/"history.csv") and all(sha(dest/n)==h for n,h in record["frames"].items())
    assert record["steps"]==steps or record["passed"]
    return record


def run(work,prep,build,resume=False,species_convection=None):
    if not resume:work.mkdir(parents=True,exist_ok=False)
    assert work.is_dir()
    with exclusive(work):
        # 接续沿用原任务登记的离散方式、体网和驱动。旧回执使用原迎风格式。
        receipt=work/"qualified.json"
        if resume and receipt.is_file():
            previous=read(receipt)
            recorded=previous.get("speciesConvection","SCALAR_UPWIND")
            if species_convection is not None and species_convection!=recorded:
                raise ValueError(f"该任务登记为 {recorded}，更换标量对流方式需新建任务。")
            species_convection=recorded
            if prep is None:prep=Path(previous["preparedPath"])
            if build is None:build=Path(previous["driverPath"])
        if species_convection is None:species_convection="SCALAR_UPWIND"
        if prep is None:
            prep=work/"prepared"
            if not prep.exists():prepare(prep)
        audit=qualify(prep)
        if build is None:
            build=ROOT/"outputs/v317_validation/build01"
            if not all((build/(k+".exe")).is_file() and sha(build/(k+".exe"))==h for k,h in DRIVER.items()):
                from buildThermistorWireDriver import build as compile_driver
                build=work/"driver"
                if not build.exists():compile_driver(build)
        assert all(sha(build/(k+".exe"))==h for k,h in DRIVER.items())
        audit.update(preparedPath=str(prep),driverPath=str(build))
        if species_convection!="SCALAR_UPWIND":audit["speciesConvection"]=species_convection
        if receipt.exists():assert read(receipt)==audit
        else:atomic_write_json(receipt,audit)
        kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build"))
        os.environ["PATH"]=str(kit/"toolchain/w64devkit/bin")+os.pathsep+os.environ["PATH"]
        start=prep/"initial.csv";trend=[]
        try:
            for i in range(1,13):
                dest=work/f"steady{i:02d}"
                rec=checked_stage(prep,dest,start,build,200,"solve",species_convection)
                worst=max(rec["window"].values());trend.append(worst)
                print(dest.name,rec["steps"],worst,flush=True)
                start=dest/"restart.csv"
                if rec["passed"]:break
                if len(trend)>=3 and trend[-1]>=trend[-2]>=trend[-3]:
                    raise RuntimeError("连续两段末窗不改善")
            else:raise RuntimeError("已用完登记的稳态步数")
            for name,seed,n,mode in [("reaction",start,1,"release"),("plus20",start,20,"solve"),
                    ("reaction20",dest/"plus20/restart.csv",1,"release"),("control",start,1,"solve"),("observer",start,1,"observer")]:
                rec=checked_stage(prep,dest/name,seed,build,n,mode,species_convection)
                if name=="plus20":assert rec["passed"]
            atomic_write_json(work/"endpoints.json",dict(stage=dest.name,frame=sha(start),case=sha(CASE),residuals=trend))
        except Exception as error:
            failure=work/"failure.json"
            if not failure.exists():atomic_write_json(failure,dict(stage=dest.name,acceptedWindows=trend,error=repr(error)))
            raise


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work",type=Path);parser.add_argument("--prepared",type=Path);parser.add_argument("--driver",type=Path)
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--species-convection",choices=("SCALAR_UPWIND","BOUNDED_SCALAR"),
                        help="新算默认 SCALAR_UPWIND；接续沿用原任务的选择。")
    args=parser.parse_args()
    work=args.work.resolve();work.relative_to(ROOT/"outputs")
    run(work,args.prepared.resolve() if args.prepared else None,args.driver.resolve() if args.driver else None,args.resume,args.species_convection)
