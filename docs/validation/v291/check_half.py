"""核对热敏珠半步的残差、完整时间层与热量账；只将合格端点放行给下一步。"""
import argparse
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from run_half import ENGINE_SHA, HERE, MAX_INNER, OFFICIAL, ROOT, RUN, SOURCE, config, input_file, observe, read, sha
from observe_half import observer_config, release_config
from thermistorStartupSurface import checked, columns, map_surface, step_energy, write

sys.path[:0]=[str(HERE.parent/version) for version in ("v287","v288","v289")]
from audit_step import RMS, frame
from check_fields import rhs_at_endpoint


def state(index,source,out,settings=config,inner_limit=MAX_INNER):
    rec=checked(source)
    assert rec["command"]==[str(OFFICIAL),"case.cfg","--threads","1"] and sha(OFFICIAL)==ENGINE_SHA
    assert (source/"case.cfg").read_text("ascii")==settings(index)
    hist=source/f"history_{index:05}.csv"
    end=source/f"restart_{index:05}.csv"
    old=source/f"seed_{index-1:05}.csv"
    assert end.name in (source/"debugger.log").read_text(errors="replace")
    h=columns(hist,["Time_Iter","Inner_Iter"]+RMS)
    assert 20<=len(h)<=inner_limit and np.all(h[:,0]==index) and np.isfinite(h).all()
    np.testing.assert_array_equal(h[:,1],np.arange(len(h)))
    with np.load(input_file("nodes.npz")) as grid:
        ids,xyz=grid["ids"],grid["xyz"]
    before=frame(old,xyz)
    after=frame(end,xyz)
    assert not np.array_equal(before[:,4:],after[:,4:]) and np.array_equal(ids[0],ids[-1])
    geom=read(input_file("geometry.json"))
    pressure=geom["airPa"]+geom["rho"]*geom["speedMps"]**2*after[:,4]
    theta=after[:,9]
    wall=ids[:-1,:65]
    slip=float(np.max(abs(after[wall,5:8])))
    wall_error=float(np.max(abs(theta[wall]-1)))
    rhs=rhs_at_endpoint(end,xyz,h[-1,2:])
    window=h[-20:,2:].max(axis=0)
    good=np.all(h[:,2:]<-8,axis=1)
    passed=(window.max()<-8 and pressure.min()>0 and theta.min()>=-1e-12 and theta.max()<=1+1e-12
            and slip<=1e-12 and wall_error<=1e-12 and rhs["maxLogDifferenceFromHistory"]<1e-7)
    write(out/"state.json",dict(step=index,passed=bool(passed),innerRows=len(h),
          last20MaxLogRMS=dict(zip(RMS,window)),finalLogRMS=dict(zip(RMS,h[-1,2:])),
          trailingPassingRows=int(len(h)-1-np.flatnonzero(~good)[-1]) if not good.all() else len(h),
          maxWallVelocityOverU=slip,wallScalarError=wall_error,scalarRange=[float(theta.min()),float(theta.max())],
          pressureRangePa=[float(pressure.min()),float(pressure.max())],nodeRHS=rhs,
          dtSeconds=geom["dtSeconds"]/2,elapsedPhysicalTimeS=index*geom["dtSeconds"]/2,
          oldFrameSHA256=sha(old),frameSHA256=sha(end),historySHA256=sha(hist),
          configSHA256=sha(source/"case.cfg"),meshSHA256=sha(source/"mesh.su2"),seconds=rec["elapsedSeconds"]))


def energy(index,source,out,proof_path,observer,release,settings=config,inner_limit=MAX_INNER):
    proof=read(proof_path)
    assert proof["passed"] and proof["step"]==index
    official,obs,ctrl=[checked(folder) for folder in (source,observer,release)]
    old=f"seed_{index-1:05}.csv"
    end=f"restart_{index:05}.csv"
    assert official["command"]==ctrl["command"]==[str(OFFICIAL),"case.cfg","--threads","1"]
    assert sha(OFFICIAL)==ENGINE_SHA and sha(source/end)==proof["frameSHA256"]==sha(release/old)
    assert sha(source/old)==proof["oldFrameSHA256"]==sha(observer/old)
    assert all(sha(f/"mesh.su2")==sha(input_file("mesh.su2")) for f in (source,observer,release))
    assert (source/"case.cfg").read_text("ascii")==settings(index)
    assert (observer/"case.cfg").read_text("ascii")==observer_config(index,proof["innerRows"],settings,inner_limit)
    assert (release/"case.cfg").read_text("ascii")==release_config(index,settings,inner_limit)
    with np.load(input_file("nodes.npz")) as grid:
        ids,xyz=grid["ids"],grid["xyz"]
    geom=read(input_file("geometry.json"))
    geom["dtSeconds"]=proof["dtSeconds"]
    assert geom["dtSeconds"]==.025*geom["transverseDiameterM"]/geom["speedMps"]
    values,reaction,balance=step_energy(source,observer,release,geom,ids,xyz,obs,index,proof["innerRows"],.025)
    report=dict(accepted=True,step=index,acceptedFrame=str(source/end),elapsedPhysicalTimeS=proof["elapsedPhysicalTimeS"],energy=balance)
    if index==2:
        report["surface"]=map_surface(SOURCE,geom,ids,xyz,values,reaction,balance,out/"wall.csv")
    write(out/"energy.json",report)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=["state","energy"])
    cli.add_argument("step",type=int,choices=[1,2])
    cli.add_argument("source",type=Path)
    cli.add_argument("output",type=Path)
    cli.add_argument("--proof",type=Path)
    cli.add_argument("--observer",type=Path)
    cli.add_argument("--release",type=Path)
    cli.add_argument("--worker",action="store_true")
    a=cli.parse_args()
    folder=a.source.resolve();out=a.output.resolve()
    folder.relative_to(RUN);out.relative_to(RUN)
    if a.worker:
        if a.stage=="state":
            state(a.step,folder,out)
        else:
            energy(a.step,folder,out,a.proof.resolve(),a.observer.resolve(),a.release.resolve())
    else:
        out.mkdir(parents=True,exist_ok=False)
        files=[HERE/"PLAN.md",Path(__file__).resolve(),HERE/"observe_half.py",HERE/"run_half.py",
               ROOT/"scripts/thermistorStartupSurface.py",ROOT/"scripts/beadSurfaceResult.py",
               HERE.parent/"v289/check_fields.py",HERE.parent/"v287/audit_step.py",
               input_file("nodes.npz"),input_file("body.json"),input_file("geometry.json"),
               folder/"case.cfg",folder/"debugger_process.json",folder/f"seed_{a.step-1:05}.csv",
               folder/f"restart_{a.step:05}.csv",folder/f"history_{a.step:05}.csv"]
        if a.stage=="energy":
            files.append(a.proof.resolve())
            files.extend(a.observer.resolve()/n for n in ("debugger_process.json","before_wall.csv","after_wall.csv","stages.csv"))
            files.extend(a.release.resolve()/n for n in ("debugger_process.json",f"restart_{a.step:05}.csv"))
        command=[sys.executable,str(Path(__file__).resolve()),a.stage,str(a.step),str(folder),str(out),"--worker"]
        for key in ("proof","observer","release"):
            path=getattr(a,key)
            if path:command.extend(["--"+key,str(path.resolve())])
        rec=observe(command,out,{str(p):sha(p) for p in files},120,100*1024**2)
        assert not rec["stopReason"] and rec["debuggerExitCode"]==0,out/"debugger.log"
