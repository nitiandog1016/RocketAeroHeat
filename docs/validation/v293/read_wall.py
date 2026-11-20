"""复用热敏珠的固定端点装配，对加密网逐节点核对壁反力、储热率和流动状态。"""
import argparse
import gzip
import json
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
from run_step import HERE, ROOT, RUN, ENGINE_SHA, OFFICIAL, BIN, LIMIT, configuration, inputs, observe, read, sha, bound_source
from observe_heat import BUILD, invoke, replace
from thermistorStartupSurface import checked, columns, map_surface, step_energy, write
from beadSurfaceResult import projection


def build(out,grid):
    info=inputs(grid)
    prior=ROOT/"outputs/v292_validation/endpoint_build02"
    original=read(HERE.parent/"v292/evidence/endpoint_build.json")
    file=prior/"CIntegration.cpp"
    assert sha(file)==original["sourceSHA256"] and info["nodes"]==172482
    text=file.read_text("utf8")
    assert text.count("geometry->GetnPointDomain()==148274")==2
    text=text.replace("geometry->GetnPointDomain()==148274","geometry->GetnPointDomain()==172482")
    old="config->GetTimeIter()==1 && config->GetInnerIter()==260 &&\n      config->GetnInner_Iter()==261"
    new="(config->GetTimeIter()==1 || config->GetTimeIter()==2) &&\n      config->GetInnerIter()+1==config->GetnInner_Iter()"
    (out/file.name).write_text(replace(text,old,new),encoding="utf8")
    cached=json.loads(gzip.decompress((HERE.parent/"v272/evidence/build.json.gz").read_bytes()))["cachedBuildInputs"]
    assert all(sha(Path(p))==h for p,h in cached.items())
    recipes=json.loads(gzip.decompress((HERE.parent/"v292/evidence/processes.json.gz").read_bytes()))
    for stage,limit in (("compile",120),("link",60)):
        recipe=recipes[f"endpoint_build02/{stage}_process.json"]
        assert Path(recipe["cwd"])==BUILD and recipe["exitCode"]==0
        args=[str(out/Path(s).name) if s.startswith(str(prior)+os.sep) else s for s in recipe["args"]]
        assert sha(Path(args[0]))==original["compilerSHA256"]
        invoke(stage,args,BUILD,out,limit,3)
    assert all(sha(Path(p))==h for p,h in cached.items())
    assert sum(p.stat().st_size for p in out.iterdir())<300*1024**2
    write(out/"build.json",dict(observerSHA256=sha(out/"observe_bead.exe"),sourceSHA256=sha(out/file.name),
          parentSourceSHA256=sha(file),compilerSHA256=original["compilerSHA256"],nodes=info["nodes"],
          cachedInputsVerified=len(cached)))


def diag_config(index,n,stage,settings=configuration,inner_cap=LIMIT):
    cfg=settings(index)
    count=n+1 if stage=="observe" else 1
    cfg=replace(cfg,f"INNER_ITER= {inner_cap}\n",f"INNER_ITER= {count}\n")
    if stage=="release":return replace(cfg,"(wall, VALUE, 1.0)","(wall, FLUX, 0.0)")
    if n<inner_cap:cfg=replace(cfg,"CONV_STARTITER= 50\n",f"CONV_STARTITER= {n+2}\n")
    return cfg


def sample(index,stage,area,out,grid,src,proof,settings=configuration,inner_cap=LIMIT,built=None,replay_seconds=900,extra_inputs=()):
    info=inputs(grid);s=read(proof);checked(src)
    assert s["step"]==index and s["passed"] and s["gridSHA256"]==info["gridSHA256"]
    name=f"seed_{index-1:05}.csv";end=src/f"restart_{index:05}.csv"
    assert sha(end)==s["frameSHA256"] and sha(src/name)==s["oldFrameSHA256"]
    locks=[HERE/"PLAN.md",HERE/"run_step.py",Path(__file__).resolve(),grid/"grid.json",proof,
           src/"case.cfg",src/"debugger_process.json",ROOT/"scripts/runSU2NativeDebugger.py",*extra_inputs]
    if stage=="observe":
        if built is None:built=area/"build01"
        exe=built/"observe_bead.exe"
        assert sha(exe)==read(built/"build.json")["observerSHA256"]
        locks.extend([built/"build.json",built/"CIntegration.cpp"]);seed=src/name;limit,mb=replay_seconds,400
    else:
        exe=OFFICIAL;assert sha(exe)==ENGINE_SHA
        seed=end;limit,mb=100,200
    os.link(seed,out/name);os.link(grid/"mesh.su2",out/"mesh.su2")
    cfg=out/"case.cfg";cfg.write_text(diag_config(index,s["innerRows"],stage,settings,inner_cap),encoding="ascii")
    locks.extend([exe,cfg,out/name,out/"mesh.su2"])
    os.environ["PATH"]=str(BIN)+os.pathsep+os.environ["PATH"]
    observe([str(exe),"case.cfg","--threads","1"],out,{str(p):sha(p) for p in locks},limit,mb*1024**2)
    checked(out)


def energy(index,grid,src,obs,ctrl,proof,out,settings=configuration,inner_cap=LIMIT,built=None):
    info=inputs(grid);s=read(proof)
    assert s["passed"] and s["step"]==index
    official,witness,release=[checked(f) for f in (src,obs,ctrl)]
    assert official["command"]==release["command"]==[str(OFFICIAL),"case.cfg","--threads","1"] and sha(OFFICIAL)==ENGINE_SHA
    if built is None:built=src.parent/"build01"
    assert sha(Path(witness["command"][0]))==read(built/"build.json")["observerSHA256"]
    assert all(sha(f/"mesh.su2")==info["gridSHA256"] for f in (src,obs,ctrl))
    name=f"seed_{index-1:05}.csv";end=f"restart_{index:05}.csv";n=s["innerRows"]
    assert sha(src/end)==sha(ctrl/name)==s["frameSHA256"] and sha(src/name)==sha(obs/name)==s["oldFrameSHA256"]
    assert (src/"case.cfg").read_text("ascii")==settings(index)
    assert (obs/"case.cfg").read_text("ascii")==diag_config(index,n,"observe",settings,inner_cap)
    assert (ctrl/"case.cfg").read_text("ascii")==diag_config(index,n,"release",settings,inner_cap)
    skipped=(obs/"readout_updates.csv").read_text("ascii").splitlines()
    assert skipped==[f"{index},{n},flow",f"{index},{n},scalar"]
    with np.load(grid/"nodes.npz") as z:
        ids,xyz,dual=z["ids"],z["xyz"],z["dualVolume"]
    raw=columns(obs/"before_wall.csv",["id","p","u","v","w","theta","volume"])
    raw=raw[np.argsort(raw[:,0])]
    scratch=columns(obs/end,["PointID","Pressure","Velocity_x","Velocity_y","Velocity_z","Species_0"])
    assert np.array_equal(raw[:,0],scratch[:,0]) and np.array_equal(raw[:,0],np.arange(len(xyz)))
    scratch_error=np.max(abs(raw[:,1:6]-scratch[:,1:]),axis=0)
    dual_error=float(np.max(abs(raw[:,6]/dual-1)))
    assert scratch_error.max()<2e-12 and dual_error<2e-10
    geom=read(bound_source("geometry.json"));geom["dtSeconds"]=s["dtSeconds"]
    state,reaction,balance=step_energy(src,obs,ctrl,geom,ids,xyz,witness,index,n,.025)
    result=dict(accepted=True,step=index,gridSHA256=info["gridSHA256"],acceptedFrame=str(src/end),energy=balance,
                elapsedPhysicalTimeS=s["elapsedPhysicalTimeS"],nativeDualMaxRelativeDifference=dual_error,
                scratchFrameMaxDifference=scratch_error.tolist(),skippedAppendedUpdates=skipped)
    if "verifiedPrefixRows" in s:result["verifiedPrefixRows"]=s["verifiedPrefixRows"]
    if index==2:
        # 第一列温度属于加密网首层；另取原第一层的同一 XYZ，再作原面内积分。
        os.link(bound_source("body.json"),out/"body.json");os.link(grid/"mesh.su2",out/"mesh.su2")
        result["surface"]=map_surface(out,geom,ids,xyz,state,reaction,balance,out/"wall.csv")
        body=read(out/"body.json");w,wall=projection(np.array(body["V"]),np.array(body["F"])-1,ids[:,:65],np.linspace(0,2*np.pi,49))
        with np.load(bound_source("nodes.npz")) as base:
            assert np.array_equal(ids[:,130:195],base["ids"][:,65:130])
            original=base["xyz"][base["ids"][:,65:130]]
            assert np.array_equal(xyz[ids[:,130:195]],original)
        near=np.zeros(len(xyz),int);near[ids[:,:65]]=ids[:,130:195]
        matched=geom["airK"]+(geom["nominalWallK"]-geom["airK"])*(w@state[near[wall],9])
        assert np.all((matched>=geom["airK"])&(matched<=geom["nominalWallK"]))
        np.savetxt(out/"matched_temperature.csv",np.column_stack((np.arange(len(matched)),matched)),fmt=["%d","%.17g"],
                   delimiter=",",header="faceID,temperatureAtOriginalFirstLayerK",comments="")
        result["matchedTableSHA256"]=sha(out/"matched_temperature.csv")
    write(out/"accepted.json",result)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("build","observe","release","energy"));cli.add_argument("step",type=int,choices=[1,2])
    cli.add_argument("--root",type=Path,default=RUN);cli.add_argument("--worker",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN)
    label="first" if a.step==1 else "second";grid=area/"grid01";src=area/(label+"01")
    proof=area/(label+"_state")/"state.json";out=area/("build01" if a.stage=="build" else label+"_"+a.stage)
    obs=area/(label+"_observe");ctrl=area/(label+"_release")
    if a.worker:
        assert a.stage=="energy"
        energy(a.step,grid,src,obs,ctrl,proof,out)
    else:
        out.mkdir(exist_ok=False)
        if a.stage=="build":build(out,grid)
        elif a.stage in ("observe","release"):sample(a.step,a.stage,area,out,grid,src,proof)
        else:
            files=[HERE/"PLAN.md",Path(__file__).resolve(),HERE/"run_step.py",proof,grid/"nodes.npz",grid/"grid.json",
                   ROOT/"scripts/thermistorStartupSurface.py",ROOT/"scripts/beadSurfaceResult.py",src/f"seed_{a.step-1:05}.csv",
                   src/f"restart_{a.step:05}.csv",ctrl/f"restart_{a.step:05}.csv",obs/f"restart_{a.step:05}.csv"]
            files.extend(obs/n for n in ("stages.csv","before_wall.csv","after_wall.csv","readout_updates.csv","debugger_process.json"))
            rec=observe([sys.executable,str(Path(__file__).resolve()),"energy",str(a.step),"--root",str(area),"--worker"],out,
                        {str(p):sha(p) for p in files},180,200*1024**2)
            assert not rec["stopReason"] and rec["debuggerExitCode"]==0,out/"debugger.log"
