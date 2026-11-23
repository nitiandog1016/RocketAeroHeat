"""从通过残差的真实端点执行定温延续、只读通量观察与独立释放壁行。"""
import argparse
import importlib.util
import os
from pathlib import Path

import numpy as np

from run import BASE, CASE, HERE, PREP, ROOT, configuration, observe, read, sha, write_json
from build import SCALAR, SPECIES, arguments, invoke, once
from shortSphere65 import DRIVER

DRIVER_DIR=ROOT/"outputs/v317_validation/build01"
KIT=ROOT/"tools/runtime/su2-thermal-build"
OLD=HERE.parent/"v326"


def observer_build(out,production):
    out.mkdir(parents=True,exist_ok=False)
    spec=importlib.util.spec_from_file_location("short326_probe",OLD/"build_probe.py")
    probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    pro=read(production/"build.json")
    for name in ("silver.hpp","capacity.hpp","CSpeciesSolver.cpp","CScalarSolver.inl"):
        assert sha(production/name)==pro["sources"][name]
        (out/name).write_bytes((production/name).read_bytes())
    header=out/"silver.hpp"
    header.write_text(once(header.read_text("utf8"),"  double lo=1e30,",r"""
  if(!predicted && step==0) {
    auto* f=file("nodes.bin","wb");
    for(const auto& n:nodesData) {
      const double a[]={double(n.global),n.old,n.low,n.final,double(n.wall),n.volume,n.dt,n.native,n.raw,
          n.correction,n.limited,n.farMass,n.farFlux,n.divergence,n.budget.give,n.budget.take};
      if(std::fwrite(a,sizeof(double),16,f)!=16) std::abort();
    }
    std::fclose(f); f=file("edges.bin","wb");
    for(const auto& e:edges) {
      const double a[]={double(nodesData[e.i].global),double(nodesData[e.j].global),e.mass,e.high,e.correction,e.alpha};
      if(std::fwrite(a,sizeof(double),6,f)!=6) std::abort();
    }
    std::fclose(f);
  }
  double lo=1e30,"""),encoding="utf8")
    compiler=KIT/"toolchain/w64devkit/bin/g++.exe";work=KIT/"build-v117";src=KIT/"source"
    (out/"probe.hpp").write_bytes((OLD/"probe.hpp").read_bytes())
    assert sha(src/probe.FLOW)==read(ROOT/"outputs/v326_validation/build01/build.json")["sources"][probe.FLOW]
    (out/"CIncEulerSolver.cpp").write_text(probe.flow_source((src/probe.FLOW).read_text("utf8")),encoding="utf8")
    entries=read(work/"compile_commands.json");objs=[]
    for rel in (SPECIES,probe.FLOW):
        entry,=[e for e in entries if e["file"].replace("\\","/").endswith("/"+rel)]
        args=arguments(entry["command"]);args[0]=compiler;obj=out/(Path(rel).stem+".obj")
        for flag,path in (("-MQ",obj),("-MF",str(obj)+".d"),("-o",obj)):args[args.index(flag)+1]=path
        args[-1:]=["-iquote",(src/rel).parent,"-iquote",(src/SCALAR).parent,out/Path(rel).name]
        invoke(obj.stem,args,work,out,180,3,bin_directory=compiler.parent);objs.append(obj)
    link=read(DRIVER_DIR/"wire_process.json")["args"];exe=out/"readout.exe"
    invoke("readout",[compiler,"-o",exe,link[3],*objs,*link[4:]],work,out,180,3,bin_directory=compiler.parent)
    write_json(out/"build.json",dict(exe=sha(exe),base=pro["exe"],observer=sha(Path(__file__)),
               generated={p.name:sha(p) for p in out.iterdir() if p.suffix in (".cpp",".hpp",".inl",".obj")}))


def advance(folder,start,exe,steps,release=False):
    folder.mkdir()
    meta=read(PREP/"prepared.json");inputs={"mesh.su2":PREP/"mesh.su2","initial_flow.csv":start,
            "wire_boundary.dat":PREP/("release.dat" if release else "profile.dat")}
    for name,source in inputs.items():os.link(source,folder/name)
    assert sha(folder/"mesh.su2")==meta["files"]["mesh.su2"] and sha(folder/"wire_boundary.dat")==meta["files"][inputs["wire_boundary.dat"].name]
    cfg=folder/"case.cfg";cfg.write_text(configuration(meta["condition"],steps,release,"BOUNDED_SCALAR"),encoding="ascii")
    paths=[exe,cfg,*[folder/n for n in inputs],HERE/"ENDPOINT_PLAN.md",Path(__file__).resolve(),PREP/"prepared.json",CASE]
    if exe.name=="readout.exe":
        point=folder/"probe.in";point.write_text("0 -1 0.46839506172839518\n",encoding="ascii");paths.append(point)
    os.environ["PATH"]=str(KIT/"toolchain/w64devkit/bin")+os.pathsep+os.environ["PATH"]
    r=observe([str(exe),"case.cfg","--threads","4"],folder,{str(p):sha(p) for p in paths},150,450*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"])
    assert "Exit Success (SU2_CFD)" in (folder/"debugger.log").read_text(errors="replace")
    history=np.loadtxt(folder/"history.csv",delimiter=",",skiprows=1,ndmin=2)
    assert len(history)==steps and np.isfinite(history).all()
    assert sha(folder/"restart.csv")==sha(folder/f"restart_{steps-1:06d}.csv")
    return dict(start=sha(start),end=sha(folder/"restart.csv"),exe=sha(exe),history=sha(folder/"history.csv"),
                seconds=r["elapsedSeconds"],peakBytes=r["peakFamilyRSSBytes"],window=history[-20:,3:9].max(axis=0).tolist())


def sample(out,stage,production,observer):
    end=stage/"restart.csv";rec=read(stage/"stage.json");assert rec["passed"] and sha(end)==rec["frames"]["restart_000199.csv"]
    out.mkdir(exist_ok=False)
    src=production/"silver.exe";probe=observer/"readout.exe";release=DRIVER_DIR/"wire.exe"
    assert sha(src)==read(production/"build.json")["exe"]==read(observer/"build.json")["base"]
    assert sha(probe)==read(observer/"build.json")["exe"] and sha(release)==DRIVER["wire"]
    records={}
    for name,seed,driver,steps,let in (("plus20",end,src,20,False),("control",end,src,1,False),
            ("early",end,probe,1,False),("releaseEarly",end,release,1,True),
            ("late",out/"plus20/restart.csv",probe,1,False),("releaseLate",out/"plus20/restart.csv",release,1,True)):
        records[name]=advance(out/name,seed,driver,steps,let)
        write_json(out/name/"sample.json",records[name])
    write_json(out/"samples.json",records)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("action",choices=("build","sample"))
    p.add_argument("name");p.add_argument("--observer",default="readBuild01");p.add_argument("--stage",default="repeat03")
    a=p.parse_args();assert all(x.isidentifier() for x in (a.name,a.observer,a.stage))
    if a.action=="build":observer_build(BASE/a.name,BASE/"build01")
    else:sample(BASE/a.name,BASE/a.stage,BASE/"build01",BASE/a.observer)
