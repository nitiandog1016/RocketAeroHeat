"""在相同的 50 µs 工况中开关只读观察，绑定终帧、父初场和原生运行身份。"""
import argparse
import os
from pathlib import Path
import shutil
import sys

from build_observer import BASE,BIN,DOC,ROOT,sha
sys.path.insert(0,str(DOC.parent/"v366"))
from physical_step import inspect,observe,read,save,terminal

PARENT=ROOT/"outputs/v366_validation/time03/full"
KINDS=("off","on")


def replay(folder,build,kind,threads):
    folder=folder.resolve();folder.relative_to(BASE);build=build.resolve();build.relative_to(BASE)
    prior=inspect(PARENT,"full","newton")
    assert not prior["innerPassed"] and prior["innerIterations"]==250
    receipt=read(build/"build.json")
    exe=build/"tail_observe.exe"
    assert sha(exe)==receipt["executableSHA256"]
    folder.mkdir(parents=True,exist_ok=False)
    locks={}
    for name in ("mesh.su2","case.cfg","seed_00000.csv"):
        src=PARENT/name
        dst=folder/name
        shutil.copyfile(src,dst)
        assert sha(dst)==sha(src)
        locks[str(dst)]=sha(src)
    for p in (exe,build/"build.json",DOC/"PLAN.md",DOC/"assembly.hpp",Path(__file__).resolve(),DOC/"build_observer.py",
              DOC.parent/"v366/physical_step.py",ROOT/"scripts/runSU2NativeDebugger.py"):
        locks[str(p)]=sha(p)
    # 本进程只为其子进程装载本地编译运行库和观察开关，退出时恢复。
    key="RAH_M030_TAIL_BUDGET"
    old=os.environ.get(key);path=os.environ["PATH"]
    try:
        os.environ["PATH"]=str(BIN)+os.pathsep+path
        if kind=="on":os.environ[key]="1"
        else:os.environ.pop(key,None)
        result=observe([str(exe),"case.cfg","--threads",str(threads)],folder,locks,150,64*1024**2)
    finally:
        os.environ["PATH"]=path
        if old is None:os.environ.pop(key,None)
        else:os.environ[key]=old
    print(kind,result["debuggerExitCode"],result["elapsedSeconds"],flush=True)
    terminal(folder)
    products=["restart_00001.csv","surface_00001.vtk","history_00001.csv"]
    traces=sorted(folder.glob("[0-9]*_*.csv"))
    expected={f"{i}_{s}.csv" for i in (0,1,249) for s in
              ("start","transport","afterViscousCall","source","time","boundary","rhs","edges")}
    assert {p.name for p in traces}==(expected if kind=="on" else set())
    save(folder/"replay.json",dict(kind=kind,threads=threads,physicsQualified=False,
        assemblyIterations=[0,1,249] if traces else [],parentSeedSHA256=prior["seedSHA256"],
        outputs={p.name:sha(p) for p in [*(folder/n for n in products),*traces]}))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path);parser.add_argument("build",type=Path);parser.add_argument("kind",choices=KINDS)
    parser.add_argument("--threads",type=int,choices=(1,2),default=2)
    a=parser.parse_args();replay(a.output,a.build,a.kind,a.threads)
