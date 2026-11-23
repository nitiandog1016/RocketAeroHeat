"""在一个工作台目录管理引线算例、内部均温初猜及可复用的原生驱动。"""
import argparse
from pathlib import Path
import uuid

from buildThermistorWireDriver import build
from runEmbeddedSU2 import atomic_write_json
from runThermistorMeasuredCold import SavedPause, exclusive, execute
from runThermistorWireCase import run
from thermistorBead import ROOT, read, sha, write_json
from thermistorMeasuredCold import grid


def calculate(action, work, prepared, base=None, invocation=None):
    grid(prepared)
    token=str(uuid.UUID(invocation)) if invocation else str(uuid.uuid4())
    flag=work.with_name(work.name+"."+token+".pause")
    if action=="solve":work.mkdir(parents=True,exist_ok=False)
    with exclusive(work):
        path=work/"wire_job.json"
        if action=="solve":
            write_json(path,dict(kind="B-wire-3ms-12K",directory=str(work),externalSeed=str(base) if base else "",
                                 meshSHA256=sha(prepared/"mesh.su2"),caseSHA256=sha(prepared/"case.json")))
        job=read(path)
        assert job["kind"]=="B-wire-3ms-12K" and job["directory"]==str(work)
        assert job["meshSHA256"]==sha(prepared/"mesh.su2") and job["caseSHA256"]==sha(prepared/"case.json")
        progress=dict(invocation=token,stage="uniform",status="running")
        def checkpoint(stage):
            progress["stage"]=stage
            if flag.exists():raise SavedPause
            atomic_write_json(work/"progress.json",progress)
        try:
            checkpoint("uniform")
            if job["externalSeed"]:
                base=Path(job["externalSeed"]).resolve();base.relative_to(ROOT/"outputs")
            else:
                base=work/"uniform"
                if not (base/"completed.json").exists():
                    execute("resume" if base.exists() else "cold",base,prepared,parent_pause=flag)
            checkpoint("driver")
            installed=work/"driver"
            if not installed.exists():build(installed)
            checked=read(installed/"build.json")
            assert checked["driverSourceSHA256"]==sha(ROOT/"scripts/buildThermistorWireDriver.py")
            for name,digest in checked["programs"].items():assert sha(installed/(name+".exe"))==digest
            result=work/"run"
            run(base,result,installed,resume=result.exists(),checkpoint=checkpoint)
            progress["status"]="complete"
        except SavedPause:
            progress["status"]="paused"
        atomic_write_json(work/"progress.json",progress)
        write_json(work/("receipt-"+token+".json"),progress)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("action",choices=("solve","resumeWire"))
    parser.add_argument("folder",type=Path);parser.add_argument("prepared",type=Path);parser.add_argument("--base",type=Path)
    parser.add_argument("--invocation")
    a=parser.parse_args();work,prepared=a.folder.resolve(),a.prepared.resolve()
    for p in (work,prepared):p.relative_to(ROOT/"outputs")
    base=a.base.resolve() if a.base else None
    if base:base.relative_to(ROOT/"outputs")
    calculate(a.action,work,prepared,base,a.invocation)
