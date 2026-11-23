"""接回侧栏真实创建的任务，从首个暂停点准备原生输入并在同态控制步前再次保存。"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha,write_json


def check(ui,original,prepared,out):
    work=Path(read(ui/"ui.json")["pausedNewJob"])
    first=read(work/"progress.json");assert first["status"]=="paused" and first["stage"]=="uniform"
    base=Path(read(work/"wire_job.json")["externalSeed"])
    # 冷启动目录保留计算体网；侧栏用的建网合同来自同内容的独立准备目录。
    from thermistorMeasuredCold import GRID_KEY
    assert sha(prepared/"mesh.su2")==read(GRID_KEY)["meshSHA256"]==sha(base/"cold/mesh.su2")
    token=str(uuid.uuid4());cmd=[sys.executable,str(ROOT/"scripts/thermistorWireWorkbench.py"),"resumeWire",str(work),str(prepared),"--invocation",token]
    proc=subprocess.Popen(cmd,cwd=ROOT,creationflags=subprocess.CREATE_NO_WINDOW if os.name=="nt" else 0)
    flag=work.with_name(work.name+"."+token+".pause")
    started=time.monotonic();sent=False
    while proc.poll() is None:
        state=read(work/"progress.json")
        if state["invocation"]==token and state["stage"]=="driver" and not sent:flag.touch(exist_ok=False);sent=True
        if time.monotonic()-started>360 and not sent:flag.touch(exist_ok=False);sent=True
        time.sleep(.1)
    receipt=read(work/("receipt-"+token+".json"))
    assert proc.returncode==0 and sent and receipt["status"]=="paused" and receipt["stage"]=="uniform_control"
    assert not (work/"run/uniform_control").exists() and not (work/"run/completed.json").exists()
    paths=("profile.dat","release.dat","initial_flow.csv")
    checks={n:sha(work/"run/prepared"/n) for n in paths}
    assert all(checks[n]==sha(original/"prepared"/n) for n in paths)
    write_json(out,dict(first=first,resume=receipt,baseFrame=read(base/"completed.json")["frameSHA256"],
               inputs=checks,built=read(work/"driver/build.json"),elapsedSeconds=time.monotonic()-started))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for n in ("ui","original","prepared","out"):p.add_argument(n,type=Path)
    a=p.parse_args();check(a.ui.resolve(),a.original.resolve(),a.prepared.resolve(),a.out.resolve())
