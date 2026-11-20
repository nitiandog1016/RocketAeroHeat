"""按实际验收结果串行运行民用热敏珠近壁对照，也可接管已经启动的第一半步。"""
import argparse
from pathlib import Path
import subprocess
import sys
import psutil

from advance import HERE,RUN,OFFICIAL,read


def await_first(area):
    folder=area/"first01";path=folder/"debugger_process.json";r=read(path)
    assert r["command"]==[str(OFFICIAL),"case.cfg","--threads","1"]
    if r["status"]=="running":
        parent=psutil.Process(r["supervisorProcessID"])
        assert any(Path(s).name=="advance.py" for s in parent.cmdline()),"监督进程身份不符"
        child=psutil.Process(r["processID"])
        assert child.ppid()==parent.pid and Path(child.exe())==OFFICIAL
        assert Path(child.cwd())==folder
        while parent.is_running():
            try:parent.wait(timeout=20)
            except psutil.TimeoutExpired:continue
        # 读取结束回执。进程退出自身不足以放行后处理。
    r=read(path)
    assert r["status"]=="terminal" and not r["stopReason"] and r["debuggerExitCode"]==0 and not r["survivingProcessIDs"]


def run(area,existing):
    def call(script,*args):
        cmd=[sys.executable,str(HERE/script),*map(str,args),"--root",str(area)]
        print(" ".join(cmd),flush=True);subprocess.run(cmd,check=True)
    if existing:await_first(area)
    else:call("advance.py",1)
    for i in (1,2):
        if i==2:call("advance.py",i)
        for stage in ("state","observe","release","energy"):call("verify.py",stage,i)
    call("compare.py")


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument("--root",type=Path,default=RUN)
    cli.add_argument("--existing-first",action="store_true")
    a=cli.parse_args();area=a.root.resolve();area.relative_to(RUN);run(area,a.existing_first)
