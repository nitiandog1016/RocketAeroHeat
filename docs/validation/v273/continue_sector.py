"""从已核对的完整 13000 帧接续同一个低速三维算例，独立记录预算和进程资源。"""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

import psutil

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from cm011Tunnel import write_json
from readCM011TunnelSector import state_at_frame
from runEmbeddedSU2 import sha


def main():
    base=ROOT/"outputs/v273_validation"
    parent,out=base/"sector01",base/"sector02"
    record=state_at_frame(parent)
    assert record["timedOut"] and record["savedIndex"]==13000
    receipt=parent/"engine_process.json"
    assert sha(receipt)=="1461b47ec49cdce1e6ec0ef169c561e9c0a96e901d594b47254c345b2eb9570a"
    initial=parent/"restart_013000.csv"
    assert sha(initial)=="38acf6c69405b7c81201e4905bf90a1a4d45c31ac54c71a5b6638e4375706c75"
    out.mkdir(exist_ok=False)
    for n in ("geometry.json","nodes.npz","mesh.su2"): shutil.copyfile(parent/n,out/n)
    shutil.copyfile(initial,out/"initial_flow.csv")
    cfg=(parent/"case.cfg").read_text("ascii")
    for a,b in (("ITER= 16000\n","ITER= 5000\n"),("(1000, 1000, 16000)","(1000, 1000, 5000)")):
        assert cfg.count(a)==1
        cfg=cfg.replace(a,b)
    (out/"case.cfg").write_text(cfg,encoding="ascii")
    prior=json.loads(receipt.read_bytes())["request"]
    extra=[receipt,initial,out/"geometry.json",out/"nodes.npz",Path(__file__)]
    request=dict(prior,workDirectory=str(out),maximumSeconds=320,jobToken=uuid.uuid4().hex,
                 sourceFiles=prior["sourceFiles"]+[p.resolve().relative_to(ROOT).as_posix() for p in extra])
    write_json(out/"engine_request.json",request)
    monitor=base/"monitor02"
    monitor.mkdir(exist_ok=False)
    start=time.monotonic()
    peak,size_peak,seen,reason=0,0,{},None
    with (monitor/"driver.log").open("xb") as log:
        proc=subprocess.Popen([sys.executable,str(ROOT/"scripts/runEmbeddedSU2.py"),"--request",str(out/"engine_request.json")],
                              stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
        root=psutil.Process(proc.pid)
        while proc.poll() is None:
            try: family=[root]+root.children(recursive=True)
            except psutil.NoSuchProcess: family=[]
            mem=0
            for child in family:
                try:
                    mem+=child.memory_info().rss
                    seen[child.pid]=child
                except psutil.NoSuchProcess: pass
            size=sum(p.stat().st_size for p in out.iterdir() if p.is_file())
            peak,size_peak=max(peak,mem),max(size_peak,size)
            if reason is None and (peak>2*1024**3 or size_peak>64*1024**2):
                reason="memory" if peak>2*1024**3 else "disk"
                write_json(out/"engine_cancel.json",dict(jobToken=request["jobToken"]))
            time.sleep(.25)
    _,alive=psutil.wait_procs(list(seen.values()),timeout=3)
    write_json(monitor/"resources.json",dict(elapsedSeconds=time.monotonic()-start,exitCode=proc.returncode,
               peakRSSBytes=peak,peakDirectoryBytes=size_peak,limitReason=reason,observedPids=list(seen),livePids=[x.pid for x in alive]))
    return proc.returncode


if __name__=="__main__":
    sys.exit(main())
