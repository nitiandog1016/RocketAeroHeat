"""限制本轮只读后处理、复核与工作台验收的耗时和内存，记录实际终态。"""
import argparse
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from zipfile import ZipFile

import psutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT/"outputs/v331_validation"


def run(mode, name, data, grid):
    if not name.isidentifier():
        raise ValueError(name)
    out = BASE/name; out.mkdir(parents=True,exist_ok=False)
    ui = mode=="ui"
    if ui:
        front=out/"frontend"
        keep=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
        with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v330",*keep],cwd=ROOT))) as z:
            z.extractall(front)
        def q(p): return "'"+p.as_posix().replace("'","''")+"'"
        cmd = [shutil.which("matlab"),"-wait","-singleCompThread","-batch",
               f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_direction_ui({q(out)},'{grid}')"]
        product = out/"ui.json"
    elif mode=="read":
        cmd = [sys.executable,str(HERE/"direction.py"),str(out/"data")]; product = out/"data/supportedDirection65.json"
    else:
        cmd = [sys.executable,str(HERE/"verify_direction.py"),str(data),str(out/"verified.json")]
        product = out/"verified.json"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    seconds = 180 if ui else 120; memory = (4 if ui else 2)*1024**3
    budget = 350*1024**2
    start = time.monotonic(); peak = peak_disk = 0; stopped = None
    with (out/"run.log").open("xb") as log:
        p = psutil.Popen(cmd,cwd=out,env=env,stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
        family = {p.pid:p}
        try:
            while p.poll() is None:
                try: family.update({c.pid:c for c in p.children(recursive=True)})
                except psutil.NoSuchProcess: pass
                rss = 0
                for c in family.values():
                    try: rss += c.memory_info().rss
                    except psutil.NoSuchProcess: pass
                peak = max(peak,rss)
                size = sum(f.stat().st_size for f in BASE.rglob("*") if f.is_file())
                peak_disk = max(peak_disk,size)
                if time.monotonic()-start>seconds or rss>memory or size>budget:
                    stopped="time" if time.monotonic()-start>seconds else ("memory" if rss>memory else "disk")
                    break
                time.sleep(.10)
        finally:
            if p.poll() is None:
                subprocess.run(["taskkill","/PID",str(p.pid),"/T","/F"],stdout=log,stderr=log)
            for c in reversed(list(family.values())):
                try:
                    if c.is_running(): c.kill()
                except psutil.NoSuchProcess: pass
            _, alive = psutil.wait_procs(list(family.values()), timeout=5)
            code=p.wait(5)
    receipt=dict(mode=mode,command=cmd,exitCode=code,stopReason=stopped,survivors=[c.pid for c in alive],
                 seconds=time.monotonic()-start,peakRSS=peak,peakOutputBytes=peak_disk,
                 budgetSeconds=seconds,budgetRSS=memory,product=product.name,productPresent=product.is_file())
    (out/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf8")
    print(receipt)
    if code or stopped or alive or not product.is_file():
        raise SystemExit(1)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=("read","verify","ui"))
    parser.add_argument("name");parser.add_argument("--data",type=Path,default=ROOT/"cases/short_sphere")
    parser.add_argument("--grid",choices=("coarse","fine"),default="fine")
    a=parser.parse_args();run(a.mode,a.name,a.data.resolve(),a.grid)
