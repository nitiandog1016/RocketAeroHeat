"""从隔离安装的 MATLAB 工程入口执行新算或在新窗口重开。"""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
from zipfile import ZipFile

from clean_install import ROOT,HERE,AREA
sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import invoke


def run(app,name,job=None):
    app.relative_to(AREA);assert name.isidentifier()
    proof=AREA/"windows";proof.mkdir(exist_ok=True);out=proof/name;out.mkdir()
    with ZipFile(app/"cases/short_sphere/supported65.zip") as z:stl=z.read("supported_m.stl")
    fixture=proof/"supported_m.stl"
    if fixture.exists():assert fixture.read_bytes()==stl
    else:fixture.write_bytes(stl)
    # 仅使用本安装的 MATLAB 路径与本次命令。原生进程分别由各自监督器核实。
    matlab=Path(shutil.which("matlab"));mode="reopen" if job else "solve"
    for k in ("MATLABPATH","RAH_THERMISTOR_SU2_BUILD","PYTHONPATH"):os.environ.pop(k,None)
    def quote(v):return "'"+str(v).replace("\\","/").replace("'","''")+"'"
    cmd=(f"restoredefaultpath;cd({quote(app)});addpath({quote(app/'docs/validation/v336')});"
         f"assert(strcmp(which('rocketAeroHeatWorkbench'),fullfile(pwd,'rocketAeroHeatWorkbench.m')));"
         f"check_clean_ui({quote(out)},'{mode}',{quote(job or '')})")
    invoke("window",[matlab,"-wait","-singleCompThread","-batch",cmd],app,out,6000 if not job else 240,6,bin_directory=matlab.parent)
    rec=json.loads((out/"window_process.json").read_bytes());ui=json.loads((out/"ui.json").read_bytes())
    assert rec["exitCode"]==0 and not(rec["livePids"] or rec["stopReason"]) and ui["action"]==mode
    origin=Path(ui["directory"]).resolve();origin.relative_to(app/"outputs")
    assert ui["faces"]==32480 and (origin/"completed.json").is_file()
    print(json.dumps(dict(action=mode,job=str(origin),frame=ui["source"]["sha256"],power=ui["outW"])))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("app",type=Path);p.add_argument("name");p.add_argument("--job",type=Path)
    a=p.parse_args();run(a.app.resolve(),a.name,a.job.resolve() if a.job else None)
