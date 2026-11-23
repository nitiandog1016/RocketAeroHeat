"""在发布前端的同一工作台运行银球已算结果和真实冷启动。子阶段各有原生预算。"""
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe

if __name__=="__main__":
    name=sys.argv[1];assert name.isidentifier();out=ROOT/"outputs/v321_validation"/name;out.mkdir(exist_ok=False)
    front=out/"frontend";fixed=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v320",*fixed],cwd=ROOT))) as z:z.extractall(front)
    def q(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_ui({q(out)})"
    files=[*front.rglob("*"),*ROOT.glob("*ShortSphere*.m"),ROOT/"rocketAeroHeatWorkbench.m",ROOT/"buildAerothermalViewerData.m",
           ROOT/"exportWallPairEvidence.m",*ROOT.glob("scripts/shortSphere*.py"),*(ROOT/"cases/short_sphere").glob("*"),
           HERE/"check_ui.m",HERE/"PLAN.md",Path(__file__)]
    locks={str(p):sha(p) for p in files if p.is_file()}
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,locks,240,200*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"]) and (out/"ui.json").is_file()
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<200*1024**2
