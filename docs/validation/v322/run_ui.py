"""只使用发布前端检查银球真实尾迹计算、快照身份及导出。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe,sha

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("name");p.add_argument("kind",choices=("compute","read"));a=p.parse_args()
    assert a.name.isidentifier();out=ROOT/"outputs/v322_validation"/a.name;out.mkdir(exist_ok=False);front=out/"frontend"
    fixed=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v321",*fixed],cwd=ROOT))) as z:z.extractall(front)
    def q(f):return "'"+f.as_posix().replace("'","''")+"'"
    cmd=f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_ui({q(out)},'{a.kind}')"
    files=[*front.rglob("*"),*ROOT.glob("*ShortSphere*.m"),ROOT/"rocketAeroHeatWorkbench.m",ROOT/"buildAerothermalViewerData.m",
           ROOT/"exportWallPairEvidence.m",*ROOT.glob("scripts/*hortSphere*.py"),*(ROOT/"cases/short_sphere").glob("*"),HERE/"check_ui.m",HERE/"PLAN.md",Path(__file__)]
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(f):sha(f) for f in files if f.is_file()},240,200*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values()) and (out/"ui.json").is_file()
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
