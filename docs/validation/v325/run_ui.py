"""从已发布的离线前端启动 MATLAB 工作台，验证双支管网格切换与实际重建。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe,sha


if __name__=="__main__":
    assert sys.argv[1].isidentifier()
    out=ROOT/"outputs/v325_validation"/sys.argv[1];out.mkdir(parents=True,exist_ok=False);front=out/"frontend"
    keep=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v324",*keep],cwd=ROOT))) as z:
        z.extractall(front)
    def q(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_fluid_ui({q(out)})"
    paths=[*front.rglob("*"),*ROOT.glob("*ShortSupport*.m"),ROOT/"rocketAeroHeatWorkbench.m",HERE/"PLAN.md",HERE/"check_fluid_ui.m",
           *ROOT.glob("scripts/shortSupport*.py"),*(ROOT/"cases/short_sphere").glob("*")]
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
              {str(p):sha(p) for p in paths if p.is_file()},180,100*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"]) and (out/"ui.json").is_file()
