"""用发布前端从统一侧栏核对银球的新增剖面页。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe, sha


if __name__ == "__main__":
    assert sys.argv[1].isidentifier()
    out = ROOT/"outputs/v323_validation"/sys.argv[1]
    out.mkdir(parents=True,exist_ok=False)
    front = out/"frontend"
    keep = ["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v322",*keep],cwd=ROOT))) as z:
        z.extractall(front)
    def quoted(path):
        return "'"+path.as_posix().replace("'","''")+"'"
    cmd = f"cd({quoted(out)});addpath({quoted(ROOT)},{quoted(HERE)});addpath({quoted(front)},'-begin');check_resolution_ui({quoted(out)})"
    paths = [*front.rglob("*"),*ROOT.glob("*ShortSphere*.m"),ROOT/"rocketAeroHeatWorkbench.m",
             *(ROOT/"cases/short_sphere").glob("*"),HERE/"check_resolution_ui.m",HERE/"PLAN.md"]
    r = observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(p):sha(p) for p in paths if p.is_file()},180,100*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not (r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"]) and (out/"ui.json").is_file()
