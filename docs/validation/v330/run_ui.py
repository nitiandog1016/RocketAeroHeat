"""在已发布的离线前端运行实际近壁对照页面，避开并发编辑。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

from check_grid import HERE, ROOT
from runSU2NativeDebugger import observe, sha


if __name__=="__main__":
    assert len(sys.argv)==2 and sys.argv[1].isidentifier()
    out=ROOT/"outputs/v330_validation"/sys.argv[1];out.mkdir(exist_ok=False)
    front=out/"frontend";keep=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v329",*keep],cwd=ROOT))) as z:z.extractall(front)
    def q(path):return "'"+path.as_posix().replace("'","''")+"'"
    command=f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_ui({q(out)})"
    files=[*front.rglob("*"),*ROOT.glob("*ShortSupport*.m"),ROOT/"buildAerothermalViewerData.m",ROOT/"rocketAeroHeatWorkbench.m",
           Path(__file__),HERE/"check_ui.m",*(ROOT/"cases/short_sphere").glob("*")]
    proc=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",command],out,
                 {str(p):sha(p) for p in files if p.is_file()},180,350*1024**2)
    assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
    assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"]) and (out/"ui.json").is_file()
