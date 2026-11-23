"""隔离并发网页，分别从已发布前端检查两个物理时间步分支的真实 MATLAB 入口。"""
import io
import shutil
import subprocess
import sys
from zipfile import ZipFile
from build_time import ROOT,HERE,sha
from runSU2NativeDebugger import observe

if __name__=="__main__":
    assert len(sys.argv)==3 and sys.argv[1].isidentifier() and sys.argv[2] in ("full","half")
    branch=sys.argv[2];out=ROOT/"outputs/v333_validation"/sys.argv[1]/branch;out.mkdir(parents=True,exist_ok=False)
    front=out/"frontend";keep=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v332",*keep],cwd=ROOT))) as z:z.extractall(front)
    def q(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_time_ui({q(out)},'{branch}')"
    files=[*front.rglob("*"),*ROOT.glob("*ShortSupport*.m"),ROOT/"shortSupportDirectionReasons.m",ROOT/"rocketAeroHeatWorkbench.m",
           ROOT/"buildAerothermalViewerData.m",HERE/"check_time_ui.m",*(ROOT/"cases/short_sphere").glob("supportedTime65*")]
    proc=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(p):sha(p) for p in files if p.is_file()},180,350*1024**2)
    assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
    assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"]) and (out/"ui.json").is_file()
