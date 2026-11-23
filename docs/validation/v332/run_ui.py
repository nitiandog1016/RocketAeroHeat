"""隔离并发网页，用 v331 已发布的前端测试双支管周向加密的真实侧栏入口。"""
import io
import shutil
import subprocess
import sys
from zipfile import ZipFile
from sector_grid import ROOT, HERE, sha
from runSU2NativeDebugger import observe

if __name__=="__main__":
    assert len(sys.argv)==2 and sys.argv[1].isidentifier()
    out=ROOT/"outputs/v332_validation"/sys.argv[1];out.mkdir(exist_ok=False)
    front=out/"frontend";keep=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v331",*keep],cwd=ROOT))) as z:z.extractall(front)
    def q(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({q(out)});addpath({q(ROOT)},{q(HERE)});addpath({q(front)},'-begin');check_sector_ui({q(out)})"
    files=[*front.rglob("*"),*ROOT.glob("*ShortSupport*.m"),ROOT/"rocketAeroHeatWorkbench.m",ROOT/"buildAerothermalViewerData.m",
           HERE/"check_sector_ui.m",*(ROOT/"cases/short_sphere").glob("supportedAngular*")]
    proc=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                 {str(p):sha(p) for p in files if p.is_file()},180,350*1024**2)
    assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
    assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"]) and (out/"ui.json").is_file()
