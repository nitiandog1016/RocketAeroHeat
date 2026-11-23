"""从同一侧栏检验重剖民用试件，冻结独立发布的 WebGL 及本轮求解入口。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe,sha

if __name__=="__main__":
    name=sys.argv[1];assert name.isidentifier()
    out=ROOT/"outputs/v320_validation"/name;out.mkdir(exist_ok=False);front=out/"frontend"
    files=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v319",*files],cwd=ROOT))) as z:z.extractall(front)
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_ui({quote(out)})"
    files=[*front.rglob("*"),*ROOT.glob("*Thermistor*.m"),ROOT/"scripts/thermistorSubfaces.py",
           ROOT/"buildAerothermalViewerData.m",ROOT/"exportWallPairEvidence.m",HERE/"check_ui.m",HERE/"PLAN.md",Path(__file__),
           *(ROOT/"cases/thermistor_b/wallNodes96").glob("*"),*(out.parent/"fixture").glob("*")]
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(p):sha(p) for p in files if p.is_file()},420,300*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["survivingProcessIDs"] or r["supervisorError"]) and (out/"ui.json").is_file()
