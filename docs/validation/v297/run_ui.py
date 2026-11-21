"""运行带引线试件的实际五步工作台，以已发布的网页资源隔离并发前端修改。"""
from pathlib import Path
import io
import shutil
import subprocess
import sys
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe


def run(out,mode,other=None):
    out.relative_to(ROOT/"outputs/v297_validation")
    out.mkdir(exist_ok=False)
    front=out/"frontend"
    resources=["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    data=subprocess.check_output(["git","archive","--format=zip","v296",*resources],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(data)) as z:z.extractall(front)
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    extra=","+quote(other) if other else ""
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_ui({quote(out)},'{mode}'{extra})"
    sources=[HERE/"check_ui.m",Path(__file__).resolve(),HERE/"PLAN.md",*front.rglob("*")]
    sources.extend(ROOT/p for p in ("rocketAeroHeatWorkbench.m","prepareThermistorLeadMesh.m","runThermistorLeadWorkbench.m",
                   "buildAerothermalViewerData.m","aeroHeatWorkflowState.m","runAeroHeatSingleCase.m","scripts/thermistorLeadSurface.py",
                   "scripts/thermistorLeadWorkbench.py","scripts/runThermistorLeadCase.py","scripts/thermistorLeadMesh.py",
                   "cases/thermistor_b/leads3.zip","cases/thermistor_b/leads3.json","attachThermistorLeadNearWall.m",
                   "scripts/thermistorLeadNearWall.py","docs/validation/v282/observe_bead.py"))
    if other:sources.append(other)
    seconds=180 if mode=="archive" else 180+180+720+3*60+120+60+2*60
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(p):sha(p) for p in sources if p.is_file()},seconds,80*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and (out/"ui.json").exists(),out/"debugger.log"


if __name__=="__main__":
    assert sys.argv[2] in ("archive","solve")
    run(Path(sys.argv[1]).resolve(),sys.argv[2],Path(sys.argv[3]).resolve() if len(sys.argv)>3 else None)
