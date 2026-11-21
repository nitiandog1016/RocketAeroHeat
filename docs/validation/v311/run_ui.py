"""以发布的前端验证热敏珠全周侧吹侧栏和周向截面。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe, sha


def run(name):
    out=ROOT/"outputs/v311_validation"/name
    out.mkdir(parents=True,exist_ok=False)
    front=out/"frontend"
    snapshot=subprocess.check_output(["git","archive","--format=zip","v310","web","createAerothermalWebTab.m",
                                      "createAerothermalGeometryView.m","exportAerothermalWeb.m"],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(snapshot)) as z:
        z.extractall(front)
    exporter=ROOT/"createAerothermalEvidenceView.m"
    shutil.copyfile(exporter,front/exporter.name)
    def quote(p):
        return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_sidewind_ui({quote(out)})"
    files=[*front.rglob("*"),exporter,HERE/"WORKBENCH.md",HERE/"check_sidewind_ui.m",Path(__file__).resolve(),
           ROOT/"cases/thermistor_b/sidewind3.json",ROOT/"cases/thermistor_b/sidewind3.zip"]
    names=("rocketAeroHeatWorkbench","thermistorWorkbenchBoundary","aeroHeatWorkflowState","prepareThermistorLeadMesh",
           "runAeroHeatSingleCase","runThermistorSidewindWorkbench","buildAerothermalViewerData","openThermistorArchive")
    files.extend(ROOT/(n+".m") for n in names)
    files.extend(ROOT/"scripts"/(n+".py") for n in ("thermistorSidewindWorkbench","thermistorLeadCrossflow","thermistorCrossflowSurface"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(p):sha(p) for p in files if p.is_file()},180,180*1024**2)
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
    assert all(rec["inputsStable"].values()) and (out/"ui.json").is_file()
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<180*1024**2


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    args=parser.parse_args()
    if not args.name.isidentifier():parser.error("请使用新的目录名")
    run(args.name)
