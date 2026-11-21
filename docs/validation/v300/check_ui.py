"""在已发布 Web 资源上检查生产 MATLAB 工作台；不改动工作区的并发前端。"""
import argparse
import io
import shutil
import subprocess
import zipfile

from mesh_audit import HERE, NEW, ROOT, sha
from runSU2NativeDebugger import observe


def check(out):
    out.mkdir(exist_ok=False)
    front=out/"frontend"
    resources=["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    archived=subprocess.check_output(["git","archive","--format=zip","v299",*resources],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(archived)) as z:z.extractall(front)
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_wired_ui({quote(out)})"
    files=[HERE/"check_wired_ui.m",HERE/"PLAN.md",*front.rglob("*")]
    files.extend(ROOT/n for n in ("rocketAeroHeatWorkbench.m","aeroHeatWorkflowState.m","runAeroHeatSingleCase.m",
        "runThermistorLeadWorkbench.m","attachThermistorLeadNearWall.m","attachThermistorLeadGridComparison.m",
        "attachThermistorPairedControl.m","buildAerothermalViewerData.m","prepareThermistorLeadMesh.m","openThermistorArchive.m"))
    for label in ("leads3","paired3","midbody3","leadsMidbody3"):
        files.extend((ROOT/"cases/thermistor_b").glob(label+".*"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(p):sha(p) for p in files if p.is_file()},180,200*1024**2)
    assert rec["debuggerExitCode"]==0 and not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]),out/"debugger.log"
    assert all(rec["inputsStable"].values()) and (out/"ui.json").is_file()


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("directory")
    name=parser.parse_args().directory
    if not name.isidentifier():parser.error("请使用单个目录名")
    check(NEW.parent/name)
