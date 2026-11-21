"""在 v303 已发布的前端上运行第三档带线试件的生产侧栏和当前导出器。"""
import argparse
import io
import shutil
import subprocess
import zipfile

from audit_third_mesh import FINE, HERE, ROOT, sha
from runSU2NativeDebugger import observe


def run(name):
    out=FINE.parent/name
    out.mkdir(exist_ok=False)
    front=out/"frontend"
    assets=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    snapshot=subprocess.check_output(["git","archive","--format=zip","v303",*assets],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(snapshot)) as z:
        z.extractall(front)
    shutil.copyfile(ROOT/"createAerothermalEvidenceView.m",front/"createAerothermalEvidenceView.m")
    def quote(p):
        return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_three_ui({quote(out)})"
    files=[HERE/"check_three_ui.m",HERE/"PLAN.md",*front.rglob("*")]
    files.extend(ROOT/n for n in ("rocketAeroHeatWorkbench.m","aeroHeatWorkflowState.m","runAeroHeatSingleCase.m",
        "thermistorLeadGridChoice.m","attachThermistorThreeGrids.m","prepareThermistorLeadMesh.m","thermistorNetPowerLimit.m",
        "buildAerothermalViewerData.m","openThermistorArchive.m","scaleThermistorWallField.m","runThermistorLeadWorkbench.m",
        "scripts/thermistorLeadWorkbench.py","scripts/runThermistorLeadCase.py","scripts/thermistorMidbodyMesh.py","scripts/thermistorLeadMesh.py"))
    for key in ("leads3","leadsMidbody3","leadsThird3"):
        files.extend((ROOT/"cases/thermistor_b").glob(key+".*"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(p):sha(p) for p in files if p.is_file()},210,220*1024**2)
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]),out/"debugger.log"
    assert all(rec["inputsStable"].values()) and (out/"ui.json").is_file()


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("name");a=p.parse_args()
    if not a.name.isidentifier():p.error("请使用单个新目录名")
    run(a.name)
