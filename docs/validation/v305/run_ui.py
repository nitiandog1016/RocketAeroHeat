"""以已发布 Web 前端和本轮 MATLAB 文件核验热敏珠分带页面。"""
import argparse
import io
import shutil
import subprocess
import zipfile

from common_layer import BASE, HERE, ROOT, sha
from runSU2NativeDebugger import observe


def run(name):
    out=BASE/name;out.mkdir(exist_ok=False);front=out/"frontend"
    assets=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    snapshot=subprocess.check_output(["git","archive","--format=zip","v304",*assets],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(snapshot)) as z:z.extractall(front)
    shutil.copyfile(ROOT/"createAerothermalEvidenceView.m",front/"createAerothermalEvidenceView.m")
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_budget_ui({quote(out)})"
    files=[HERE/"check_budget_ui.m",HERE/"PLAN.md",*front.rglob("*"),ROOT/"cases/thermistor_b/layerBands.json"]
    files.extend(ROOT/n for n in ("rocketAeroHeatWorkbench.m","aeroHeatWorkflowState.m","runAeroHeatSingleCase.m",
        "attachThermistorLayerBudget.m","attachThermistorThreeGrids.m","runThermistorLeadWorkbench.m",
        "prepareThermistorLeadMesh.m","buildAerothermalViewerData.m","openThermistorArchive.m","createAerothermalEvidenceView.m"))
    for key in ("leads3","leadsMidbody3","leadsThird3"):files.extend((ROOT/"cases/thermistor_b").glob(key+".*"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(p):sha(p) for p in files if p.is_file()},210,200*1024**2)
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]) and (out/"ui.json").is_file()


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("directory");a=p.parse_args()
    if not a.directory.isidentifier():p.error("请填写一个新目录名")
    run(a.directory)
