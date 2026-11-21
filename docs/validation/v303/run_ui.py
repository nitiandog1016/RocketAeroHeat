"""在发布版前端上执行两种热边界的实际侧栏、WebGL 与导出检查。"""
import argparse
import io
import shutil
import subprocess
import zipfile

from check_state import HERE,NEW,ROOT,sha
from runSU2NativeDebugger import observe


def run(name):
    out=NEW.parent/name;out.mkdir(exist_ok=False);front=out/"frontend"
    assets=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    snapshot=subprocess.check_output(["git","archive","--format=zip","v302",*assets],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(snapshot)) as z:z.extractall(front)
    shutil.copyfile(ROOT/"createAerothermalEvidenceView.m",front/"createAerothermalEvidenceView.m")
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_domain_ui({quote(out)})"
    files=[HERE/"check_domain_ui.m",HERE/"PLAN.md",*front.rglob("*")]
    files.extend(ROOT/n for n in ("rocketAeroHeatWorkbench.m","aeroHeatWorkflowState.m","runAeroHeatSingleCase.m","attachThermistorOuterDomain.m",
        "buildAerothermalViewerData.m","openThermistorArchive.m","scaleThermistorWallField.m","runThermistorLeadWorkbench.m"))
    for key in ("outer30","leadsMidbody3","paired3","leads3","midbody3"):files.extend((ROOT/"cases/thermistor_b").glob(key+".*"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(p):sha(p) for p in files if p.is_file()},180,200*1024**2)
    assert rec["debuggerExitCode"]==0 and not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]),out/"debugger.log"
    assert all(rec["inputsStable"].values()) and (out/"ui.json").is_file()


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("directory");a=p.parse_args()
    if not a.directory.isidentifier():p.error("请使用单个新目录名")
    run(a.directory)
