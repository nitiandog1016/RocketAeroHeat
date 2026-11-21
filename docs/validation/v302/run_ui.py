"""用已发布的 Web 资源运行真实功率工作台，并复核共用换算对光滑珠体的兼容性。"""
import argparse
import io
import shutil
import subprocess
import zipfile

from half_case import HERE,ROOT,sha
from runSU2NativeDebugger import observe


def run(name,smooth):
    root=ROOT/"outputs/v302_validation";out=root/name;out.mkdir(exist_ok=False)
    front=out/"frontend"
    resources=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","createAerothermalEvidenceView.m","exportAerothermalWeb.m"]
    release=subprocess.check_output(["git","archive","--format=zip","v301",*resources],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(release)) as z:z.extractall(front)
    script=(HERE.parent/"v286/check_power.m") if smooth else HERE/"check_power_ui.m"
    def arg(p):return "'"+p.as_posix().replace("'","''")+"'"
    call=f"check_power({arg(out)},4)" if smooth else f"check_power_ui({arg(out)},{arg(root/'native01')})"
    cmd=f"cd({arg(out)});addpath({arg(ROOT)},{arg(script.parent)});addpath({arg(front)},'-begin');{call}"
    sources=[script,HERE/"PLAN.md",root/"native01/half.json",root/"native01/direct_wall.csv",*front.rglob("*")]
    sources.extend(ROOT/n for n in ("rocketAeroHeatWorkbench.m","runAeroHeatSingleCase.m","aeroHeatWorkflowState.m","thermistorNetPowerLimit.m",
        "thermistorWorkbenchBoundary.m","scaleThermistorWallField.m","runThermistorLeadWorkbench.m","runThermistorBeadCase.m",
        "buildAerothermalViewerData.m","attachThermistorLeadNearWall.m","openThermistorArchive.m"))
    sources.extend((ROOT/"cases/thermistor_b").glob("*.zip"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(p):sha(p) for p in sources if p.is_file()},180,200*1024**2)
    assert rec["debuggerExitCode"]==0 and not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"]),out/"debugger.log"
    assert all(rec["inputsStable"].values()) and (out/("result.json" if smooth else "ui.json")).is_file()


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("name");p.add_argument("--smooth",action="store_true");a=p.parse_args()
    if not a.name.isidentifier():p.error("请用一个新的目录名")
    run(a.name,a.smooth)
