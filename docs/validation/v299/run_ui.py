"""在发布用 Web 资源上验证生产工作台，保留并发的前端修改。"""
import argparse
import io
import shutil
import subprocess
import zipfile

from check_grid import HERE, ROOT, sha
from runSU2NativeDebugger import observe


def run(name, wired=False):
    out = ROOT/"outputs/v299_validation"/name
    out.mkdir(exist_ok=False)
    front = out/"frontend"
    files = ["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    snapshot = subprocess.check_output(["git","archive","--format=zip","v298",*files], cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(snapshot)) as z:
        z.extractall(front)
    script = HERE.parent/"v298/check_pair_ui.m" if wired else HERE/"check_midbody_ui.m"
    def quote(value):
        return "'"+value.as_posix().replace("'", "''")+"'"
    call = f"check_pair_ui({quote(out)},'wired')" if wired else f"check_midbody_ui({quote(out)})"
    matlab = f"cd({quote(out)});addpath({quote(ROOT)},{quote(script.parent)});addpath({quote(front)},'-begin');{call}"
    sources = [script,HERE/"PLAN.md",*front.rglob("*")]
    sources.extend(ROOT/name for name in ("rocketAeroHeatWorkbench.m","aeroHeatWorkflowState.m","runAeroHeatSingleCase.m",
         "runThermistorLeadWorkbench.m","attachThermistorLeadNearWall.m","attachThermistorPairedControl.m",
         "buildAerothermalViewerData.m","prepareThermistorLeadMesh.m","openThermistorArchive.m","thermistorWorkbenchModel.m"))
    for case in ("midbody3","paired3","leads3"):
        sources.extend((ROOT/"cases/thermistor_b").glob(case+".*"))
    rec = observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",matlab],out,
                  {str(p):sha(p) for p in sources if p.is_file()},180,200*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and not rec["supervisorError"],out/"debugger.log"
    assert all(rec["inputsStable"].values()) and (out/"ui.json").is_file()


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("name", help="本轮 outputs 内的新目录名")
    cli.add_argument("--wired",action="store_true",help="回归检查带线试件")
    a=cli.parse_args()
    if not a.name.isidentifier():
        cli.error("请使用单个目录名")
    run(a.name,a.wired)
