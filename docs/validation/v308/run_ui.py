"""用已发布的离线前端、当前 MATLAB 代码和冻结参照包运行指定网格的侧栏检查。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from runSU2NativeDebugger import observe, sha


def run(level, name):
    out = ROOT / "outputs/v308_validation" / name
    out.mkdir(parents=True, exist_ok=False)
    front = out / "frontend"
    assets = ["web", "createAerothermalWebTab.m", "createAerothermalGeometryView.m", "exportAerothermalWeb.m"]
    snapshot = subprocess.check_output(["git", "archive", "--format=zip", "v307", *assets], cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(snapshot)) as z:
        z.extractall(front)
    shutil.copyfile(ROOT/"createAerothermalEvidenceView.m", front/"createAerothermalEvidenceView.m")
    def quote(path):
        return "'" + path.as_posix().replace("'", "''") + "'"
    cmd = f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_wall_ui({quote(out)},{level})"
    inputs = [*front.rglob("*"), HERE/"check_wall_ui.m", HERE/"PLAN.md", Path(__file__).resolve(),
              ROOT/"cases/thermistor_b/layerBands.json", ROOT/"createAerothermalEvidenceView.m"]
    names = ["aeroHeatWorkflowState", "runAeroHeatSingleCase", "rocketAeroHeatWorkbench", "buildAerothermalViewerData",
             "runThermistorLeadWorkbench", "thermistorWorkbenchBoundary", "scaleThermistorWallField",
             "attachThermistorLeadGridComparison", "attachThermistorLayerBudget", "attachThermistorThreeGrids",
             "attachThermistorOuterDomain", "attachThermistorPairedControl", "prepareThermistorLeadMesh", "openThermistorArchive",
             "exportWallPairEvidence"]
    inputs.extend(ROOT/(k+".m") for k in names)
    for key in ("leads3", "leadsMidbody3", "leadsThird3", "midbody3", "outer30", "paired3"):
        inputs.extend((ROOT/"cases/thermistor_b").glob(key+".*"))
    locks = {str(p): sha(p) for p in inputs if p.is_file()}
    rec = observe([shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", cmd], out, locks, 300, 450*1024**2)
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
    assert (out/"ui.json").is_file() and sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<450*1024**2


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("level", type=int, choices=[1,2]);parser.add_argument("name")
    args=parser.parse_args()
    if not args.name.isidentifier():parser.error("请使用未占用的目录名")
    run(args.level,args.name)
