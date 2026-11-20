"""在独立目录核对热敏珠净功率边界，用当前发布的 Web 资源验收。"""
import argparse
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe


def run(cells, name):
    here = Path(__file__).resolve().parent
    out = ROOT / "outputs" / "v286_validation" / name
    out.mkdir(parents=True, exist_ok=False)
    # 使用已发布前端。并行工作区内尚未发布的网页不进入本轮回执。
    release = out / "release"
    resources = ["web", "createAerothermalWebTab.m", "createAerothermalEvidenceView.m", "exportAerothermalWeb.m"]
    archive = subprocess.check_output(["git", "archive", "--format=zip", "v285", *resources], cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(archive)) as files:
        files.extractall(release)
    inputs = [here / "PLAN.md", here / "check_power.m", Path(__file__).resolve(),
              ROOT / "docs/validation/v285/evidence/native.json", *release.rglob("*"),
              *(ROOT / "cases/thermistor_b").glob("*.zip")]
    inputs += [ROOT / (name + ".m") for name in (
        "thermistorNetPowerLimit", "thermistorWorkbenchBoundary", "openThermistorArchive",
        "runThermistorBeadCase", "runAeroHeatSingleCase", "rocketAeroHeatWorkbench",
        "aeroHeatWorkflowState", "buildAerothermalViewerData")]
    locks = {str(p): sha(p) for p in inputs if p.is_file()}
    def arg(path):
        return "'" + path.as_posix().replace("'", "''") + "'"
    call = (f"cd({arg(out)});addpath({arg(ROOT)},{arg(here)});addpath({arg(release)},'-begin');"
            f"check_power({arg(out)},{cells})")
    record = observe([shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", call],
                     out, locks, 180, 80 * 1024**2)
    assert record["debuggerExitCode"] == 0 and record["stopReason"] is None, out / "debugger.log"
    proof = json.loads((out / "result.json").read_text("utf8"))
    assert proof["cellAngleDeg"] == 30 / cells and proof["workbench"]["faces"] == 5216
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) < 80 * 1024**2
    print(json.dumps({"seconds": record["elapsedSeconds"], "cases": len(proof["cases"]), "output": str(out)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cells", type=int, choices=[2, 4])
    parser.add_argument("directory")
    args = parser.parse_args()
    if Path(args.directory).name != args.directory or args.directory in {".", ".."}:
        parser.error("需要新的单层输出目录名")
    run(args.cells, args.directory)
