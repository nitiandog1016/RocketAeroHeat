"""从正常五步工作台验证民用热敏珠的读取或新算，前端固定在已发布版本。"""
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
from runSU2NativeDebugger import observe
from runEmbeddedSU2 import sha


def quoted(path):
    return "'" + path.as_posix().replace("'", "''") + "'"


def run(out, mode, job=None):
    out.relative_to(ROOT / "outputs/v295_validation")
    out.mkdir(exist_ok=False)
    front = out / "frontend"
    resources = ["web", "createAerothermalWebTab.m", "createAerothermalEvidenceView.m",
                 "createAerothermalGeometryView.m", "exportAerothermalWeb.m"]
    data = subprocess.check_output(["git", "archive", "--format=zip", "v294", *resources], cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        z.extractall(front)
    old = ROOT / "outputs/v294_validation" if job is None else job.resolve()
    command = (f"cd({quoted(out)});addpath({quoted(ROOT)},{quoted(HERE)});addpath({quoted(front)},'-begin');"
               f"check_workbench({quoted(out)},'{mode}',{quoted(old)})")
    files = [HERE / "PLAN.md", HERE / "check_workbench.m", Path(__file__).resolve(),
             ROOT / "scripts/runThermistorStartupJob.py", *front.rglob("*")]
    files.extend(ROOT / (name + ".m") for name in (
        "rocketAeroHeatWorkbench", "runThermistorStartupWorkbench", "prepareThermistorWorkbenchMesh",
        "verifyThermistorWorkbenchMesh", "thermistorWorkbenchBoundary", "readThermistorStartup",
        "runAeroHeatSingleCase", "aeroHeatWorkflowState", "runCM011Python"))
    seconds = 180 if mode == "read" else 180 + 2*(2*1100 + 100 + 2*180) + 2*180
    rec = observe([shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", command],
                  out, {str(p): sha(p) for p in files if p.is_file()}, seconds, 80*1024**2)
    assert rec["debuggerExitCode"] == 0 and not rec["stopReason"] and (out / "workflow.json").is_file(), out / "debugger.log"


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("out", type=Path)
    cli.add_argument("mode", choices=("read", "run", "resume"))
    cli.add_argument("--job", type=Path)
    a = cli.parse_args()
    if a.mode == "resume" and a.job is None:
        cli.error("接续须指定已有任务目录")
    run(a.out.resolve(), a.mode, a.job)
