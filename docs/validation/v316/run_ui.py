"""用固定发布前端实际操作 MATLAB 侧栏，验证冷态暂停或原面结果接续。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from runSU2NativeDebugger import observe, sha


def run(name, mode, job):
    out = ROOT / "outputs/v316_validation" / name
    out.mkdir(parents=True, exist_ok=False)
    front = out / "frontend"
    files = ["web", "createAerothermalWebTab.m", "createAerothermalGeometryView.m", "createAerothermalEvidenceView.m", "exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git", "archive", "--format=zip", "v315", *files], cwd=ROOT))) as bundle:
        bundle.extractall(front)
    def quote(path):
        return "'" + path.as_posix().replace("'", "''") + "'"
    cmd = f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_cold_ui({quote(out)},'{mode}',{quote(job)})"
    paths = [*front.rglob("*"), Path(__file__).resolve(), HERE / "check_cold_ui.m", HERE / "PLAN.md"]
    paths.extend(ROOT / name for name in ("rocketAeroHeatWorkbench.m","runThermistorMeasuredPoint.m","aeroHeatWorkflowState.m",
                 "runAeroHeatSingleCase.m","exportWallPairEvidence.m","buildAerothermalViewerData.m","scripts/thermistorMeasuredCold.py",
                 "scripts/thermistorMeasuredCollect.py","scripts/thermistorCrossflowSurface.py","scripts/runThermistorMeasuredCold.py"))
    result = observe([shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", cmd], out,
                     {str(p): sha(p) for p in paths if p.is_file()}, 240, 200*1024**2)
    assert result["status"] == "terminal" and result["debuggerExitCode"] == 0 and all(result["inputsStable"].values())
    assert not (result["stopReason"] or result["supervisorError"] or result["survivingProcessIDs"])
    assert (out / "ui.json").is_file() and sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) < 200*1024**2


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("name");p.add_argument("mode", choices=("controls", "result"));p.add_argument("job", type=Path)
    args = p.parse_args()
    if not args.name.isidentifier():p.error("使用新的输出目录名")
    run(args.name, args.mode, args.job.resolve())
