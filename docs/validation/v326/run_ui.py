"""从 v325 已发布前端启动 MATLAB 侧栏，只核对发生变化的试件页面。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/"scripts"))
from runSU2NativeDebugger import observe, sha


def matlab_path(path):
    return "'" + path.as_posix().replace("'", "''") + "'"


if __name__ == "__main__":
    assert len(sys.argv) == 2 and sys.argv[1].isidentifier()
    out = ROOT/"outputs/v326_validation"/sys.argv[1]
    out.mkdir(exist_ok=False)
    front = out/"frontend"
    keep = ["web", "createAerothermalWebTab.m", "createAerothermalGeometryView.m", "exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git", "archive", "--format=zip", "v325", *keep], cwd=ROOT))) as z:
        z.extractall(front)
    cmd = (f"cd({matlab_path(out)});addpath({matlab_path(ROOT)},{matlab_path(HERE)});"
           f"addpath({matlab_path(front)},'-begin');check_flux_ui({matlab_path(out)})")
    paths = [*front.rglob("*"), *ROOT.glob("*ShortSupport*.m"), ROOT/"rocketAeroHeatWorkbench.m",
             HERE/"check_flux_ui.m", Path(__file__), *(ROOT/"cases/short_sphere").glob("*")]
    result = observe([shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", cmd],
                     out, {str(p): sha(p) for p in paths if p.is_file()}, 180, 100*1024**2)
    assert result["status"] == "terminal" and result["debuggerExitCode"] == 0 and all(result["inputsStable"].values())
    assert not (result["stopReason"] or result["supervisorError"] or result["survivingProcessIDs"])
    assert (out/"ui.json").is_file()
