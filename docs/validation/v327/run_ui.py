"""用 v326 的已发布前端检查实际 MATLAB 侧栏，保留其他前端工作。"""
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


def quoted(path):
    return "'" + path.as_posix().replace("'", "''") + "'"


if __name__ == "__main__":
    assert len(sys.argv) == 2 and sys.argv[1].isidentifier()
    out = ROOT/"outputs/v327_validation"/sys.argv[1]
    out.mkdir(parents=True, exist_ok=False)
    front = out/"frontend"
    keep = ["web", "createAerothermalWebTab.m", "createAerothermalGeometryView.m", "exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git", "archive", "--format=zip", "v326", *keep], cwd=ROOT))) as z:
        z.extractall(front)
    cmd = (f"cd({quoted(out)});addpath({quoted(ROOT)},{quoted(HERE)});"
           f"addpath({quoted(front)},'-begin');check_diffusion_ui({quoted(out)})")
    files = [*front.rglob("*"), *ROOT.glob("*ShortSupport*.m"), ROOT/"rocketAeroHeatWorkbench.m",
             *HERE.glob("*.m"), Path(__file__), *(ROOT/"cases/short_sphere").glob("*")]
    r = observe([shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", cmd], out,
                {str(p):sha(p) for p in files if p.is_file()}, 180, 100*1024**2)
    assert r["status"] == "terminal" and r["debuggerExitCode"] == 0 and all(r["inputsStable"].values())
    assert not (r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
    assert (out/"ui.json").is_file()
