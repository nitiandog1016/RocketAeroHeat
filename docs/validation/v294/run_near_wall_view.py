"""在实际 MATLAB 侧栏中核对热敏珠近壁对照；已发布前端与并发开发文件相互独立。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import zipfile

from advance import HERE, ROOT, RUN, observe, read, sha


def matlab_path(path):
    return "'" + path.as_posix().replace("'", "''") + "'"


def verify(out):
    out.relative_to(RUN)
    area = out.parent
    accepted = read(area / "second_energy/accepted.json")
    table = area / "second_energy"
    package = ROOT / "cases/thermistor_b/startup_6ms_nearwall56.zip"
    assert accepted["accepted"] and sha(package) == sha(area / "compare" / package.name)
    out.mkdir(exist_ok=False)
    release, live = out / "release", out / "live"
    resources = ["web", "createAerothermalWebTab.m", "createAerothermalEvidenceView.m", "exportAerothermalWeb.m"]
    archive = subprocess.check_output(["git", "archive", "--format=zip", "v293", *resources], cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        z.extractall(release)
    shutil.copytree(release / "web", live / "web")
    for name in resources[1:3]:
        shutil.copyfile(release / name, live / name)
    witness = HERE.parent / "v290/observe_display.js"
    shutil.copyfile(witness, live / "web/observe_display.js")
    index = live / "web/index.html"
    html = index.read_text("utf8")
    assert html.count("</body>") == 1
    index.write_text(html.replace("</body>", '<script src="observe_display.js"></script></body>'), encoding="utf8")
    command = (f"cd({matlab_path(out)});addpath({matlab_path(ROOT)},{matlab_path(HERE)});"
               f"addpath({matlab_path(release)},'-begin');addpath({matlab_path(live)},'-begin');"
               f"check_near_wall_view({matlab_path(out)},{matlab_path(table)})")
    files = [HERE / "PLAN.md", HERE / "check_near_wall_view.m", Path(__file__).resolve(), witness,
             table / "wall.csv", table / "matched_temperature.csv", package, *release.rglob("*"), *live.rglob("*")]
    files.extend(ROOT / (name + ".m") for name in (
        "rocketAeroHeatWorkbench", "runThermistorBeadCase", "buildAerothermalViewerData",
        "openThermistorArchive", "readThermistorStartup", "attachThermistorStepComparison"))
    files.extend(package.parent / name for name in ("startup_6ms.zip", "startup_6ms_halfsteps.zip"))
    receipt = observe([shutil.which("matlab"), "-wait", "-singleCompThread", "-batch", command], out,
                      {str(p): sha(p) for p in files if p.is_file()}, 180, 80 * 1024**2)
    assert receipt["debuggerExitCode"] == 0 and not receipt["stopReason"], out / "debugger.log"
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) < 80 * 1024**2
    assert read(out / "near_wall_ui.json")["surfaceRows"] == 5216


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    verify(parser.parse_args().output.resolve())
