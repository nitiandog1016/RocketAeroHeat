"""从指定 Git 版本运行 RFZ 边界工作台验收，只复用已冻结的米制 STL。"""
import argparse
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "docs/validation/v272"))
from observe_flux import invoke


def run(out, revision):
    out = out.resolve()
    out.relative_to(ROOT / "outputs")
    out.mkdir(parents=True, exist_ok=False)
    app = out / "app"
    app.mkdir()
    tree = subprocess.check_output(["git","rev-parse","--verify","--end-of-options",revision+"^{tree}"],cwd=ROOT,text=True).strip()
    paths = [":(top,glob)*.m", "web", "cases/rocket.stl", "cases/rfz_rec",
             "docs/validation/v343/check_boundary_ui.m", "docs/validation/v342/check_workbench.m"]
    data = subprocess.check_output(["git","archive","--format=zip",tree,"--",*paths],cwd=ROOT)
    with ZipFile(io.BytesIO(data)) as archive:
        archive.extractall(app)
    relative = Path("outputs/rfz_rec_geometry/RFZ-REC_m.stl")
    (app/relative).parent.mkdir(parents=True)
    os.link(ROOT/relative, app/relative)
    prior = out / "previous"
    prior.mkdir()
    def quote(path):
        return "'" + path.as_posix().replace("'","''") + "'"
    expr = (f"restoredefaultpath;cd({quote(app)});addpath('docs/validation/v343','docs/validation/v342');"
            f"check_boundary_ui({quote(out)});check_workbench({quote(prior)})")
    matlab = Path(shutil.which("matlab"))
    invoke("window",[matlab,"-wait","-singleCompThread","-batch",expr],app,out,240,4,bin_directory=matlab.parent)
    assert (out/"ui.json").is_file() and (prior/"workbench.json").is_file()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--revision",default="HEAD")
    args = parser.parse_args()
    run(args.output,args.revision)
