"""在冻结前端和本次 MATLAB 文件的独立目录中验证实窗及离线导出的七站对照。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
AREA=ROOT/"outputs/v337_validation"
sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import invoke


def run(name):
    assert name.isidentifier()
    out=AREA/name;out.mkdir();app=out/"app";app.mkdir()
    base="139bbe3fd7fd840b53634ea6d3b962a82b543612"
    raw=subprocess.check_output(["git","archive","--format=zip",base,"--",":(top,glob)*.m","web","cases/short_sphere","cases/rocket.stl"],cwd=ROOT)
    with ZipFile(io.BytesIO(raw)) as z:z.extractall(app)
    for f in ("rocketAeroHeatWorkbench.m","runShortSupportImportedCase.m","attachShortStationGradients.m",
              "addShortSupportSurfaceTab.m","exportAerothermalWeb.m","loadShortSupportSurface.m","cases/short_sphere/stationGradients65.json"):
        shutil.copyfile(ROOT/f,app/f)
    with ZipFile(app/"cases/short_sphere/supported65.zip") as z:z.extract("supported_m.stl",out)
    matlab=Path(shutil.which("matlab"))
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"restoredefaultpath;cd({quote(app)});addpath({quote(HERE)});check_ui({quote(out)},{quote(AREA/'reread01')})"
    invoke("window",[matlab,"-wait","-singleCompThread","-batch",cmd],app,out,300,4,bin_directory=matlab.parent)
    assert (out/"ui.json").is_file() and (out/"workbench/index.html").is_file()


if __name__=="__main__":run(sys.argv[1])
