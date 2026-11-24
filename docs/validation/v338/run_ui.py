"""用冻结的三维前端和本次曲线组件，在真实桌面窗口重开带重新编号的银球。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];AREA=ROOT/"outputs/v338_validation"
sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import invoke


def run(name,archived=False):
    assert name.isidentifier()
    out=AREA/name;out.mkdir();app=out/"app";app.mkdir()
    base="66ba1dcabf0d8081f6a7e403be0666514f24d372"
    raw=subprocess.check_output(["git","archive","--format=zip",base,"--",":(top,glob)*.m","web","cases/short_sphere","cases/rocket.stl"],cwd=ROOT)
    with ZipFile(io.BytesIO(raw)) as z:z.extractall(app)
    for f in ("createAerothermalWebTab.m","updateAerothermalWebTab.m","loadShortSupportSurface.m","exportAerothermalWeb.m",
              "runShortSupportImportedCase.m","rocketAeroHeatWorkbench.m","readShortSupportAzimuth.m",
              "syncShortAzimuthTab.m","shortAzimuthPayload.m","web/short65-azimuth.html","cases/short_sphere/stationAzimuth65.zip"):
        shutil.copyfile(ROOT/f,app/f)
    with ZipFile(app/"cases/short_sphere/supported65.zip") as z:b=z.read("supported_m.stl")
    n=int.from_bytes(b[80:84],"little");assert n==32480 and len(b)==84+50*n
    (out/"input.stl").write_bytes(b if archived else b[:84]+b"".join(b[84+i*50:134+i*50] for i in reversed(range(n))))
    matlab=Path(shutil.which("matlab"))
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"restoredefaultpath;cd({quote(app)});addpath({quote(HERE)});check_ui({quote(out)},{quote(AREA/'reread01')},{int(archived)})"
    invoke("window",[matlab,"-wait","-singleCompThread","-batch",cmd],app,out,300,4,bin_directory=matlab.parent)
    assert (out/"ui.json").is_file()


if __name__=="__main__":run(sys.argv[1],"--archive" in sys.argv[2:])
