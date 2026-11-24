"""以冻结网页和本轮 MATLAB 完整文件，实窗验收两个时刻、两条步长与 STL 编号。"""
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];AREA=ROOT/"outputs/v339_validation"
sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import invoke


def run(name,mode):
    assert name.isidentifier() and mode in ("archive","local","early")
    out=AREA/name;out.mkdir();app=out/"app";app.mkdir()
    base="3935850822bf743f197350e70b5a7d2c27603311"
    raw=subprocess.check_output(["git","archive","--format=zip",base,"--",":(top,glob)*.m","web","cases/short_sphere","cases/rocket.stl"],cwd=ROOT)
    with ZipFile(io.BytesIO(raw)) as z:z.extractall(app)
    for f in ("aeroHeatWorkflowState.m","attachShortStationGradients.m","buildAerothermalViewerData.m","exportWallPairEvidence.m",
              "loadShortSupportSurface.m","readShortSupportAzimuth.m","readShortSupportTimeDifference.m","readShortSupportTimePair.m",
              "rocketAeroHeatWorkbench.m","runShortSupportImportedCase.m","shortAzimuthPayload.m","shortSupportImportBoundary.m",
              "shortSupportJobFrame.m","validateShortTimeChain.m","web/short65-azimuth.html"):
        shutil.copyfile(ROOT/f,app/f)
    later=ROOT/"cases/short_sphere/later1"
    if later.is_dir():shutil.copytree(later,app/"cases/short_sphere/later1")
    with ZipFile(app/"cases/short_sphere/supported65.zip") as z:stl=z.read("supported_m.stl")
    n=int.from_bytes(stl[80:84],"little");assert n==32480
    (out/"input.stl").write_bytes(stl if mode!="local" else stl[:84]+b"".join(stl[84+i*50:134+i*50] for i in reversed(range(n))))
    matlab=Path(shutil.which("matlab"))
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"restoredefaultpath;cd({quote(app)});addpath({quote(HERE)});check_ui({quote(out)},'{mode}',{quote(AREA/'later01')})"
    invoke("window",[matlab,"-wait","-singleCompThread","-batch",cmd],app,out,300,4,bin_directory=matlab.parent)
    assert (out/"ui.json").is_file()


if __name__=="__main__":run(*sys.argv[1:])
