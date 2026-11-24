"""从普通米制 STL 运行工作台实窗，分别检验发布归档与反向面序的本地计算。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];AREA=ROOT/"outputs/v340_validation"
sys.path.insert(0,str(HERE.parent/"v272"))
from observe_flux import invoke


def run(name,mode,job):
    assert name.isidentifier() and mode in ("archive","local")
    out=AREA/name;out.mkdir();app=out/"app";app.mkdir()
    files=[":(top,glob)*.m","web","cases/short_sphere","cases/rocket.stl"]
    frozen=subprocess.check_output(["git","archive","--format=zip","eeff6cce24c5d8e07428fbadd6190f8be728ea51","--",*files],cwd=ROOT)
    with ZipFile(io.BytesIO(frozen)) as z:z.extractall(app)
    for path in ("aeroHeatWorkflowState.m","buildAerothermalViewerData.m","exportAerothermalWeb.m","exportWallPairEvidence.m",
                 "loadShortSupportSurface.m","readShortSupportPowerCase.m","rocketAeroHeatWorkbench.m","runShortSupportImportedCase.m",
                 "shortSupportImportBoundary.m","shortSupportJobFrame.m","cases/short_sphere/calorimetry65.json","cases/short_sphere/calorimetry65.zip"):
        shutil.copyfile(ROOT/path,app/path)
    with ZipFile(app/"cases/short_sphere/supported65.zip") as z:stl=z.read("supported_m.stl")
    n=int.from_bytes(stl[80:84],"little");assert n==32480
    if mode=="local":stl=stl[:84]+b"".join(stl[84+50*i:134+50*i] for i in reversed(range(n)))
    (out/"input.stl").write_bytes(stl)
    matlab=Path(shutil.which("matlab"))
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"restoredefaultpath;cd({quote(app)});addpath({quote(HERE)});check_power_ui({quote(out)},'{mode}',{quote(job.resolve())})"
    invoke("window",[matlab,"-wait","-singleCompThread","-batch",cmd],app,out,300,4,bin_directory=matlab.parent)
    assert (out/"ui.json").is_file()


if __name__=="__main__":run(sys.argv[1],sys.argv[2],Path(sys.argv[3]))
