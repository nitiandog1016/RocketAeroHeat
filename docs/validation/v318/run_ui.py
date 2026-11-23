"""用 v317 的固定 WebGL 前端实测本轮计算、独立核对、测温表与导出。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe,sha


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("name");a=p.parse_args()
    if not a.name.isidentifier():p.error("使用新目录名")
    out=ROOT/"outputs/v318_validation"/a.name;out.mkdir(parents=True,exist_ok=False);front=out/"frontend"
    fixed=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v317",*fixed],cwd=ROOT))) as z:z.extractall(front)
    def quote(p):return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_ui({quote(out)})"
    paths=[*front.rglob("*"),HERE/"check_ui.m",HERE/"PLAN.md",Path(__file__).resolve()]
    paths.extend(ROOT/n for n in ("rocketAeroHeatWorkbench.m","runThermistorWireWorkbench.m","thermistorWorkbenchBoundary.m",
                 "aeroHeatWorkflowState.m","applyThermistorThermometry.m","thermistorResistanceModel.m","exportWallPairEvidence.m",
                 "buildAerothermalViewerData.m","createAerothermalEvidenceView.m","openThermistorArchive.m"))
    for stem in ("longWire12","resistance12","windPoint12"):paths.extend((ROOT/"cases/thermistor_b").glob(stem+".*"))
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(p):sha(p) for p in paths if p.is_file()},240,200*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
    assert (out/"ui.json").is_file() and sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<200*1024**2
