"""隔离并发前端改动，使用发布的 WebGL 页面验收热敏珠风洞测点两种计算入口。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe,sha


def run(name):
    out=ROOT/"outputs/v315_validation"/name
    out.mkdir(parents=True,exist_ok=False)
    front=out/"frontend"
    files=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v314",*files],cwd=ROOT))) as z:
        z.extractall(front)
    exporter=ROOT/"createAerothermalEvidenceView.m"
    shutil.copyfile(exporter,front/exporter.name)
    def quote(p):
        return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_point_ui({quote(out)})"
    paths=[*front.rglob("*"),exporter,HERE/"PLAN.md",HERE/"check_point_ui.m",Path(__file__).resolve()]
    for stem in ("windPoint12","sidewindAngular96","leadsThird3"):
        paths.extend((ROOT/"cases/thermistor_b").glob(stem+".*"))
    paths.extend(ROOT/(k+".m") for k in ("rocketAeroHeatWorkbench","runThermistorSidewindWorkbench","runThermistorMeasuredPoint","thermistorBeadMeanConvection","thermistorWorkbenchBoundary",
                 "aeroHeatWorkflowState","exportWallPairEvidence","buildAerothermalViewerData","openThermistorArchive","prepareThermistorLeadMesh"))
    result=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                   {str(p):sha(p) for p in paths if p.is_file()},240,200*1024**2)
    assert result["status"]=="terminal" and result["debuggerExitCode"]==0 and all(result["inputsStable"].values())
    assert not (result["stopReason"] or result["supervisorError"] or result["survivingProcessIDs"])
    assert (out/"ui.json").is_file() and sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<200*1024**2


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("name");a=p.parse_args()
    if not a.name.isidentifier():p.error("使用新的输出目录名")
    run(a.name)
