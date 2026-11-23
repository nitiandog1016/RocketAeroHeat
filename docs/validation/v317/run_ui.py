"""使用已发布的 WebGL 前端与本轮 MATLAB 页面，实际操作长引线热边界。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe,sha


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("name");p.add_argument("job",type=Path);a=p.parse_args()
    if not a.name.isidentifier():p.error("使用新的输出目录名")
    out=ROOT/"outputs/v317_validation"/a.name;out.mkdir(parents=True,exist_ok=False);front=out/"frontend"
    fixed=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v316",*fixed],cwd=ROOT))) as z:z.extractall(front)
    def quote(path):return "'"+path.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_wire_ui({quote(out)},{quote(a.job.resolve())})"
    paths=[*front.rglob("*"),Path(__file__).resolve(),HERE/"check_wire_ui.m",HERE/"PLAN.md"]
    paths.extend(ROOT/n for n in ("rocketAeroHeatWorkbench.m","runThermistorWireWorkbench.m","thermistorWorkbenchBoundary.m","runAeroHeatSingleCase.m",
                                  "exportWallPairEvidence.m","createAerothermalEvidenceView.m","buildAerothermalViewerData.m","aeroHeatWorkflowState.m",
                                  "thermistorWireBands.m","scripts/thermistorWireWorkbench.py","scripts/runThermistorWireCase.py",
                                  "scripts/runThermistorMeasuredCold.py","cases/thermistor_b/longWire12.json","cases/thermistor_b/longWire12.zip"))
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,{str(k):sha(k) for k in paths if k.is_file()},240,200*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"])
    assert (out/"ui.json").is_file() and sum(f.stat().st_size for f in out.rglob("*") if f.is_file())<200*1024**2
