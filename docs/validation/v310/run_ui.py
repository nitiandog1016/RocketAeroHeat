"""用 已发布前端和实际侧栏核对 B 型热敏珠式（28）的显示和导出。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe, sha


def run(level,name):
    out=ROOT/"outputs/v310_validation"/name
    out.mkdir(parents=True,exist_ok=False)
    front=out/"frontend"
    snapshot=subprocess.check_output(["git","archive","--format=zip","v309","web","createAerothermalWebTab.m",
                                      "createAerothermalGeometryView.m","exportAerothermalWeb.m"],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(snapshot)) as z:
        z.extractall(front)
    exporter=ROOT/"createAerothermalEvidenceView.m"
    shutil.copyfile(exporter,front/exporter.name)
    def quote(p):
        return "'"+p.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_mean_ui({quote(out)},{level})"
    files=[*front.rglob("*"),exporter,HERE/"PLAN.md",HERE/"check_mean_ui.m",Path(__file__).resolve(),
           ROOT/"cases/thermistor_b/beadExperimentGap.json",ROOT/"runThermistorLeadWorkbench.m",ROOT/"buildAerothermalViewerData.m",
           ROOT/"rocketAeroHeatWorkbench.m",ROOT/"thermistorBeadMeanConvection.m",ROOT/"phaseAwareAirState.m",
           ROOT/"data/thermo/phase_aware_dry_air_v1.json",ROOT/"aeroHeatWorkflowState.m",ROOT/"cases/thermistor_b/layerBands.json"]
    for key in ("leads3","leadsMidbody3","leadsThird3","midbody3","outer30","paired3"):
        files.extend((ROOT/"cases/thermistor_b").glob(key+".*"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(p):sha(p) for p in files if p.is_file()},180,180*1024**2)
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
    assert all(rec["inputsStable"].values()) and (out/"ui.json").is_file()
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<180*1024**2


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("level",type=int,choices=(1,2));parser.add_argument("name")
    args=parser.parse_args()
    if not args.name.isidentifier():parser.error("请使用新的目录名")
    run(args.level,args.name)
