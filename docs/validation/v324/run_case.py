"""监督银球几何、长度对照或工作台验收，保存输入指纹和原生回执。"""
import argparse
import io
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe, sha


def quote(p):
    return "'"+p.as_posix().replace("'","''")+"'"


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    parser.add_argument("action",choices=("main","length","ui"))
    args=parser.parse_args()
    assert args.name.isidentifier()
    out=ROOT/"outputs/v324_validation"/args.name
    out.mkdir(parents=True,exist_ok=False)
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});"
    paths=[ROOT/"buildShortSupportedSphere.m",ROOT/"exportShortSupportedSphere.m",HERE/"PLAN.md",Path(__file__).resolve(),
           *(ROOT/"cases/short_sphere"/f for f in ("test65.json","surface65.zip","surface65.json"))]
    if args.action=="ui":
        front=out/"frontend"
        keep=["web","createAerothermalWebTab.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
        with ZipFile(io.BytesIO(subprocess.check_output(["git","archive","--format=zip","v323",*keep],cwd=ROOT))) as z:
            z.extractall(front)
        paths.extend([*front.rglob("*"),ROOT/"rocketAeroHeatWorkbench.m",ROOT/"createShortSupportReview.m",HERE/"check_support_ui.m"])
        cmd+=f"addpath({quote(front)},'-begin');check_support_ui({quote(out)})"
    else:
        length=1 if args.action=="main" else .5
        cmd+=f"g=buildShortSupportedSphere({length});exportShortSupportedSphere(g,{quote(out/'supported.zip')});"
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(f):sha(f) for f in paths if f.is_file()},180 if args.action=="ui" else 120,100*1024**2)
    assert rec["status"]=="terminal" and rec["debuggerExitCode"]==0 and all(rec["inputsStable"].values())
    assert not (rec["stopReason"] or rec["supervisorError"] or rec["survivingProcessIDs"])
    assert (out/("ui.json" if args.action=="ui" else "supported.zip")).is_file()
