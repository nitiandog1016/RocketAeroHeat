"""限制一个低速风洞工作台回归的时间、线程、内存和输出，记录实际 MATLAB 终态。"""
import argparse
from pathlib import Path
import shutil
import sys

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
import runSU2RocketPhaseV222 as runner


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project",type=Path)
    parser.add_argument("name")
    parser.add_argument("action",choices=("solve","archive"))
    parser.add_argument("--without-comparison",action="store_true")
    parser.add_argument("--momentum-study",action="store_true",help="复核 v357 真实动量网格对照场及其来源")
    args=parser.parse_args()
    project=args.project.resolve();phase="v357" if args.momentum_study else "v356"
    runner.BASE=ROOT/"outputs"/(phase+"_validation")
    dest=runner.BASE/args.name;dest.mkdir(parents=True,exist_ok=False);runner.CASE=dest
    files=["main.m","buildAerothermalViewerData.m","rocketAeroHeatWorkbench.m","runSimpleLauncherWindTunnel.m",
           "scripts/simpleLauncherWakePrimary.py","scripts/runSimpleLauncherWakeCase.py",
           "docs/validation/v354/finalize.py","docs/validation/v354/heat_step.py","docs/validation/v355/refined.py",
           "docs/validation/v355/compare.py","docs/validation/v354/flux_transfer.py","scripts/simpleLauncherAngularHeat.py",
           "cases/simple_launcher/wake/condition.json","cases/simple_launcher/wake/flow_index.json"]
    sources=[DOC/"PLAN.md",DOC/"check_primary.m"]+[project/file for file in files]
    if args.momentum_study:sources+=[DOC.parent/"v357/WORKBENCH.md",project/"cases/simple_launcher/wake/momentum_grid.zip"]
    for file in ("receipt.json","fine_primary.json","wind_tunnel.zip","wake_heat_grid.zip"):
        path=project/"cases/simple_launcher/wake"/file
        if path.is_file():
            sources.append(path)
    cmd=[shutil.which("matlab"),"-batch",f"maxNumCompThreads(1); addpath('{DOC.as_posix()}'); "+
         f"check_primary('{dest.as_posix()}','{args.action}',{str(not args.without_comparison).lower()}"+
         (",true)" if args.momentum_study else ")")]
    limit=240 if args.momentum_study and args.action=="archive" else 300
    sys.exit(not runner.supervise(args.name,cmd,limit,300,sources,project,"parent"))
