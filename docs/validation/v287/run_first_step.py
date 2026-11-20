"""分别监督热敏珠全周建网和官方求解器的首个物理步，保留两个进程的真实资源回执。"""
import argparse
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe
from thermistorBead import ENGINE, ENGINE_SHA, sha
from thermistorFullRing import PARENT, SOURCE, config


def run(stage, name, prepared=None):
    out=ROOT/"outputs/v287_validation"/name
    out.mkdir(parents=True,exist_ok=False)
    here=Path(__file__).resolve().parent
    src=ROOT/"scripts/thermistorFullRing.py"
    files=[here/"PLAN.md",Path(__file__).resolve(),src,ROOT/"scripts/cm011TunnelSector.py",
           ROOT/"scripts/cm011Tunnel.py",ROOT/"scripts/thermistorBead.py",ROOT/"scripts/runSU2NativeDebugger.py"]
    if stage=="prepare":
        files.extend(PARENT/key for key in SOURCE)
        command=[sys.executable,str(src),str(out)]
        seconds=300
    else:
        folder=ROOT/"outputs/v287_validation"/prepared
        import json
        receipt=json.loads((folder/"debugger_process.json").read_bytes())
        assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and receipt["stopReason"] is None
        assert all(receipt["inputsStable"].values()) and sha(ENGINE)==ENGINE_SHA
        for key in ("geometry.json","nodes.npz","body.json","mesh.su2","seed_00000.csv"):
            os.link(folder/key,out/key)
            files.append(out/key)
        # 非定常输出由物理步编号命名。使用独立配置保留原网与失败现场。
        geom=json.loads((out/"geometry.json").read_bytes())
        (out/"case.cfg").write_text(config(geom["reynoldsD"],geom["pecletD"]),encoding="ascii")
        files.extend([out/"case.cfg",PARENT/"case.cfg"])
        files.extend([folder/"debugger_process.json",ENGINE])
        command=[str(ENGINE),"case.cfg","--threads","1"]
        seconds=360
    for key in ("MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[key]="1"
    result=observe(command,out,{str(p):sha(p) for p in files},seconds,180*1024**2)
    print({k:result[k] for k in ("status","debuggerExitCode","stopReason","elapsedSeconds","peakFamilyRSSBytes")})
    return result["debuggerExitCode"]==0 and result["stopReason"] is None


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=["prepare","solve"])
    parser.add_argument("name")
    parser.add_argument("--prepared")
    args=parser.parse_args()
    if args.stage=="solve" and not args.prepared:
        parser.error("求解需要已核对的建网目录")
    for component in (args.name,args.prepared):
        if component is not None and (Path(component).name!=component or component in ("",".","..")):
            parser.error("请输入单层目录名")
    raise SystemExit(0 if run(args.stage,args.name,args.prepared) else 1)
