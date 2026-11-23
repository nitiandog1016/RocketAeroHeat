"""在同一发布前端检验侧栏原生冷算，父任务预算包含真实求解时间。"""
from pathlib import Path
import shutil
import sys

from run_ui import HERE, ROOT, quoted, observe, sha


if __name__=="__main__":
    assert len(sys.argv)==3 and all(name.isidentifier() for name in sys.argv[1:])
    base=ROOT/"outputs/v329_validation";out=base/sys.argv[1];out.mkdir(exist_ok=False)
    front=base/sys.argv[2]/"frontend"
    prior=__import__("json").loads((front.parent/"debugger_process.json").read_bytes())
    assert prior["status"]=="terminal" and prior["debuggerExitCode"]==0 and all(prior["inputsStable"].values())
    frontFiles={str(p):sha(p) for p in front.rglob("*") if p.is_file()}
    assert frontFiles and all(prior["inputSHA256"][p]==digest for p,digest in frontFiles.items())
    files=[*ROOT.glob("*ShortSupport*.m"),ROOT/"buildAerothermalViewerData.m",ROOT/"rocketAeroHeatWorkbench.m",
           ROOT/"scripts/shortSupportShared65.py",HERE/"check_cold_ui.m",HERE/"GUI_PLAN.md",*HERE.glob("*.py")]
    cmd=(f"cd({quoted(out)});addpath({quoted(ROOT)},{quoted(HERE)});addpath({quoted(front)},'-begin');check_cold_ui({quoted(out)})")
    r=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
              {**frontFiles,**{str(p):sha(p) for p in files}},660,350*1024**2)
    assert r["status"]=="terminal" and r["debuggerExitCode"]==0 and all(r["inputsStable"].values())
    assert not(r["stopReason"] or r["supervisorError"] or r["survivingProcessIDs"]) and (out/"ui.json").is_file()
