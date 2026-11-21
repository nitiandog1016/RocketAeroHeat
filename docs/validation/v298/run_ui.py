"""在已发布的前端资源上运行 MATLAB 五步侧栏与真实截面导出测试。"""
import io
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

from check_mesh import HERE,ROOT
from runEmbeddedSU2 import sha
from runSU2NativeDebugger import observe


def run(out,kind):
    out.relative_to(ROOT/"outputs/v298_validation");out.mkdir(exist_ok=False)
    front=out/"frontend"
    resources=["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","createAerothermalGeometryView.m","exportAerothermalWeb.m"]
    data=subprocess.check_output(["git","archive","--format=zip","v297",*resources],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(data)) as z:z.extractall(front)
    def quote(path):return "'"+path.as_posix().replace("'","''")+"'"
    cmd=f"cd({quote(out)});addpath({quote(ROOT)},{quote(HERE)});addpath({quote(front)},'-begin');check_pair_ui({quote(out)},'{kind}')"
    sources=[HERE/"check_pair_ui.m",HERE/"PLAN.md",Path(__file__).resolve(),*front.rglob("*")]
    sources.extend(ROOT/p for p in ("rocketAeroHeatWorkbench.m","attachThermistorPairedControl.m","buildAerothermalViewerData.m","openThermistorArchive.m",
                   "runThermistorLeadWorkbench.m","runThermistorBeadCase.m","cases/thermistor_b/paired3.json","cases/thermistor_b/paired3.zip"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",cmd],out,
                {str(p):sha(p) for p in sources if p.is_file()},180,200*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"] and (out/"ui.json").is_file(),out/"debugger.log"


if __name__=="__main__":
    assert sys.argv[2] in ("wired","historical")
    run(Path(sys.argv[1]).resolve(),sys.argv[2])
