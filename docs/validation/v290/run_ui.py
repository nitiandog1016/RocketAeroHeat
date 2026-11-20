"""使用真实 MATLAB 组件与已发布网页，避免把并行前端改动混入验证。"""
import argparse
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe
from runEmbeddedSU2 import sha


def run(out,reference):
    out.mkdir(parents=True,exist_ok=False)
    release=out/"release";live=out/"live"
    resources=["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","exportAerothermalWeb.m"]
    data=subprocess.check_output(["git","archive","--format=zip","v289",*resources],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(data)) as archive:archive.extractall(release)
    for name in ("index.html","export.js","evidence.js"):
        shutil.copyfile(ROOT/"web"/name,release/"web"/name)
    shutil.copytree(release/"web",live/"web")
    for key in resources[1:3]:shutil.copyfile(release/key,live/key)
    shutil.copyfile(HERE/"observe_display.js",live/"web/observe_display.js")
    index=live/"web/index.html";html=index.read_text("utf8")
    assert html.count("</body>")==1
    index.write_text(html.replace("</body>",'<script src="observe_display.js"></script></body>'),encoding="utf8")
    def arg(p):return "'"+p.as_posix().replace("'","''")+"'"
    call=(f"cd({arg(out)});addpath({arg(ROOT)},{arg(HERE)});addpath({arg(release)},'-begin');"
          f"addpath({arg(live)},'-begin');check_snapshot({arg(out)},{arg(reference)})")
    files=[HERE/"PLAN.md",HERE/"check_snapshot.m",HERE/"observe_display.js",Path(__file__).resolve(),
           ROOT/"cases/thermistor_b/startup_6ms.zip",reference/"wall.csv",*release.rglob("*"),*live.rglob("*")]
    files.extend(ROOT/(n+".m") for n in ("openThermistorArchive","readThermistorStartup","runThermistorBeadCase",
                    "buildAerothermalViewerData","rocketAeroHeatWorkbench","runHybridSurfaceCase"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",call],out,
                {str(p):sha(p) for p in files if p.is_file()},180,80*1024**2)
    assert not rec["stopReason"] and rec["debuggerExitCode"]==0,out/"debugger.log"
    proof=json.loads((out/"ui.json").read_bytes())
    assert proof["faces"]==5216 and proof["sectionRows"]>0 and proof["webglError"]==0
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<80*1024**2


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("output",type=Path);cli.add_argument("reference",type=Path)
    args=cli.parse_args();out=args.output.resolve();out.relative_to(ROOT/"outputs/v290_validation")
    run(out,args.reference.resolve())
