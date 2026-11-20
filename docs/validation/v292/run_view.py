"""在实际 MATLAB WebGL 控件中验收热敏珠同刻对照。前端固定至发布版，保留并发修改。"""
import argparse
import io
import shutil
import subprocess
import zipfile
from pathlib import Path

from run_step import HERE,ROOT,RUN,observe,read,sha


def verify(out,table):
    out.relative_to(RUN);out.mkdir(parents=True,exist_ok=False)
    release=out/"release";live=out/"live"
    resources=["web","createAerothermalWebTab.m","createAerothermalEvidenceView.m","exportAerothermalWeb.m"]
    assets=subprocess.check_output(["git","archive","--format=zip","v291",*resources],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(assets)) as z:z.extractall(release)
    shutil.copytree(release/"web",live/"web")
    for file in resources[1:3]:shutil.copyfile(release/file,live/file)
    observer=HERE.parent/"v290/observe_display.js"
    shutil.copyfile(observer,live/"web/observe_display.js")
    index=live/"web/index.html";text=index.read_text("utf8");assert text.count("</body>")==1
    index.write_text(text.replace("</body>",'<script src="observe_display.js"></script></body>'),encoding="utf8")
    def arg(p):return "'"+p.as_posix().replace("'","''")+"'"
    call=(f"cd({arg(out)});addpath({arg(ROOT)},{arg(HERE)});addpath({arg(release)},'-begin');"
          f"addpath({arg(live)},'-begin');check_view({arg(out)},{arg(table)})")
    files=[HERE/"PLAN.md",HERE/"check_view.m",Path(__file__).resolve(),observer,table/"wall.csv",
           *release.rglob("*"),*live.rglob("*")]
    files.extend(ROOT/(n+".m") for n in ("openThermistorArchive","attachThermistorStepComparison","readThermistorStartup",
                                         "buildAerothermalViewerData","runThermistorBeadCase","rocketAeroHeatWorkbench"))
    files.extend(ROOT/"cases/thermistor_b"/n for n in ("startup_6ms.zip","startup_6ms_halfsteps.zip"))
    rec=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",call],out,
                {str(p):sha(p) for p in files if p.is_file()},180,80*1024**2)
    assert rec["debuggerExitCode"]==0 and not rec["stopReason"],out/"debugger.log"
    assert read(out/"comparison_ui.json")["surfaceRows"]==5216
    assert sum(p.stat().st_size for p in out.rglob("*") if p.is_file())<80*1024**2


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument("output",type=Path);cli.add_argument("table",type=Path)
    a=cli.parse_args();verify(a.output.resolve(),a.table.resolve())
