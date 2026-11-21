"""在单独的项目根目录验收热敏珠输入恢复和五步工作台，保留原工作区的历史输出。"""
import argparse
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
AREA = ROOT/"outputs/v296_validation"
DEST = AREA/"desktop"
sys.path.insert(0,str(ROOT/"scripts"))
from runSU2NativeDebugger import observe
from runEmbeddedSU2 import sha

OVERLAY = ["scripts/installSU2Runtime.py","scripts/installThermistorObserver.py",
           "scripts/thermistorSeedRecovery.py","scripts/thermistorStartupCache.py","scripts/thermistorStartupObserver.py",
           "scripts/runThermistorStartupJob.py","docs/validation/v272/observe_flux.py","docs/validation/v293/check_state.py",
           "docs/validation/v293/read_wall.py","docs/validation/v294/advance.py","docs/validation/v294/verify.py",
           "prepareThermistorWorkbenchMesh.m","verifyThermistorWorkbenchMesh.m","runThermistorStartupWorkbench.m",
           "rocketAeroHeatWorkbench.m"]
PROOFS = {"v272/evidence/build.json.gz","v293/PLAN.md","v293/evidence/processes.json.gz","v293/evidence/grid.json",
          "v293/evidence/readout_build.json","v293/evidence/first_state.json","v293/evidence/first_history.csv.gz",
          "v294/PLAN.md","v294/evidence/local_assets.json","v295/check_workbench.m"}


def setup():
    DEST.mkdir()
    scopes=["*.m","scripts/*.py","web","cases/thermistor_b","cases/rocket.stl"]
    scopes.extend("docs/validation/"+v for v in ("v272","v287","v290","v291","v292","v293","v294","v295"))
    archive=subprocess.check_output(["git","archive","--format=zip","v295",*scopes],cwd=ROOT)
    with zipfile.ZipFile(io.BytesIO(archive)) as z:
        for item in z.infolist():
            name=item.filename
            if item.is_dir():continue
            if name.startswith("docs/validation/"):
                if not name.endswith(".py") and name[len("docs/validation/"):] not in PROOFS:continue
            elif name.endswith(".m") and "/" in name:continue
            z.extract(item,DEST)
    for name in OVERLAY:
        shutil.copyfile(ROOT/name,DEST/name)
    shutil.copytree(ROOT/"cases/thermistor_b/startup6",DEST/"cases/thermistor_b/startup6")
    (DEST/"outputs").mkdir()
    # 用同一官方包安装运行库；绝对可执行路径属于本次安装，清单独立写入。
    rel=Path("tools/runtime/su2/8.6.0")
    runtime=json.loads((ROOT/rel/"runtime.json").read_bytes())
    for record in runtime["files"]:
        source=ROOT/rel/record["path"];target=DEST/rel/record["path"]
        assert sha(source)==record["sha256"]
        target.parent.mkdir(parents=True,exist_ok=True)
        os.link(source,target)
    runtime["executable"]=str(DEST/rel/"bin/SU2_CFD.exe")
    (DEST/rel/"runtime.json").write_text(json.dumps(runtime,indent=2),encoding="utf8")


def quote(path):
    return "'"+path.as_posix().replace("'","''")+"'"


def run_ui(out):
    out.mkdir()
    # 保留真正启动这次测试的安装脚本，后续可独立回查运行回执中的源码指纹。
    shutil.copyfile(Path(__file__).resolve(),out/"executed_install_check.py")
    # 已安装官方内核和固定编译环境，小型工作台输入来自本安装，不提供旧计算目录。
    os.environ["RAH_THERMISTOR_SU2_BUILD"]=str(ROOT/"tools/runtime/su2-thermal-build")
    command=(f"cd({quote(DEST)});addpath({quote(DEST/'docs/validation/v295')});"
             f"assert(strcmp(which('rocketAeroHeatWorkbench'),fullfile(pwd,'rocketAeroHeatWorkbench.m')));"
             f"check_workbench({quote(out)},'run',{quote(DEST/'outputs/v294_validation')})")
    files=[DEST/name for name in OVERLAY]
    files.extend([HERE/"PLAN.md",Path(__file__).resolve(),DEST/"docs/validation/v295/check_workbench.m",
                  *list((DEST/"cases/thermistor_b/startup6").glob("*"))])
    seconds=180+2*(2*1100+100+2*180)+2*180+120+60
    receipt=observe([shutil.which("matlab"),"-wait","-singleCompThread","-batch",command],out,
                    {str(p):sha(p) for p in files if p.is_file()},seconds,80*1024**2)
    assert receipt["debuggerExitCode"]==0 and not receipt["stopReason"] and (out/"workflow.json").is_file()


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("mode",choices=("setup","run"))
    cli.add_argument("--out",type=Path)
    cli.add_argument("--install",type=Path)
    a=cli.parse_args()
    if a.install is not None:
        DEST=a.install.resolve()
        assert DEST.parent==AREA
    if a.mode=="setup":
        setup()
    else:
        out=a.out.resolve()
        out.relative_to(AREA)
        run_ui(out)
