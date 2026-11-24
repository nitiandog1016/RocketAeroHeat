"""从冻结的版本导出本试件需要的项目文件；独立安装开始时不含工具链和计算输出。"""
import io
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile
from shutil import copyfile

ROOT=Path(__file__).resolve().parents[2].parent
HERE=ROOT/"docs/validation/v336"
AREA=ROOT/"outputs/v336_validation"
BASE="59c17358d992fbba80b4d5a1a61f12e8d857da96"


def install(dest):
    dest.relative_to(AREA);dest.mkdir(parents=True,exist_ok=False)
    names=subprocess.check_output(["git","ls-tree","-r","--name-only",BASE],cwd=ROOT).decode("utf8").splitlines()
    chosen=[]
    for n in names:
        p=Path(n)
        if (len(p.parts)==1 and p.suffix==".m") or n.startswith(("web/","cases/short_sphere/","third_party/gmsh/")) or n=="cases/rocket.stl":
            chosen.append(n)
        elif n.startswith("scripts/") and p.parent==Path("scripts") and p.suffix==".py":chosen.append(n)
        elif n.startswith("docs/validation/") and len(p.parts)==4 and p.parts[2].startswith("v") and p.parts[2][1:].isdigit():
            if int(p.parts[2][1:])>=270 and (p.suffix in (".py",".hpp") or p.name=="PLAN.md"):chosen.append(n)
    for start in range(0,len(chosen),80):
        archive=subprocess.check_output(["git","archive","--format=zip",BASE,"--",*chosen[start:start+80]],cwd=ROOT)
        with ZipFile(io.BytesIO(archive)) as z:z.extractall(dest)
    # 仅叠加本轮修改的文件；其他会话的共享 Web 和研究文件不进入安装。
    changed=subprocess.check_output(["git","diff","--name-only",BASE,"--"],cwd=ROOT).decode().splitlines()
    edited={"aeroHeatWorkflowState.m","rocketAeroHeatWorkbench.m","runShortSupportImportedCase.m","shortSupportJobFrame.m",
            "docs/validation/v272/observe_flux.py","docs/validation/v329/build_shared.py","docs/validation/v329/run_shared.py","scripts/shortSupportTime65.py",
            "scripts/installShort65Kernel.py","scripts/shortSupportFirstRun.py","scripts/buildShort65Programs.py",
            "docs/validation/v336/PLAN.md","docs/validation/v336/check_clean_ui.m"}
    for n in sorted(edited):
        if (ROOT/n).is_file() and (n in changed or n not in chosen):
            p=dest/n;p.parent.mkdir(parents=True,exist_ok=True);copyfile(ROOT/n,p)
    assert not (dest/"outputs").exists() and not (dest/"tools/runtime").exists()
    (dest/"outputs").mkdir()
    (dest.parent/(dest.name+"_installed.json")).write_text(json.dumps(dict(base=BASE,trackedFiles=len(chosen),initialOutputs=[],
        runtimeInstalled=False,installation=str(dest)),indent=2),encoding="utf8")
    sources={n:hashlib.sha256((dest/n).read_bytes()).hexdigest() for n in sorted(edited)}
    (dest.parent/(dest.name+"_source.json")).write_text(json.dumps(sources,indent=2),encoding="utf8")
    print(dest,flush=True)


if __name__=="__main__":
    name=sys.argv[1];assert name.isidentifier()
    install(AREA/name)
