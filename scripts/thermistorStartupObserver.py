"""只为 172,482 节点热敏珠构建固定端点观察器，支持独立安装的 SU2 构建工具目录。"""
import gzip
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

from thermistorStartupCache import BUILD, CACHE, CASE, OBSERVER_SPEC, PREVIOUS, ROOT, installation_lock
from thermistorSeedRecovery import identical
from thermistorStartupSurface import read, write

sys.path.insert(0,str(ROOT/"docs/validation/v272"))
from observe_flux import invoke


def check(directory):
    build=read(directory/"build.json")
    for key,name in (("observerSHA256","observe_bead.exe"),("sourceSHA256","CIntegration.cpp")):
        if build[key] != OBSERVER_SPEC[key]:
            raise ValueError("观察器构建回执与本工况不符")
        identical(directory/name,build[key])


def build_observer(out, kit):
    """将原始编译命令中的安装和输出路径重新定位，并按原来源逐项验证。"""
    recipes=json.loads(gzip.decompress((PREVIOUS/"evidence/processes.json.gz").read_bytes()))
    compile_rec=recipes["build01/compile_process.json"]
    link_rec=recipes["build01/link_process.json"]
    old_build=Path(compile_rec["cwd"])
    old_output=Path(link_rec["args"][2]).parent
    build=kit/old_build.name
    compiler=kit/Path(compile_rec["args"][0]).relative_to(old_build.parent)
    identical(compiler,OBSERVER_SPEC["compilerSHA256"])
    cache_file=ROOT/"docs/validation/v272/evidence/build.json.gz"
    inputs=json.loads(gzip.decompress(cache_file.read_bytes()))["cachedBuildInputs"]
    for file,digest in inputs.items():
        identical(kit/Path(file).relative_to(old_build.parent),digest)
    shutil.copyfile(identical(CASE/"CIntegration.cpp",OBSERVER_SPEC["sourceSHA256"]),out/"CIntegration.cpp")
    def relocate(arg):
        path=Path(arg)
        if not path.is_absolute():
            return arg
        if path.parent==old_output:
            return str(out/path.name)
        return str(kit/path.relative_to(old_build.parent))
    for name,rec,seconds in (("compile",compile_rec,120),("link",link_rec,60)):
        assert Path(rec["cwd"])==old_build and rec["exitCode"]==0 and not rec["livePids"]
        args=[relocate(v) for v in rec["args"]]
        invoke(name,args,build,out,seconds,3,bin_directory=compiler.parent)
    identical(out/"observe_bead.exe",OBSERVER_SPEC["observerSHA256"])
    write(out/"build.json",dict(OBSERVER_SPEC,origin="rebuilt",buildKit=str(kit),verifiedInputs=len(inputs)))


def ensure_observer():
    with installation_lock("observer"):
        if not BUILD.exists():
            legacy=ROOT/"outputs/v293_validation/build01"
            use_legacy=(legacy/"observe_bead.exe").is_file()
            kit=Path(os.environ.get("RAH_THERMISTOR_SU2_BUILD",ROOT/"tools/runtime/su2-thermal-build")).resolve()
            if not use_legacy and not kit.is_dir():
                raise FileNotFoundError(
                    "缺少本工况的热流观察器。可用 python scripts/installThermistorObserver.py "
                    "安装已核验的 observe_bead.exe，或将 RAH_THERMISTOR_SU2_BUILD 指向已安装的构建工具目录。"
                    "详见 docs/validation/v296/INSTALL.md。"
                )
            work=Path(tempfile.mkdtemp(prefix="observer-",dir=CACHE))
            if use_legacy:
                for name,key in (("observe_bead.exe","observerSHA256"),("CIntegration.cpp","sourceSHA256")):
                    shutil.copyfile(identical(legacy/name,OBSERVER_SPEC[key]),work/name)
                write(work/"build.json",dict(observerSHA256=OBSERVER_SPEC["observerSHA256"],
                      sourceSHA256=OBSERVER_SPEC["sourceSHA256"],origin="verified-local-migration"))
            else:
                build_observer(work,kit)
            publish=Path(tempfile.mkdtemp(prefix="publish-observer-",dir=CACHE))
            for name in ("observe_bead.exe","CIntegration.cpp"):
                os.link(work/name,publish/name)
            write(publish/"build.json",dict(read(work/"build.json"),sourceDirectory=str(work)))
            publish.rename(BUILD)
        check(BUILD)


if __name__=="__main__":
    if not __debug__:
        raise RuntimeError("热敏珠构建审核需要开启 Python 断言")
    ensure_observer()
