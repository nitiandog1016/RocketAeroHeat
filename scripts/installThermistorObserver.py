"""在没有编译工具的 Windows 安装中恢复指纹已登记的热敏珠专用观察器。"""
import argparse
from pathlib import Path
import shutil
import tempfile

from thermistorStartupCache import BUILD, CACHE, CASE, OBSERVER_SPEC, installation_lock
from thermistorSeedRecovery import identical
from thermistorStartupObserver import check
from thermistorStartupSurface import write


def install(executable):
    exe=identical(executable.resolve(),OBSERVER_SPEC["observerSHA256"])
    source=identical(CASE/"CIntegration.cpp",OBSERVER_SPEC["sourceSHA256"])
    with installation_lock("observer"):
        if not BUILD.exists():
            work=Path(tempfile.mkdtemp(prefix="publish-observer-",dir=CACHE))
            shutil.copyfile(exe,work/"observe_bead.exe")
            shutil.copyfile(source,work/"CIntegration.cpp")
            write(work/"build.json",dict(observerSHA256=OBSERVER_SPEC["observerSHA256"],
                  sourceSHA256=OBSERVER_SPEC["sourceSHA256"],origin="verified-installed-binary"))
            check(work)
            work.rename(BUILD)
        check(BUILD)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("executable",type=Path)
    install(cli.parse_args().executable)
