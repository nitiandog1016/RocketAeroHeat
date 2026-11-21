"""核对官方程序清单迁移、无编译器安装专用观察器，以及错误输入的拒绝路径。"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha


def run():
    install=ROOT/"outputs/v296_validation/desktop"
    extra=install/"tools/runtime/su2/8.6.0"
    original=ROOT/"tools/runtime/su2/8.6.0/runtime.json"
    manifest=extra/"runtime.json"
    # 此测试安装尚未运行求解，使用有旧绝对路径的真实清单，检验安装器对硬链接的隔离。
    assert not (install/"outputs/v294_validation").exists()
    manifest.rename(extra/"initial-test-manifest.json")
    os.link(original,manifest)
    shared=sha(original)
    for name in ("installSU2Runtime.py","installThermistorObserver.py","thermistorStartupObserver.py","thermistorStartupCache.py"):
        shutil.copyfile(ROOT/"scripts"/name,install/"scripts"/name)
    def invoke(name,*args,ok=True):
        proc=subprocess.run([sys.executable,str(install/"scripts"/name),*map(str,args)],
                            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,cwd=install)
        assert (proc.returncode==0)==ok,proc.stdout
        return proc.stdout
    invoke("installSU2Runtime.py")
    assert sha(original)==shared and not os.path.samefile(original,manifest)
    assert Path(json.loads(manifest.read_bytes())["executable"])==extra/"bin/SU2_CFD.exe"
    new=sha(manifest);invoke("installSU2Runtime.py");assert sha(manifest)==new
    cache=install/"outputs/thermistor_inputs/B6-56"
    bad=install/"cases/thermistor_b/startup6/base6.cfg"
    assert not (cache/"initial").exists() and not (cache/"observer").exists()
    invoke("thermistorStartupCache.py","--install-seed",bad,ok=False)
    invoke("installThermistorObserver.py",bad,ok=False)
    assert not (cache/"initial").exists() and not (cache/"observer").exists()
    rebuilt=ROOT/"outputs/v296_validation/install/outputs/thermistor_inputs/B6-56/observer/observe_bead.exe"
    invoke("installThermistorObserver.py",rebuilt)
    invoke("thermistorStartupObserver.py")
    assert not (install/"tools/runtime/su2-thermal-build").exists()
    info=json.loads((cache/"observer/build.json").read_bytes())
    assert info["origin"]=="verified-installed-binary" and info["observerSHA256"]==sha(rebuilt)==sha(cache/"observer/observe_bead.exe")
    return dict(originalManifestSHA256=shared,relocatedManifestSHA256=new,sourceManifestUnchanged=True,
                secondInstallIdempotent=True,incorrectInitialStateRejected=True,incorrectObserverRejected=True,
                observerRestoredWithoutCompiler=True,observerSHA256=info["observerSHA256"])


if __name__=="__main__":
    proof=run()
    with Path(sys.argv[1]).open("x",encoding="utf8") as f:
        json.dump(proof,f,ensure_ascii=False,indent=2)
