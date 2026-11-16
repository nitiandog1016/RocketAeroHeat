"""模块职责：只编译单壳零孔薄入口并链接准确TetGen旧对象，不重编成熟核心。"""
import subprocess
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v242_validation';DEST=ROOT/'tools/runtime/tetgen/1.5.1'


def build():
    assert read(BASE/'surface01/report.json')['passed']
    expected=read(ROOT/'docs/validation/v142/storage_and_sources.json')['runtimeFiles']
    for name in ('tetgen.o','predicates.o','source/tetgen.h'):
        assert sha(DEST/name)==expected[name]['sha256'],name
    out=BASE/'build01';out.mkdir(exist_ok=False)
    cc=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe';src=ROOT/'third_party/tetgen/single_shell_plc_v242.cpp';exe=out/'single_shell_plc.exe'
    cmd=[str(cc),'-std=c++11','-O2','-DTETLIBRARY','-I'+str(DEST/'source'),'-static',str(src),str(DEST/'tetgen.o'),str(DEST/'predicates.o'),'-o',str(exe)]
    save(out/'command.json',cmd);result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (out/'build.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout
    save(out/'runtime.json',dict(passed=True,version='1.5.1',license='AGPL-3.0-or-later',executable=str(exe),executableSHA256=sha(exe),
      numericalCoreObjectsRecompiled=False,sourceSHA256={str(p):sha(p) for p in [src,DEST/'source/tetgen.h',DEST/'tetgen.o',DEST/'predicates.o',cc]},physicalQualification=False))


if __name__=='__main__':build()
