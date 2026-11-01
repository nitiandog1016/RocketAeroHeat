"""复用准确核心对象，构建独立PLC入口；不重编/修改已交付TetGen数值源码。"""
from pathlib import Path
import json
from prepareTetGenRuntime import ROOT,DEST,sha
from buildCGALRemeshingRuntime import run
PLC=ROOT/'tools/runtime/tetgen/plc-v143'

def main():
    assert not PLC.exists();PLC.mkdir(parents=True)
    previous=json.loads((ROOT/'docs/validation/v142/storage_and_sources.json').read_text())['runtimeFiles']
    for n,r in previous.items():assert sha(DEST/n)==r['sha256']
    cc=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe';src=ROOT/'third_party/tetgen/fixed_surface_plc.cpp'
    paths=[cc,src,Path(__file__).resolve(),ROOT/'scripts/buildCGALRemeshingRuntime.py',DEST/'source/tetgen.h',
        DEST/'tetgen.o',DEST/'predicates.o',DEST/'runtime.json',ROOT/'docs/validation/v143/plan.md']
    locked={str(p.resolve()):sha(p) for p in paths};exe=PLC/'fixed_surface_plc.exe'
    cmd=[str(cc),'-std=c++11','-O2','-DTETLIBRARY','-I'+str(DEST/'source'),'-static',str(src),
        str(DEST/'tetgen.o'),str(DEST/'predicates.o'),'-o',str(exe)]
    r=run(cmd,ROOT/'docs/validation/v143/build',60,locked)
    assert r['exitCode']==0 and all(r['sourceStableAfter'].values())
    meta=dict(version='1.5.1',status='PLC-built-not-tested',executable=str(exe),executableSHA256=sha(exe),
        sourceSHA256=locked,buildReceipt=r,physicsQualified=False,license='AGPL-3.0-or-later')
    for dst in [PLC/'runtime.json',ROOT/'docs/validation/v143/runtime.json']:dst.write_text(json.dumps(meta,indent=2)+'\n')
    size=sum(p.stat().st_size for p in PLC.iterdir());assert size<10*1024**2
    print(r['exitCode'],r['seconds'],size,flush=True)

if __name__=='__main__':main()
