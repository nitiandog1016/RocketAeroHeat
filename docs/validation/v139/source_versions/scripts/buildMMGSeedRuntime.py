"""复用准确MMG静态库构建保留内部点的独立入口，不覆写v138运行时。"""
from pathlib import Path
import json,shutil
from prepareCGALRemeshingRuntime import ROOT,sha
from prepareMMGLocalRuntime import DEST as BASE
from buildCGALRemeshingRuntime import run
DEST=ROOT/'tools/runtime/mmg/5.8.0-seeds'

def main():
    assert not DEST.exists();DEST.mkdir(parents=True)
    evidence=ROOT/'docs/validation/v139/runtime';evidence.mkdir(parents=True,exist_ok=True)
    old=json.loads((BASE/'runtime.json').read_text());lib=BASE/'build/lib/libmmg3d.a'
    assert sha(lib)==old['files']['build/lib/libmmg3d.a']['sha256']
    compiler=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe'
    src=ROOT/'third_party/mmg/local_tetra_remesh.cpp';exe=DEST/'local_tetra_optimize.exe'
    paths=[lib,compiler,src,Path(__file__).resolve(),BASE/'runtime.json']+list((BASE/'build/include').rglob('*'))
    locked={str(p.resolve()):sha(p) for p in paths if p.is_file()}
    r=run([str(compiler),'-std=c++17','-O2','-static','-DRAH_MMG_KEEP_SEEDS=1','-I'+str(BASE/'build/include'),str(src),str(lib),'-o',str(exe)],evidence/'build',60,locked)
    assert r['exitCode']==0 and all(r['sourceStableAfter'].values())
    shutil.copyfile(str(src),str(DEST/'local_tetra_remesh.cpp'))
    record=dict(version='5.8.0',flavor='fixed-boundary-keep-seeds',executable=str(exe),executableSHA256=sha(exe),
        upstreamLibrary=str(lib),upstreamLibrarySHA256=sha(lib),sourceSHA256=locked,build=r,
        noinsert=True,nomove=False,noswap=False,physicsQualified=False,
        fileLengthBytes=sum(p.stat().st_size for p in DEST.iterdir() if p.is_file()))
    assert record['fileLengthBytes']<15*1024**2
    (DEST/'runtime.json').write_text(json.dumps(record,indent=2)+'\n');(evidence/'runtime.json').write_bytes((DEST/'runtime.json').read_bytes())
    print(r['exitCode'],r['seconds'],record['fileLengthBytes'])
if __name__=='__main__':main()
