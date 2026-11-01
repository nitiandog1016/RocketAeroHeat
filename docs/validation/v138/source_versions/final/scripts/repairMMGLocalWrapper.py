"""仅修正未编入SCOTCH时的参数调用；不重编上游库、不重做旧数值。"""
from pathlib import Path
import json
from prepareCGALRemeshingRuntime import ROOT,sha
from prepareMMGLocalRuntime import DEST,EVID
from buildCGALRemeshingRuntime import run

def main():
    old=json.loads((DEST/'runtime.json').read_text());compiler=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe'
    src=ROOT/'third_party/mmg/local_tetra_remesh.cpp';lib=DEST/'build/lib/libmmg3d.a';exe=DEST/'local_tetra_remesh.exe'
    assert sha(lib)==old['files']['build/lib/libmmg3d.a']['sha256']
    locked={p.resolve().as_posix():sha(p) for p in [src,lib,compiler,Path(__file__).resolve()]}
    r=run([str(compiler),'-std=c++17','-O2','-static','-I'+str(DEST/'build/include'),str(src),str(lib),'-o',str(exe)],EVID/'wrapper02',60,locked)
    assert r['exitCode']==0 and all(r['sourceStableAfter'].values())
    files={p.relative_to(DEST).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in DEST.rglob('*') if p.is_file() and p.name!='runtime.json'}
    old.update(executableSHA256=sha(exe),wrapperCorrection=r,files=files,fileLengthBytes=sum(v['bytes'] for v in files.values()),
        originalBuildSourceSHA256=old['sourceSHA256'],sourceSHA256=dict(old['sourceSHA256'],**{src.resolve().as_posix():sha(src)}),
        correctionReason='renum parameter exists only with USE_SCOTCH; runtime built without SCOTCH',
        priorRuntimeSHA256=sha(EVID/'runtime_before_interface_fix.json'))
    assert old['fileLengthBytes']<350*1024**2
    (DEST/'runtime.json').write_text(json.dumps(old,indent=2)+'\n');(EVID/'runtime.json').write_bytes((DEST/'runtime.json').read_bytes())
    print(r['exitCode'],r['seconds'],old['fileLengthBytes'])
if __name__=='__main__':main()
