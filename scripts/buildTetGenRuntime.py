"""限时构建TetGen原版数值核心与独立接口，记录每一步进程和准确字节。"""
from pathlib import Path
import json,time
from prepareTetGenRuntime import ROOT,DEST,EVID,sha
from buildCGALRemeshingRuntime import run

def main():
    start=time.monotonic();meta=json.loads((DEST/'source.json').read_text())
    for n,r in meta['sourceGitBlobChecks'].items():assert sha(DEST/'source'/n)==r['sha256']
    cc=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe'
    wrapper=ROOT/'third_party/tetgen/fixed_tetra_repair.cpp'
    paths=[cc,wrapper,Path(__file__).resolve(),ROOT/'scripts/buildCGALRemeshingRuntime.py',DEST/'source.json']
    paths.extend((DEST/'source').iterdir());locked={str(p.resolve()):sha(p) for p in paths if p.is_file()}
    common=[str(cc),'-std=c++11','-DTETLIBRARY','-I'+str(DEST/'source')]
    commands=[common+['-O0','-c',str(DEST/'source/predicates.cxx'),'-o',str(DEST/'predicates.o')],
        common+['-O2','-c',str(DEST/'source/tetgen.cxx'),'-o',str(DEST/'tetgen.o')],
        common+['-O2','-static',str(wrapper),str(DEST/'tetgen.o'),str(DEST/'predicates.o'),'-o',str(DEST/'fixed_tetra_repair.exe')]]
    receipts=[]
    for i,cmd in enumerate(commands):
        remaining=120-(time.monotonic()-start);assert remaining>0
        r=run(cmd,EVID/('build%02d'%(i+1)),remaining,locked);receipts.append(r)
        print(i+1,r['exitCode'],r['seconds'],flush=True)
        assert r['exitCode']==0 and all(r['sourceStableAfter'].values())
    exe=DEST/'fixed_tetra_repair.exe';size=sum(p.stat().st_size for p in DEST.rglob('*') if p.is_file());assert size<=30*1024**2
    result=dict(status='built-not-mesh-or-physics-qualified',version='1.5.1',executable=str(exe),
        executableSHA256=sha(exe),sourceSHA256=locked,receipts=receipts,seconds=time.monotonic()-start,
        runtimeFileLengthBytes=size,license='AGPL-3.0-or-later',physicsQualified=False)
    assert result['seconds']<=120
    for p in [DEST/'runtime.json',EVID/'runtime.json']:p.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
