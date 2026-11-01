"""有预算地构建固定CGAL局部重网格入口，保存进程与准确源码身份。"""
from pathlib import Path
import json,os,subprocess,time
import psutil
from prepareCGALRemeshingRuntime import ROOT,DEST,EVID,sha

def run(command,folder,timeout,locked):
    folder.mkdir(parents=True,exist_ok=False);start=time.monotonic();peak=0;stopped=None
    with (folder/'stdout.log').open('wb') as log:
        p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
        while p.poll() is None:
            try:
                proc=psutil.Process(p.pid);family=[proc]+proc.children(recursive=True)
                memory=sum(q.memory_info().rss for q in family if q.is_running());peak=max(peak,memory)
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
            if time.monotonic()-start>timeout or peak>6*1024**3:
                stopped='time' if time.monotonic()-start>timeout else 'memory'
                for q in reversed(family):
                    try:q.kill()
                    except psutil.NoSuchProcess:pass
                p.wait();break
            time.sleep(.2)
    receipt=dict(command=command,exitCode=p.returncode,seconds=time.monotonic()-start,peakFamilyRSSBytes=peak,
        stopReason=stopped,sourceSHA256=locked,sourceStableAfter={n:sha(ROOT/n)==h for n,h in locked.items()},
        outputSHA256={q.name:sha(q) for q in folder.iterdir() if q.is_file()},status='terminal')
    (folder/'process.json').write_text(json.dumps(receipt,indent=2)+'\n');return receipt

def main():
    meta=json.loads((DEST/'headers.json').read_text());assert meta['status']=='headers-prepared'
    for name,entry in meta['files'].items():assert sha(DEST/name)==entry['sha256']
    compiler=Path(meta['compiler']);assert sha(compiler)==meta['compilerSHA256']
    src=ROOT/'third_party/cgal/local_tetra_remesh.cpp';exe=DEST/'local_tetra_remesh.exe'
    command=[str(compiler),'-std=c++17','-O1','-frounding-math','-DCGAL_DISABLE_GMP=1','-static',
        '-I'+str(DEST/'cgal/include'),'-I'+str(DEST/'boost'),str(src),'-o',str(exe)]
    locked={p.relative_to(ROOT).as_posix():sha(p) for p in [src,Path(__file__).resolve(),compiler,DEST/'headers.json']}
    r=run(command,EVID/'build01',240,locked);print(json.dumps(r,indent=2))
    assert r['exitCode']==0 and all(r['sourceStableAfter'].values())
    for name,entry in meta['files'].items():assert sha(DEST/name)==entry['sha256']
    meta.update(status='built-not-physics-qualified',executable=str(exe),executableSHA256=sha(exe),
        sourceSHA256=locked,buildReceipt=r,filesStableAfter=True,numberOfIterations=1,physicsQualified=False)
    (DEST/'runtime.json').write_text(json.dumps(meta,indent=2)+'\n')
    (EVID/'runtime.json').write_bytes((DEST/'runtime.json').read_bytes())
if __name__=='__main__':main()
