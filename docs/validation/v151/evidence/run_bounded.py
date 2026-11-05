"""本轮命令监督：保存真实终态，限制墙钟和本进程族内存，不操作其他任务。"""
import argparse,hashlib,json,os,subprocess,time
from pathlib import Path
import psutil

root=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('request');args=ap.parse_args()
    request_path=Path(args.request).resolve();request=json.loads(request_path.read_text(encoding='utf-8'))
    folder=root/request['receiptFolder'];folder.resolve().relative_to(root/'outputs/v151_validation');assert not folder.exists();folder.mkdir()
    locked={str(p):sha(p) for p in [request_path,Path(__file__).resolve()]+[root/n for n in request.get('sources',[])]}
    start=time.monotonic();peak=0;reason=None
    with (folder/'stdout.log').open('wb') as log:
        p=subprocess.Popen(request['command'],cwd=str(root),stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'))
        (folder/'process.json').write_text(json.dumps(dict(status='running',processID=p.pid,command=request['command'])))
        while p.poll() is None:
            try:
                family=[psutil.Process(p.pid)]+psutil.Process(p.pid).children(recursive=True)
                rss=sum(q.memory_info().rss for q in family if q.is_running());peak=max(peak,rss)
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
            if time.monotonic()-start>request['maximumSeconds']:reason='wall-clock budget'
            if peak>request.get('maximumRSSBytes',8*1024**3):reason='family memory budget'
            if reason:
                if p.poll() is None:subprocess.run(['taskkill','/PID',str(p.pid),'/T','/F'],stdout=log,stderr=log,creationflags=subprocess.CREATE_NO_WINDOW)
                break
            time.sleep(.25)
        code=p.wait()
    stable={path:sha(Path(path))==h for path,h in locked.items()}
    result=dict(status='terminal',exitCode=code,seconds=time.monotonic()-start,peakFamilyRSSBytes=peak,
        stopReason=reason,command=request['command'],sourceSHA256=locked,sourceStableAfter=stable,
        stdoutSHA256=sha(folder/'stdout.log'))
    (folder/'process.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','exitCode','seconds','peakFamilyRSSBytes','stopReason']}))
    return code if code else (0 if all(stable.values()) and reason is None else 1)
if __name__=='__main__':raise SystemExit(main())
