"""模块：只监督本轮两个固定参考任务，限制时间/内存/输出，不自动重试。"""
from pathlib import Path
import json,os,subprocess,sys,time,shutil
import psutil
ROOT=Path(__file__).resolve().parents[1]
REFERENCE_PYTHON=Path.home()/'.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'

def main():
    base=ROOT/'outputs/v204_validation';out=base/'reference_supervision01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');started=time.monotonic();receipts=[]
    for kind in ['mixture','air']:
        command=[str(REFERENCE_PYTHON),str(ROOT/'scripts/validateNistPressureDomain.py'),'--kind',kind,'--out',str(base/(kind+'01'))]
        env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');start=time.monotonic();peak=0;stop=None
        with (out/(kind+'.log')).open('wb') as f:
            p=subprocess.Popen(command,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,env=env,creationflags=subprocess.CREATE_NO_WINDOW)
            while p.poll() is None:
                try:
                    family=[psutil.Process(p.pid)];family+=family[0].children(recursive=True)
                    rss=sum(x.memory_info().rss for x in family if x.is_running());peak=max(peak,rss)
                except psutil.NoSuchProcess:pass
                size=sum(x.stat().st_size for x in base.rglob('*') if x.is_file() and (kind+'01' in x.parts or out in x.parents))
                if time.monotonic()-start>60 or time.monotonic()-started>120 or peak>1024**3 or size>16*1024**2:
                    stop='time-memory-output-budget'
                    for child in reversed(family):
                        try:child.kill()
                        except psutil.NoSuchProcess:pass
                    break
                time.sleep(.05)
            code=p.wait()
        receipt=dict(command=command,nativeExitCode=code,stopReason=stop,peakObservedRSSBytes=peak,elapsedSeconds=time.monotonic()-start,
            maximumSeconds=60,totalMaximumSeconds=120,memoryBytes=1024**3,outputBytes=16*1024**2)
        receipts.append(receipt);(out/'processes.json').write_text(json.dumps(receipts,indent=2)+'\n',encoding='utf-8');print(json.dumps(receipt),flush=True)
        if code!=0 or stop:raise RuntimeError('参考任务失败，保留输出，不自动重复')
    (out/'report.json').write_text(json.dumps(dict(passed=True,processes=receipts,elapsedSeconds=time.monotonic()-started),indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
