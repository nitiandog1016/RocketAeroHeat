"""模块职责：监督既有数据的有限只读设计/判断；不包含原生/网格/CFD入口。"""
from pathlib import Path
import json,sys,os,time,subprocess,hashlib,shutil
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v216_validation'
ENTRIES={'inventory':'inspectCompleteRegionInputsV216.py','design':'designCompleteRegionsV216.py',
         'numerical':'assessFrozenDifferenceResolutionV216.py','final':'finalizeCompleteRegionDecisionV216.py'}
def windows_run(argv,record,seconds,locks):
    """Windows复用已有SciPy；采样实际RSS/虚拟字节并回收本任务进程树，无硬RLIMIT声明。"""
    import psutil
    record.mkdir(exist_ok=False);hashfile=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    before={str(p):hashfile(p) for p in locks};shutil.copyfile(__file__,record/'executed_entry.py')
    shutil.copyfile(argv[-1],record/'executed_task.py')
    start=time.monotonic();peak=0;vpeak=0;stop=None;family={}
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONIOENCODING='utf8')
    with (record/'stdout.txt').open('wb') as log:
        p=subprocess.Popen(argv,cwd=str(BASE),stdout=log,stderr=log,env=env,creationflags=subprocess.CREATE_NO_WINDOW)
        parent=psutil.Process(p.pid);parent.cpu_affinity(parent.cpu_affinity()[:1]);family[p.pid]=parent
        try:
            while p.poll() is None:
                try:
                    for child in parent.children(recursive=True):family[child.pid]=child
                except psutil.NoSuchProcess:pass
                rss=0;vms=0
                for child in family.values():
                    try:
                        m=child.memory_info();rss+=m.rss;vms+=m.vms
                    except psutil.NoSuchProcess:pass
                peak=max(peak,rss);vpeak=max(vpeak,vms)
                size=sum(q.stat().st_size for q in BASE.rglob('*') if q.is_file())
                if time.monotonic()-start>seconds or max(rss,vms)>4*1024**3 or size>128*1024**2:
                    stop='resource-limit';break
                time.sleep(.05)
        finally:
            if p.poll() is None:
                for child in reversed(list(family.values())):
                    try:child.kill()
                    except psutil.NoSuchProcess:pass
            code=p.wait(timeout=5);_,alive=psutil.wait_procs(list(family.values()),timeout=5)
    stable=all(hashfile(p)==h for p,h in before.items())
    r=dict(status='terminal',argv=argv,passed=code==0 and stop is None and stable and not alive,nativeExitCode=code,
        elapsedSeconds=time.monotonic()-start,peakObservedGroupRSSBytes=peak,peakObservedGroupVirtualBytes=vpeak,
        inputsStable=stable,inputSHA256=before,stopReason=stop,liveGroupAfterCleanup=[p.pid for p in alive],
        cpuCountLimit=1,limitSeconds=seconds,memoryMiB=4096,hardAddressSpaceLimit=False,CFDAdvanceRuns=0)
    (record/'report.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf8');return r
def main(stage):
    BASE.mkdir(exist_ok=True);used=0
    for p in BASE.glob('process_*/report.json'):
        r=json.loads(p.read_text());assert r['status']!='running';used+=r['elapsedSeconds']
    name='inventory' if stage in ('inventory2','inventory3') else stage
    assert name in ENTRIES and used<120
    source=ROOT/'scripts'/ENTRIES[name]
    locks=[source,Path(__file__).resolve(),ROOT/'docs/validation/v216/PLAN.md',
        ROOT/'scripts/eilmerNKSAResourceBudget.py',ROOT/'outputs/v144_validation/actual/hybrid_volume.mat',
        ROOT/'outputs/v212_validation/design01/regions.json',ROOT/'outputs/v214_validation/control01/nk_samples.csv',
        ROOT/'outputs/v214_validation/control01/directions.csv',ROOT/'outputs/v213_validation/control03/geometry.csv']
    if os.name=='nt':r=windows_run([sys.executable,str(source)],BASE/('process_'+stage),min(90,120-used),locks)
    else:
        from runEilmerProcess import run
        r=run(['taskset','-c',str(min(os.sched_getaffinity(0))),sys.executable,str(source)],BASE,BASE/('process_'+stage),
            min(90,120-used),128,locks,4096,4096)
    print(json.dumps({k:r[k] for k in ['passed','elapsedSeconds','peakObservedGroupRSSBytes','nativeExitCode']}))
    return r['passed']
if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
