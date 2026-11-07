"""OpenFOAM原生任务监督：固定包、进程组、输入、容量与温域警告。"""
from pathlib import Path
import json,os,resource,signal,subprocess,time,shutil
from runEilmerProcess import group_members,sha

ROOT=Path(__file__).resolve().parents[1]
RUNTIME=ROOT/'tools/runtime/openfoam/13-20260624/runtime.json'


def run(argv,cwd,record,seconds,locks=(),output_mib=96,memory_mib=2048,guard=None):
    """模块一：只监督自己创建的进程组，不改变共享环境；过程成功不授物理资格。"""
    if os.name!='posix' or not 0<seconds<=900:raise ValueError('仅限WSL有限任务')
    if not 0<output_mib<=512 or not 0<memory_mib<=4096:raise ValueError('资源超出入口范围')
    cwd,record=Path(cwd).resolve(),Path(record).resolve()
    record.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,record/'executed_runner.py')
    runtime=json.loads(RUNTIME.read_text())
    assert runtime['profile']=='openfoam13-20260624-stock'
    for p,h in runtime['filesSHA256'].items():assert sha(p)==h,p
    inputs={str(Path(p).resolve()):sha(p) for p in [__file__,RUNTIME,*locks]}
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',LC_ALL='C')
    def limits():
        resource.setrlimit(resource.RLIMIT_AS,(memory_mib*1024**2,)*2)
        resource.setrlimit(resource.RLIMIT_CORE,(0,0))
        os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:1])
    # 固定脚本只加载原生环境；用户参数通过argv传递，不插入shell代码。
    command=['bash','--noprofile','--norc','-c','source /opt/openfoam13/etc/bashrc; exec "$@"','rah-openfoam',*map(str,argv)]
    start=time.monotonic();proc=None;stop=None;peak=0;code=None;extra={}
    log=record/'stdout.txt'
    def terminate():
        if proc is None:return
        try:os.killpg(proc.pid,signal.SIGTERM)
        except ProcessLookupError:return
        try:proc.wait(timeout=2)
        except subprocess.TimeoutExpired:pass
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        proc.wait(timeout=2)
    def cancelled(sig,frame):raise InterruptedError('supervisor signal '+str(sig))
    old=signal.signal(signal.SIGTERM,cancelled)
    try:
        with log.open('wb') as out:
            proc=subprocess.Popen(command,cwd=cwd,env=env,stdout=out,stderr=subprocess.STDOUT,
                start_new_session=True,preexec_fn=limits)
            next_size=0.;size=0
            while proc.poll() is None:
                elapsed=time.monotonic()-start
                members=group_members(proc.pid);rss=sum(x['rssBytes'] for x in members.values());peak=max(peak,rss)
                if elapsed>=next_size:
                    size=sum(p.stat().st_size for p in cwd.rglob('*') if p.is_file());next_size=elapsed+2
                if elapsed>seconds:stop='wall-clock-limit'
                if rss>memory_mib*1024**2:stop='group-rss-limit'
                if size>output_mib*1024**2:stop='case-output-limit'
                if stop is None and guard is not None:
                    problem=guard(log,False)
                    if problem:stop='scientific-guard';extra['guardReason']=problem
                if stop:terminate();break
                time.sleep(.15)
            code=proc.wait()
            if stop is None and guard is not None:
                problem=guard(log,True)
                if problem:stop='scientific-guard-after-exit';extra['guardReason']=problem
    except BaseException as error:
        stop='supervisor-interrupted' if isinstance(error,(KeyboardInterrupt,InterruptedError)) else 'supervisor-error'
        extra['supervisorError']=repr(error)
    finally:
        terminate();signal.signal(signal.SIGTERM,old)
        if proc is not None:code=proc.returncode
    live=group_members(proc.pid) if proc is not None else {}
    end=time.monotonic()+2
    while live and time.monotonic()<end:time.sleep(.05);live=group_members(proc.pid)
    stable={p:Path(p).is_file() and sha(p)==h for p,h in inputs.items()}
    size=sum(p.stat().st_size for p in cwd.rglob('*') if p.is_file())
    result=dict(schema='RAH.OpenFOAMProcess.v1',argv=command,cwd=str(cwd),nativeExitCode=code,
        stopReason=stop,elapsedSeconds=time.monotonic()-start,maximumSeconds=seconds,
        peakObservedGroupRSSBytes=peak,liveGroupAfterCleanup=live,inputSHA256=inputs,inputStable=stable,
        caseBytes=size,caseLimitMiB=output_mib,addressSpaceLimitMiB=memory_mib,cpuCount=1,
        physicalQualification=False,**extra)
    result['passed']=code==0 and stop is None and not live and all(stable.values()) and size<=output_mib*1024**2
    (record/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
