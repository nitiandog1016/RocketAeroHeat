"""Linux/WSL内监督一个Eilmer过程；固定预算、进程组、输入SHA与准确失败记录。"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import resource
import shutil
import signal
import subprocess
import time

PREFIX = Path('/opt/rocket-aeroheat/eilmer-v5.0.0')
COMMIT = 'f53f4609a0331d48efee69a4e4f3c3598378cc03'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def run(argv, cwd, record, seconds, output_mib=128, locks=(), memory_mib=8192):
    """模块一：明确源与运行包；输出只写新回执目录，绝不隐式接续旧求解。"""
    assert os.name == 'posix' and 0 < seconds <= 1800
    assert 0 < output_mib <= 2048 and 0 < memory_mib <= 8192
    cwd, record = Path(cwd).resolve(), Path(record).resolve()
    assert cwd.is_dir() and not record.exists()
    record.mkdir(parents=True)
    shutil.copyfile(__file__, record/'executed_runner.py')
    head = subprocess.check_output(['git','-C',str(PREFIX/'source'),'rev-parse','HEAD'], text=True).strip()
    assert head == COMMIT
    subprocess.run(['git','-C',str(PREFIX/'source'),'diff','--exit-code'], check=True, stdout=subprocess.DEVNULL)
    hashes = {str(Path(p).resolve()):sha(p) for p in locks}
    hashes[str(Path(__file__).resolve())] = sha(__file__)
    runtime = {p.name:sha(p) for p in (PREFIX/'install/bin').glob('lmr*') if p.is_file()}
    env = os.environ.copy()
    env.update(DGD=str(PREFIX/'install'), DGD_REPO=str(PREFIX/'source'),
               DGD_LUA_PATH=str(PREFIX/'install/lib/?.lua'),
               DGD_LUA_CPATH=str(PREFIX/'install/lib/?.so'),
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1')
    env['PATH'] = str(PREFIX/'install/bin')+os.pathsep+env['PATH']
    result = dict(schema='RAH.EilmerProcess.v1', argv=argv, cwd=str(cwd), sourceCommit=head,
                  runtimeSHA256=runtime, limitSeconds=seconds, outputLimitMiB=output_mib,
                  addressSpaceLimitMiB=memory_mib, cpuCountLimit=2, physicalQualification=False,
                  status='running', startedUTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))

    def limit_child():
        resource.setrlimit(resource.RLIMIT_AS, (memory_mib*1024**2, memory_mib*1024**2))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
        if hasattr(os, 'sched_getaffinity'):
            os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:2])

    def save():
        (record/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')

    save(); started=time.monotonic(); peak=0; bytes_now=0; stopped=None
    # 模块二：同组终止与容量监测；预算退出不同于原生正常退出。
    proc=None; code=None
    def stop_group():
        if proc is None:return
        try:os.killpg(proc.pid,signal.SIGTERM)
        except ProcessLookupError:return
        try:proc.wait(timeout=3)
        except subprocess.TimeoutExpired:pass
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        proc.wait(timeout=3)
    def interrupted(signum,frame):
        raise InterruptedError('supervisor signal '+str(signum))
    old_term=signal.signal(signal.SIGTERM,interrupted)
    try:
      with (record/'stdout.txt').open('wb') as out:
        proc=subprocess.Popen(argv,cwd=cwd,env=env,stdout=out,stderr=subprocess.STDOUT,
                              start_new_session=True,preexec_fn=limit_child)
        result['processId']=proc.pid;save();next_size=0.0
        while proc.poll() is None:
            elapsed=time.monotonic()-started
            if elapsed >= next_size:
                bytes_now=sum(p.stat().st_size for p in cwd.rglob('*') if p.is_file())
                next_size=elapsed+2
            try:
                for line in Path('/proc/{}/status'.format(proc.pid)).read_text().splitlines():
                    if line.startswith('VmRSS:'): peak=max(peak,int(line.split()[1])*1024)
            except (OSError,ValueError):
                pass
            if elapsed > seconds: stopped='wall-clock-limit'
            if bytes_now > output_mib*1024**2: stopped='output-size-limit'
            if stopped:
                stop_group()
                break
            time.sleep(.2)
        code=proc.wait()
    except BaseException as error:
        stopped='supervisor-interrupted' if isinstance(error,(KeyboardInterrupt,InterruptedError)) else 'supervisor-error'
        result['supervisorError']=repr(error)
    finally:
        # 外层取消、读取异常和父过程退出时都收回自己的进程组，不留下后台求解。
        stop_group()
        if proc is not None:code=proc.returncode
        signal.signal(signal.SIGTERM,old_term)
    stable={p:Path(p).is_file() and sha(p)==h for p,h in hashes.items()}
    final_bytes=sum(p.stat().st_size for p in cwd.rglob('*') if p.is_file())
    # 模块三：退出0只授予过程完成，不授予任何数值或物理资格。
    result.update(nativeExitCode=code,stopReason=stopped,elapsedSeconds=time.monotonic()-started,
                  peakObservedParentRSSBytes=peak,finalCaseBytes=final_bytes,
                  inputSHA256=hashes,inputStable=stable,inputsStable=all(stable.values()))
    result['passed']=code==0 and stopped is None and all(stable.values()) and final_bytes<=output_mib*1024**2
    result['status']='process-completed' if result['passed'] else 'failed'
    save()
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cwd',required=True);p.add_argument('--record',required=True)
    p.add_argument('--seconds',type=float,required=True)
    p.add_argument('--output-mib',type=int,default=128)
    p.add_argument('--lock',action='append',default=[])
    p.add_argument('argv',nargs=argparse.REMAINDER)
    a=p.parse_args();argv=a.argv[1:] if a.argv[:1]==['--'] else a.argv
    assert argv
    r=run(argv,a.cwd,a.record,a.seconds,a.output_mib,a.lock)
    print(json.dumps({k:r[k] for k in ('passed','nativeExitCode','stopReason','elapsedSeconds','inputsStable')}))
    return 0 if r['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
