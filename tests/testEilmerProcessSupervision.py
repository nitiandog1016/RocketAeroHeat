"""实际创建有子进程的小任务，核对超时与取消均收回本任务进程组。"""
from pathlib import Path
import argparse
import json
import os
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import run,sha


def verify(out):
    start=time.monotonic();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False)
    payload="import subprocess,time; p=subprocess.Popen(['sleep','60']); print(p.pid,flush=True); time.sleep(60)"
    r=run([sys.executable,'-c',payload],out,out/'timeout',.8,128)
    assert not r['passed'] and r['stopReason']=='wall-clock-limit'
    child=int((out/'timeout/stdout.txt').read_text().strip())
    def running(pid):
        try:return Path('/proc/%d/stat'%pid).read_text().split()[2] not in ('Z','X')
        except FileNotFoundError:return False
    assert not running(child) and not running(r['processId'])
    command=[sys.executable,str(ROOT/'scripts/runEilmerProcess.py'),'--cwd',str(out),
             '--record',str(out/'cancel'),'--seconds','30','--',sys.executable,'-c',payload]
    p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    deadline=time.monotonic()+8
    while time.monotonic()<deadline:
        f=out/'cancel/stdout.txt'
        if f.is_file() and f.stat().st_size:break
        time.sleep(.05)
    else:p.kill();raise AssertionError('测试子任务没有启动')
    child2=int(f.read_text().strip());p.send_signal(signal.SIGTERM)
    stdout=p.communicate(timeout=8)[0];(out/'cancel_launcher.txt').write_bytes(stdout)
    cancel=json.loads((out/'cancel/report.json').read_text())
    assert not cancel['passed'] and cancel['stopReason']=='supervisor-interrupted'
    assert not running(child2) and not running(cancel['processId'])
    report=dict(passed=True,timeoutProcessGroupStopped=True,cancelProcessGroupStopped=True,
                nativeCFDRuns=0,elapsedSeconds=time.monotonic()-start,
                runnerSHA256=sha(ROOT/'scripts/runEilmerProcess.py'))
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True)
    print(json.dumps(verify(p.parse_args().out)))
