"""隐藏、有限地核对新MATLAB物性接口，准确源在执行前绑定。"""
from pathlib import Path
import json,subprocess,time,shutil,os
import psutil
from prepareOpenFoamControl import ROOT,sha,save

def main():
    out=ROOT/'outputs/v206_validation/matlab01';out.mkdir(exist_ok=False);start=time.monotonic()
    names=['phaseAwareAirState.m','phaseAwareAirReferenceState.m','tests/testPhaseAwareAirState.m',
        'data/thermo/phase_aware_dry_air_v1.json','outputs/v206_validation/python_tests01/matlab_fixture.json',
        'docs/validation/v206/PLAN.md','scripts/checkPhaseAwareMatlab.py']
    hashes={n:sha(ROOT/n) for n in names}
    for n in names:
        q=out/'executed_source'/n;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,q)
    save(out/'admission.json',dict(maximumSeconds=60,memoryMiB=2048,outputMiB=32,inputSHA256=hashes,automaticRetry=0))
    command=['D:/matlab/bin/matlab.exe','-batch',"restoredefaultpath; maxNumCompThreads(1); addpath('tests'); testPhaseAwareAirState();"]
    family=[];peak=0;stop=None
    with (out/'stdout.log').open('wb') as log:
        p=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW,
            env=dict(os.environ,MATLABPATH='',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        try:
            while p.poll() is None:
                try:
                    parent=psutil.Process(p.pid);family=[parent]+parent.children(recursive=True)
                    peak=max(peak,sum(x.memory_info().rss for x in family if x.is_running()))
                except psutil.NoSuchProcess:pass
                size=sum(x.stat().st_size for x in out.rglob('*') if x.is_file())
                if time.monotonic()-start>60 or peak>2048*1024**2 or size>32*1024**2:
                    stop='resource-limit';break
                time.sleep(.1)
        finally:
            if p.poll() is None:
                for child in reversed(family):
                    try:child.kill()
                    except psutil.NoSuchProcess:pass
            code=p.wait(timeout=5)
    live=[]
    for child in family:
        try:
            if child.is_running() and child.status()!=psutil.STATUS_ZOMBIE:live.append(child.pid)
        except psutil.NoSuchProcess:pass
    stable={n:sha(ROOT/n)==h for n,h in hashes.items()}
    report=dict(passed=code==0 and stop is None and not live and all(stable.values()),exitCode=code,stop=stop,
        peakObservedRSSBytes=peak,liveAfter=live,inputStable=stable,inputSHA256=hashes,elapsedSeconds=time.monotonic()-start)
    save(out/'process.json',report);print(json.dumps(report));return report['passed']

if __name__=='__main__':raise SystemExit(0 if main() else 1)
