"""监督一次MATLAB只读入口检查；超时只结束本次子进程树，不重启。"""
from pathlib import Path
import json,subprocess,sys,time
import psutil
ROOT=Path(__file__).resolve().parents[2];base=Path(__file__).resolve().parent
receipt=base/'matlab_parent.json';assert not receipt.exists()
argv=[sys.executable,str(ROOT/'scripts/runMatlabResearch.py'),'explicit01',
      "addpath('tests'); r=testEilmerExplicitEntry(); assert(r.passed);",'--folder','outputs/v183_validation',
      '--input','outputs/v183_validation/small01/case_manifest.json',
      '--input','outputs/v183_validation/small01_recovery/report.json',
      '--input','outputs/v183_validation/small01_run/report.json']
start=time.monotonic();state=dict(status='running',secondsLimit=90,argv=argv)
with (base/'matlab_parent_stdout.txt').open('wb') as log:
    p=subprocess.Popen(argv,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
    state['pid']=p.pid;receipt.write_text(json.dumps(state,indent=2))
    try:state['exitCode']=p.wait(timeout=90)
    except subprocess.TimeoutExpired:
        state['timeout']=True
        try:
            group=psutil.Process(p.pid).children(recursive=True)
            for child in reversed(group):
                try:child.kill()
                except psutil.NoSuchProcess:pass
        except psutil.NoSuchProcess:pass
        p.kill();state['exitCode']=p.wait(timeout=10)
state.update(status='terminal',elapsedSeconds=time.monotonic()-start)
receipt.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state))
raise SystemExit(0 if state['exitCode']==0 and not state.get('timeout') else 1)
