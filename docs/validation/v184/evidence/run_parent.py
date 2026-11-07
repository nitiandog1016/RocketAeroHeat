"""记录独立父时间；timeout负责有限终止，内部监督负责收回自有求解进程组。"""
from pathlib import Path
import json,subprocess,sys,time
base=Path(__file__).resolve().parent;receipt=base/'parent.json';assert not receipt.exists()
argv=['/usr/bin/timeout','--signal=TERM','--kill-after=5','720',sys.executable,str(base/'run_real.py')]
start=time.monotonic();state=dict(status='running',parentLimitSeconds=720,argv=argv)
with (base/'parent_stdout.txt').open('wb') as log:
    p=subprocess.Popen(argv,stdout=log,stderr=subprocess.STDOUT)
    state['pid']=p.pid;receipt.write_text(json.dumps(state,indent=2)+'\n')
    code=p.wait()
state.update(status='terminal',exitCode=code,elapsedSeconds=time.monotonic()-start)
receipt.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state))
raise SystemExit(0 if code==0 else 1)
