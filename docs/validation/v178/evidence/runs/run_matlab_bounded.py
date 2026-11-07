from pathlib import Path
import subprocess,time,json,sys
base=Path('outputs/v178_validation');pfile=base/'matlab_supervision.json'
assert not pfile.exists()
t=time.monotonic()
cmd=[sys.executable,'scripts/runMatlabResearch.py','eilmer_entry01',"addpath(pwd); addpath('tests'); testEilmerMatlabEntry();",'--folder','outputs/v178_validation']
p=subprocess.Popen(cmd,creationflags=subprocess.CREATE_NO_WINDOW)
stop=None
try:code=p.wait(timeout=120)
except subprocess.TimeoutExpired:
    stop='120-second-budget';subprocess.run(['taskkill','/PID',str(p.pid),'/T','/F'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT);code=p.wait()
r=dict(exitCode=code,stopReason=stop,elapsedSeconds=time.monotonic()-t,limitSeconds=120)
pfile.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));sys.exit(code)
