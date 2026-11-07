"""120秒监督新增MATLAB发展准备/失败读取入口，不启动CFD。"""
from pathlib import Path
import subprocess,time,json,sys
out=Path('outputs/v180_validation/matlab_supervision.json');assert not out.exists()
start=time.monotonic()
cmd=[sys.executable,'scripts/runMatlabResearch.py','eilmer_development01',
     "addpath(pwd); addpath('tests'); testEilmerDevelopmentEntry();",'--folder','outputs/v180_validation']
p=subprocess.Popen(cmd,creationflags=subprocess.CREATE_NO_WINDOW);stop=None
try:code=p.wait(timeout=120)
except subprocess.TimeoutExpired:
    stop='120-second-budget';subprocess.run(['taskkill','/PID',str(p.pid),'/T','/F'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT);code=p.wait()
out.write_text(json.dumps(dict(exitCode=code,stopReason=stop,elapsedSeconds=time.monotonic()-start,limitSeconds=120),indent=2)+'\n')
print(out.read_text());raise SystemExit(code)
