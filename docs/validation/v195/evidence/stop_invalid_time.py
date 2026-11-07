"""只停止本次已登记且仍存活的原生进程组；首步巨大伪时间超出设计。"""
from pathlib import Path
import json
import os
import signal
import time
base=Path('/mnt/e/RocketAeroHeat/outputs/v195_validation/run01')
r=json.loads((base/'native/report.json').read_text());pid=r['processId']
expected=r['argv'];assert r['cwd']=='/mnt/e/RocketAeroHeat/outputs/v195_validation/rocket01'
proc=Path('/proc')/str(pid);cmd=(proc/'cmdline').read_bytes().decode().strip('\0').split('\0')
stat=(proc/'stat').read_text();pgid=os.getpgid(pid)
assert cmd==expected and pgid==pid and proc.exists()
log=(base/'native/stdout.txt').read_text()
assert 'Step=      1 t= 3.000e+98 dt= 3.000e+98' in log
report=dict(reason='manual scientific stop: observed first pseudo-time dt=3e98 incompatible with configured dt_max=.001/max_time=1',
            pid=pid,pgid=pgid,argv=cmd,procStatBefore=stat,observationUTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
            notUserCancellation=True,noRestartAuthorized=True)
path=base/'invalid_time_stop.json';assert not path.exists();path.write_text(json.dumps(report,indent=2)+'\n')
os.killpg(pid,signal.SIGTERM)
print(json.dumps(dict(signaledPID=pid,reason=report['reason'])))
