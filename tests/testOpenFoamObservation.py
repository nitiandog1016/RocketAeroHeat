"""监测IO注入：短暂读取失败可恢复，持续或非瞬时错误不被隐藏。"""
from pathlib import Path
import sys,errno,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runOpenFoamProcess import observe_io

start=time.monotonic();checks=[]
for code in [errno.ENODATA,errno.EAGAIN,errno.EINTR]:
    calls=[];events=[]
    def temporary():
        calls.append(1)
        if len(calls)<3:raise OSError(code,'injected')
        return 'same observation, no process launch'
    assert observe_io(temporary,'fixture',events).startswith('same') and len(calls)==3 and len(events)==2
    checks.append('temporary-'+str(code))
    calls=[];events=[]
    def persistent():calls.append(1);raise OSError(code,'injected')
    try:observe_io(persistent,'fixture',events)
    except OSError:assert len(calls)==3 and len(events)==3
    else:raise AssertionError('persistent failure hidden')
    checks.append('persistent-'+str(code))
events=[]
def other():raise OSError(errno.EACCES,'injected')
try:observe_io(other,'fixture',events)
except OSError:assert not events
else:raise AssertionError('unrelated error hidden')
checks.append('unrelated-error')
r=dict(passed=True,checks=checks,elapsedSeconds=time.monotonic()-start)
(ROOT/'outputs/v199_validation/observation_tests.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
