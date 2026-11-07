"""独立记录有限任务的父时间与终态；子过程监督负责清点自有进程组。"""
from pathlib import Path
import argparse
import json
import subprocess
import time

p = argparse.ArgumentParser()
p.add_argument('--seconds', type=int, required=True)
p.add_argument('--record', required=True)
p.add_argument('argv', nargs=argparse.REMAINDER)
a = p.parse_args()
record = Path(a.record); assert not record.exists()
argv = a.argv[1:] if a.argv[:1] == ['--'] else a.argv
assert argv and 0 < a.seconds <= 720
command = ['/usr/bin/timeout', '--signal=TERM', '--kill-after=5', str(a.seconds)] + argv
start = time.monotonic()
state = dict(status='running', parentLimitSeconds=a.seconds, argv=command)
with record.with_suffix('.stdout.txt').open('wb') as log:
    process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
    state['pid'] = process.pid
    record.write_text(json.dumps(state, indent=2)+'\n')
    code = process.wait()
state.update(status='terminal', exitCode=code, elapsedSeconds=time.monotonic()-start)
record.write_text(json.dumps(state, indent=2)+'\n')
print(json.dumps(state))
raise SystemExit(code)
