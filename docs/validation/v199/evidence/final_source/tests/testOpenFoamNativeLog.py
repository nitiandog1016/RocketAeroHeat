"""活动日志跨文件系统边界检查；仅合成Python输出，不调用CFD内核。"""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runOpenFoamProcess import run

out=ROOT/'outputs/v199_validation/native_log01';out.mkdir(exist_ok=False)
case=out/'case';case.mkdir();reads=[]
source='import sys,time\nfor i in range(100):\n sys.stdout.write((str(i)+" observation row\\n")*1000);sys.stdout.flush();time.sleep(.005)\nprint("SYNTHETIC-END",flush=True)\n'
worker=case/'writer.py';worker.write_text(source)
def guard(log,final):
    assert not str(log).startswith('/mnt/')
    text=log.read_text();n=text.count('observation row');assert not reads or n>=reads[-1];reads.append(n)
    if final:assert n==100000 and 'SYNTHETIC-END' in text
    return None
r=run(['python3',worker],case,out/'native',10,[__file__,worker],16,2048,guard)
assert r['passed'] and len(reads)>=2 and r['logExportPassed'] and not r['observationIORetries']
assert (out/'native/stdout.txt').read_text().count('observation row')==100000
report=dict(passed=True,reads=len(reads),lines=100000,process=r,CFDRuns=0,engineBuilds=0)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(passed=True,reads=len(reads),filesystem=r['liveLogFilesystem'],seconds=r['elapsedSeconds'])))
