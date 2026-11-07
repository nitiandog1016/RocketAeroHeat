"""正常退出但无收益的真实结果必须拒绝再次执行和生产外缘导出。"""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from hybridNKRocket import execute,BASE,CASE
from hybridRocketEuler import validate_terminal
from auditHybridNKControl import normal_terminal
from prepareHybridRocketEulerNative import save

start=time.monotonic();run=json.loads((BASE/'rocket_run01/report.json').read_text())
audit=json.loads((BASE/'rocket_audit01/report.json').read_text())
normal_terminal(run['process'],(BASE/'rocket_run01/native/stdout.txt').read_text(),120)
assert audit['passed'] and not audit['candidateScreeningPassed'] and not audit['physicalQualification']
assert audit['newSolvedHeatAreaM2']==0 and audit['qIdentity']=='unknown' and audit['solidTwIdentity']=='unknown'
reasons={}
for name,call in [('closed-execution',execute),('production-export',lambda:validate_terminal(CASE,
    BASE/'rocket_preparation01/report.json',BASE/'rocket_run01/report.json'))]:
    try:call()
    except ValueError as error:reasons[name]=str(error)
    else:raise AssertionError('实际未合格候选被放行: '+name)
out=BASE/'closure01';out.mkdir(exist_ok=False)
save(out/'report.json',dict(passed=True,normalExitDoesNotGrantQualification=True,rejections=reasons,
    CFDRuns=0,elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(passed=True,rejected=len(reasons),elapsedSeconds=time.monotonic()-start)))
