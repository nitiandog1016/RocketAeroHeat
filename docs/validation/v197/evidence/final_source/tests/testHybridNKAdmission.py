"""真实NK合同的故障拒绝；不启动CFD、不修改绑定的源与场。"""
from pathlib import Path
import copy
import json
import sys
import tempfile
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from hybridNKRocket import check_admission,BASE
from auditHybridNKControl import normal_terminal,read_diagnostics
from prepareHybridRocketEulerNative import save

start=time.monotonic();a=json.loads((BASE/'rocket_admission.json').read_text())
check_admission(a);rejections=[]
def reject(label,callback):
    try:callback()
    except (ValueError,KeyError):rejections.append(label)
    else:raise AssertionError('应拒绝但被接受: '+label)
for key,value in [('schema','wrong'),('profile','wrong'),('temperatureRangeK',[150.,6000.]),
                  ('cells',64981),('wallFaces',19751),('physicalQualification',True),
                  ('savedSteps',[0,120]),('case',str(BASE/'tetra01')),('record',str(BASE/'run01')),
                  ('runtimeSHA256','0'*64),('evidenceSHA256',{})]:
    b=copy.deepcopy(a);b[key]=value;reject(key,lambda:check_admission(b))
b=copy.deepcopy(a);b['budget']['maxSteps']=240;reject('extended-budget',lambda:check_admission(b))
r=json.loads((BASE/'run01/report.json').read_text())['process']
t=(BASE/'run01/native/stdout.txt').read_text();normal_terminal(r,t,120)
for key,value in [('nativeExitCode',1),('inputsStable',False),('stopReason','temperature'),('liveGroupAfterCleanup',{'1':{}})]:
    b=copy.deepcopy(r);b[key]=value;reject(key,lambda:normal_terminal(b,t,120))
reject('real-instead-of-complex',lambda:normal_terminal(r,t.replace('Number-type: complex','Number-type: real'),120))
reject('wrong-final-step',lambda:normal_terminal(r,t,240))
path=BASE/'tetra01/lmrsim/diagnostics/nk-diagnostics';read_diagnostics(path,120)
with tempfile.TemporaryDirectory() as directory:
    p=Path(directory)/'history';lines=path.read_text().splitlines()
    p.write_text('\n'.join(lines[:-1]));reject('missing-diagnostic-row',lambda:read_diagnostics(p,120))
out=BASE/'contract01';out.mkdir(exist_ok=False)
save(out/'report.json',dict(passed=True,rejections=rejections,rejectionCount=len(rejections),
    realControlTerminalPassed=True,actualAdmissionPassed=True,CFDRuns=0,elapsedSeconds=time.monotonic()-start))
print(json.dumps(dict(passed=True,rejectionCount=len(rejections),elapsedSeconds=time.monotonic()-start)))
