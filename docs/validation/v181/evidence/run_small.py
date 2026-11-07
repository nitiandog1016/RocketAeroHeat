"""新实数热壁小例一次准备/运行/只读对照；原复数对照绝不重算。"""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'scripts'))
from prepareEilmerDevelopmentCase import prepare
from runPreparedEilmerCase import execute
from readEilmerNative import read_field,export_loads,audit_numerics

base=ROOT/'outputs/v181_validation';source=ROOT/'outputs/v178_validation/native03';case=base/'small01'
start=time.monotonic()
condition=json.loads((base/'profile01/condition.json').read_text())
prepared=prepare(source,source.parent/'native03_prepare/report.json',source.parent/'native03_geometry.json',
                 case,base/'small01_prepare',condition,60,256)
assert prepared['passed'],prepared.get('error')
run=execute(case,'run',base/'small01_run',30,256,base/'small01_prepare/report.json')
assert run['passed'],run.get('error')
assert run['commands'][0]['argv'][0]=='lmr-run' and run['commands'][0]['nativeExitCode']==0
nk=json.loads((case/'lmrsim/nk-config').read_text());audit_numerics(nk,condition)
history=np.genfromtxt(case/'lmrsim/diagnostics/nk-diagnostics',names=True,deletechars='')
assert len(history)==8 and all(np.isfinite(history[k]).all() for k in history.dtype.names)
ref=ROOT/'outputs/v180_validation/small01'
current=read_field(case,2,4);old=read_field(ref,2,4)
errors={}
for k in ['p','T','rho','vel.x','vel.y','vel.z','tke','omega']:
    errors[k]=float(np.max(np.abs(current[k]-old[k]))/max(1.,np.max(np.abs(old[k]))))
checks=base/'small01_check';checks.mkdir(exist_ok=False)
loads=case/'lmrsim/loads/0003/blk-0000-bndry-0-wall.dat'
wall=export_loads(case,loads,checks/'wall.csv',base/'small01_run/report.json')
actual=np.genfromtxt(loads,names=True,deletechars='');prior=np.genfromtxt(ref/'lmrsim/loads/0003/blk-0000-bndry-0-wall.dat',names=True,deletechars='')
diff={}
for k in ['p','q_total']:
    diff[k]=float(np.sum(actual['area']*np.abs(actual[k]-prior[k]))/np.sum(actual['area']*np.abs(prior[k])))
passed=max(errors.values())<=1e-3 and diff['p']<=1e-3 and diff['q_total']<=1e-2
r=dict(passed=passed,fieldMaxScaledDifferences=errors,wallAreaL1Differences=diff,wallExport=wall,
       phase=history['phase'].tolist(),lastGlobalRelativeResidual=float(history['global-residual-rel'][-1]),
       diagnosticColumns=list(history.dtype.names),peakGroupRSSBytes=run['commands'][0]['peakObservedGroupRSSBytes'],
       realNativeRuns=1,complexControlReruns=0,physicalQualification=False,newQualifiedAreaM2=0,
       elapsedSeconds=time.monotonic()-start)
(checks/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
assert passed,'真实与复数小例数值相容门槛未过'
