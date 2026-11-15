"""模块职责：独立解析原生调度日程、投影与六量保持；此控制没有流动更新。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runMUSCLRampV239 as run
from runSU2RocketPhaseV222 import ROOT,read,save
BASE=run.BASE;DOC=run.DOC;CLOSED=True


def main():
    process=read(BASE/'process_control/report.json');assert process['passed'] and process['inputsStable'] and not process['liveGroupAfterCleanup']
    out=BASE/'control_audit01';out.mkdir(exist_ok=False)
    with (run.CONTROL/'ramp.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==206;errors=[]
    for i,r in enumerate(rows):
        step=i//2+1;before=i%2==0;assert int(r['update'])==step and r['stage']==('before_iterate' if before else 'after_monitor')
        expected=.5*(1-np.cos(np.pi*min(100,max(0,step-(2 if before else 1)))/100))
        assert float(r['newtonRelaxation'])==1
        errors += [abs(float(r['ramp'])-expected),abs(float(r['nativeSyntheticProjection'])+.3*expected)]
    assert max(errors)<2e-13
    a=np.loadtxt(str(run.CONTROL/'state_pairs.csv'),delimiter=',',skiprows=1)
    assert a.shape==(144,4) and np.array_equal(a[:,2].copy().view(np.uint64),a[:,3].copy().view(np.uint64))
    control=read(run.CONTROL/'control.json');assert control['flowUpdates']==control['SAUpdates']==control['phaseGuardCalls']==0
    save(out/'report.json',dict(passed=True,schedulerCalls=103,flowUpdates=0,SAUpdates=0,phaseGuardCalls=0,
      sixStateComponentsBitEqual=144,initialRampZero=True,update1and2RampZero=True,update52RampHalf=True,fullOrderFirstUpdate=102,
      monitorAppliesToNextIterate=True,nativeSyntheticProjectionChecked=True,wholeResidualReconstructionRecomputed=False,
      maxIndependentError=max(errors),physicalQualification=False,newQualifiedAreaM2=0))
    print(json.dumps(dict(passed=True,fullOrderFirstUpdate=102,maxIndependentError=max(errors),flowUpdates=0)))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v239 schedule audit closed')
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runMUSCLRampV239.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',BASE/'process_control/report.json']+list(run.CONTROL.iterdir())
        raise SystemExit(0 if run.guard.supervise('control_audit',[sys.executable,str(Path(__file__).resolve()),'task'],30,1024,locks,BASE,'audit') else 1)
