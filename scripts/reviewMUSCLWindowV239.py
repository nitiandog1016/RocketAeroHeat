"""模块职责：只用已经完成的完整二阶窗口判断既定停止谓词；不伪造300步原生决定。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runMUSCLRampV239 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC


def rows(path):
    with path.open() as f:return list(csv.DictReader(f))


def main():
    out=BASE/'window01';out.mkdir(exist_ok=False)
    audit=read(BASE/'real_audit01/report.json');assert audit['passed'] and audit['successfulUpdates']==258
    contract=read(BASE/'rocket_contract.json');assert contract['checkUpdates']==[300,400]
    stage={(int(r['update']),int(r['system'])):r for r in rows(CASE/'stages.csv') if r['stage']=='after_time'}
    logs=np.array([[np.log10(max(1e-300,float(stage[(step,0 if j<5 else 2)]['rms%d'%j]))) for j in range(6)] for step in range(151,251)])
    ramp={(int(r['update']),r['stage']):float(r['ramp']) for r in rows(CASE/'ramp.csv')}
    assert all(ramp[(step,'flow_after_sources')]==1 for step in range(151,251))
    change={}
    for step in (150,200,250):
        a=np.loadtxt(str(BASE/'real_audit01'/('original_faces_%06d.csv'%step)),delimiter=',',skiprows=1)
        if step>150:change[step]=np.array([a[:,1]@abs(a[:,j]-previous[:,j])/max(previous[:,1]@abs(previous[:,j]),1) for j in (2,3)])
        previous=a
    ratio=10**(logs[50:].mean(axis=0)-logs[:50].mean(axis=0));gm=float(np.exp(np.log(ratio).mean()))
    old=change[200];now=change[250];threshold=np.array([contract['pressureChangeScreen'],contract['heatChangeScreen']])
    factor=contract['significantWorseningRatio'];assert factor==1.25
    pqlimit=factor*np.maximum(old,threshold)
    notWorse=bool(np.all(ratio<=factor) and np.all(now<=pqlimit))
    good=notWorse and bool(gm<=contract['geometricMeanImprovementRatio'] or np.all(logs[50:].mean(axis=0)<=-8)) and bool(np.all((now<=threshold)|(now<=.95*old)))
    assert not rows(CASE/'benefit_decisions.csv')
    report=dict(passed=True,sourceAuditSHA256=sha(BASE/'real_audit01/report.json'),native300CheckExecuted=False,
      recordedCompletedUpdates=258,completedWindowEnd=250,olderUpdates=[151,200],newerUpdates=[201,250],
      fullOrderBothWindows=True,residualMeanRatios=ratio.tolist(),geometricMeanRatio=gm,
      oldPressureHeatChanges=old.tolist(),newPressureHeatChanges=now.tolist(),pressureHeatWorseningLimits=pqlimit.tolist(),
      pressureHeatChangeRatios=(now/old).tolist(),noSignificantDeterioration=notWorse,benefit=good,
      firstCompletedWindowAlreadyTriggersWorse=not notWorse,
      nativeStopPredicateUsesAnyWorseWindow=True,remainingWindowCannotClearExistingWorse=not notWorse,
      actualOriginalTerminationReasonUnknown=True,claimedNativeBenefitStop=False,
      furtherCFDNeededToEstablishThisPredicate=False,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'report.json',report);print(json.dumps(report))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runMUSCLRampV239.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
          BASE/'real_audit01/report.json',BASE/'rocket_contract.json',run.STAGE/'commonBranchV239.hpp',CASE/'stages.csv',CASE/'ramp.csv',CASE/'benefit_decisions.csv']
        locks += [BASE/'real_audit01'/('original_faces_%06d.csv'%s) for s in (150,200,250)]
        raise SystemExit(0 if run.guard.supervise('window',[sys.executable,str(Path(__file__).resolve()),'task'],30,1024,locks,BASE,'audit') else 1)
