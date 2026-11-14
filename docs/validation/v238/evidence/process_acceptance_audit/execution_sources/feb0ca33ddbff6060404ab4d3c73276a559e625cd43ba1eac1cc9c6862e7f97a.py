"""模块职责：独立核对真实候选/实际接受、具名失败与原生终点，不重新推进或补造未保存状态。"""
from pathlib import Path
"""模块职责：独立核对受限候选与实际接受状态、失败记录及原生终点保存。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runHLLCRocketV238 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC;CLOSED=False


def rows(p):
    with p.open() as f:return list(csv.DictReader(f))


def main():
    out=BASE/'acceptance_audit01';out.mkdir(exist_ok=False);air=PhaseAwareViscousAir();geometry=np.load(str(run.GEOMETRY))
    process=read(BASE/'process_rocket/report.json');assert process['status']=='terminal' and process['inputsStable'] and not process['liveGroupAfterCleanup']
    steps=rows(CASE/'steps.csv');last=len(steps);events=rows(CASE/'accepted_steps.csv');actuals=rows(CASE/'actual_accepted_steps.csv');failures=rows(CASE/'failed_step.csv')
    def valid(u):
        if not np.isfinite(u).all() or u[0]<=0:return False
        try:
            t=inverse_energy(air,np.array([u[4]/u[0]-.5*np.sum((u[1:4]/u[0])**2)-300000]));air.evaluate(u[0]*air.R*t,t);return True
        except ValueError:return False
    def vector(r,prefix):return np.array([float(r[prefix+str(j)]) for j in range(5)])
    records={};alphas=[];directionError=0
    for r in events:
        key=(int(r['update']),int(r['localID']),int(r['globalID']));assert key not in records and key[0]<=last+1 and 0<=key[2]<335726
        old=vector(r,'old');delta=vector(r,'delta');accepted=vector(r,'accepted');upper=float(r['nativeUpper']);alpha=upper;count=0
        assert valid(old) and 0<upper<=1
        while not valid(old+alpha*delta):
            alpha*=.5;count+=1;assert alpha>=1e-10 and count<=34
        assert alpha==float(r['acceptedAlpha']) and count==int(r['reductions']) and count>0
        error=float(np.max(abs(accepted-(old+alpha*delta))/np.maximum(1,abs(accepted))));directionError=max(directionError,error);assert error<2e-12
        assert valid(accepted);records[key]=accepted;alphas.append(alpha)
    actualKeys=set();actualComplete=0
    for r in actuals:
        key=(int(r['update']),int(r['localID']),int(r['globalID']));assert key in records and key not in actualKeys
        candidate=vector(r,'candidate');actual=vector(r,'actual')
        assert np.array_equal(candidate.view(np.uint64),records[key].view(np.uint64))
        assert np.array_equal(candidate.view(np.uint64),actual.view(np.uint64)),('ACTUAL_STATE',key)
        actualKeys.add(key);actualComplete+=key[0]<=last
    assert {k for k in records if k[0]<=last}<=actualKeys
    failureReports=[]
    for r in failures:
        point=int(r['globalID']);assert int(r['update'])==last+1 and 0<=point<335726
        xyz=np.array([float(r[d]) for d in ('x','y','z')]);assert np.array_equal(xyz,geometry['points'][point])
        old=vector(r,'old');delta=vector(r,'delta');upper=float(r['nativeUpper']);assert valid(old) and np.isfinite(delta).all()
        assert r['reason']=='PHASE_STEP_NO_ADMISSIBLE_PROGRESS'
        alpha=upper;found=False;roundingNoProgress=False;tested=0
        if alpha>=1e-10:
            for reduction in range(35):
                candidate=old+alpha*delta;tested+=1
                if valid(candidate):
                    roundingNoProgress=np.array_equal(candidate.view(np.uint64),old.view(np.uint64));found=not roundingNoProgress;break
                alpha*=.5
                if alpha<1e-10:break
        assert not found,('FAILURE_NOT_INDEPENDENTLY_CONFIRMED',point,alpha)
        oldT=float(inverse_energy(air,np.array([old[4]/old[0]-.5*np.sum((old[1:4]/old[0])**2)-300000])))
        failureReports.append(dict(globalID=point,update=int(r['update']),coordinatesM=xyz.tolist(),oldTemperatureK=oldT,
            nativeUpper=upper,independentNoAdmissibleProgress=True,roundingNoProgress=roundingNoProgress,testedCandidates=tested))
    completion=read(CASE/'native_completion.json') if (CASE/'native_completion.json').exists() else None
    finalName=CASE/('state_%06d.bin'%last);nativeFiles=list(CASE.glob('restart_flow*.csv'));nativeBitEqual=None
    if completion:
        summary=read(CASE/'acceptance_summary.json');assert summary['checked']==last*335726 and summary['limited']==len(events) and summary['forcedFinalNativeOutput']
        assert len(actualKeys)==len(records) and not failures and len(nativeFiles)==1 and finalName.exists()
        frame=read_frame(finalName,335726,False)['values'];frame=frame[np.argsort(frame[:,0])]
        raw=np.loadtxt(str(nativeFiles[0]),delimiter=',',skiprows=1,usecols=range(10));raw=raw[np.argsort(raw[:,0])]
        assert np.array_equal(raw[:,0],np.arange(335726))
        nativeBitEqual=np.array_equal(raw[:,1:4].copy().view(np.uint64),frame[:,1:4].copy().view(np.uint64)) and np.array_equal(raw[:,4:10].copy().view(np.uint64),frame[:,5:11].copy().view(np.uint64))
        assert nativeBitEqual,'FINAL_NATIVE_BINARY_BITS'
    rejection=rows(CASE/'rejections.csv')
    report=dict(passed=True,realCompletedUpdates=last,limitedCandidates=len(records),actualAcceptedEvents=len(actualKeys),
        actualAcceptedInCompleteFlowSAUpdates=actualComplete,partialFlowOnlyAcceptedEvents=len(actualKeys)-actualComplete,
        candidatesNotApplied=len(records)-len(actualKeys),actualAcceptedFiveStateComponentsBitEqual=5*len(actualKeys),
        independentActualDirectionMaxError=directionError,acceptedAlphaRange=[min(alphas),max(alphas)] if alphas else None,
        nativeGuardFailures=failureReports,oldStateRecoveryEvents=len(rejection),productionGuardFailureObserverTested=bool(failures),
        completeFinalBinaryPresent=finalName.exists(),nativeFinalCSVWritten=bool(nativeFiles),nativeFinalSixStateCoordinatesBitEqual=nativeBitEqual,
        nativeFinalSHA256=sha(nativeFiles[0]) if len(nativeFiles)==1 else None,nativeCompletion=completion,
        originalFaces=19752,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'report.json',report);print(json.dumps(report))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v238 rocket acceptance audit closed')
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runHLLCRocketV238.py',ROOT/'scripts/runHLLCCompatibilityV238.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'ROCKET_PLAN.md',
          ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json',
          run.GEOMETRY,BASE/'process_rocket/report.json']+[p for p in CASE.iterdir() if p.is_file()]
        raise SystemExit(0 if run.run.guard.supervise('acceptance_audit',[sys.executable,str(Path(__file__).resolve()),'task'],90,1024,locks,BASE,'audit') else 1)
