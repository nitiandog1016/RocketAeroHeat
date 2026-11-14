"""模块职责：独立核对生产入口真实线性增量、守卫轨迹和原生终点保存，不重跑CFD。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runPhaseBranchV237 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC;CLOSED=False


def rows(p):
    with p.open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit():
    out=BASE/'production_audit01';out.mkdir(exist_ok=False);air=PhaseAwareViscousAir()
    process=read(BASE/'process_rocket/report.json');assert process['inputsStable'] and not process['liveGroupAfterCleanup']
    events=rows(CASE/'accepted_steps.csv');steps=rows(CASE/'steps.csv');last=len(steps);ids=set();errors={};byStep={};conditions=[]
    geometry=np.load(str(run.GEOMETRY))
    config={k.strip():v.strip() for k,v in (line.split('=',1) for line in (CASE/'case.cfg').read_text(encoding='utf8').splitlines() if '=' in line)}
    free=air.evaluate(float(config['FREESTREAM_PRESSURE']),float(config['FREESTREAM_TEMPERATURE']))
    def valid(u):
        if not np.isfinite(u).all() or u[0]<=0:return False
        try:
            e=u[4]/u[0]-.5*np.sum((u[1:4]/u[0])**2)-300000
            t=inverse_energy(air,np.array([e]));air.evaluate(u[0]*air.R*t,t);return True
        except ValueError:return False
    def check(name,a,b,tol=2e-12):
        a=np.asarray(a);b=np.asarray(b);value=float(np.max(abs(a-b)/np.maximum(1,abs(b))))
        errors[name]=max(errors.get(name,0),value);assert value<=tol,(name,value)
    alphas=[];reductions=[]
    for r in events:
        step=int(r['update']);point=int(r['globalID']);assert 0<=point<335726 and step<=last+1
        old=np.array([float(r['old%d'%j]) for j in range(5)]);d=np.array([float(r['delta%d'%j]) for j in range(5)]);new=np.array([float(r['accepted%d'%j]) for j in range(5)])
        upper=float(r['nativeUpper']);assert valid(old) and 0<upper<=1
        alpha=upper;count=0
        while not valid(old+alpha*d):
            alpha*=.5;count+=1;assert alpha>=1e-10 and count<=34
        assert alpha==float(r['acceptedAlpha']) and count==int(r['reductions']) and count>0
        check('sameActualFiveStateDirection',new,old+alpha*d)
        assert valid(new);alphas.append(alpha);reductions.append(count);ids.add(point);byStep.setdefault(step,[]).append((point,old,d,new))
        speed2=np.sum((new[1:4]/new[0])**2);t=inverse_energy(air,np.array([new[4]/new[0]-.5*speed2-300000]));q=air.evaluate(new[0]*air.R*t,t)
        conditions.append(dict(update=step,globalID=point,coordinatesM=geometry['points'][point].tolist(),
            acceptedTemperatureK=float(t[0]),acceptedPressurePa=float(q['pressurePa'][0]),
            entropyRelativeFreestreamJkgK=float(q['entropyJkgK'][0]-free['entropyJkgK']),
            physicalTotalEnthalpyJkg=float(q['enthalpyJkg'][0]+.5*speed2),completedUpdate=step<=last))
    # 已保存完整步按全局ID逐位比对实际接受量；未完成步只称候选，不授接受身份。
    savedChecks=0;initialEqual=None;completion=None;nativeCSV=None;nativeBits=None;unchangedFrames=[]
    frames=sorted(CASE.glob('state_*.bin'))
    for p in frames:
        step=int(p.stem.split('_')[-1])
        if step not in byStep and step not in (0,1,last):continue
        frame=read_frame(p,335726,step==0);a=frame['values'];a=a[np.argsort(a[:,0])]
        for point,old,d,new in byStep.get(step,[]):
            assert np.array_equal(a[point,5:10].view(np.uint64),new.view(np.uint64)),('SAVED_ACCEPTANCE',step,point)
            savedChecks+=5
        if step==0:
            old=read_frame(ROOT/'outputs/v236_validation/rocket01/state_000000.bin',335726,True)['values'];old=old[np.argsort(old[:,0])]
            initialEqual=np.array_equal(a[:,1:11].copy().view(np.uint64),old[:,1:11].copy().view(np.uint64));assert initialEqual
        if step and (not byStep or step<min(byStep)):
            prior=ROOT/'outputs/v236_validation/rocket01'/p.name
            if prior.exists():
                old=read_frame(prior,335726,False)['values'];old=old[np.argsort(old[:,0])]
                assert np.array_equal(a.view(np.uint64),old.view(np.uint64)),('PRE_GUARD_TRAJECTORY',step)
                unchangedFrames.append(step)
        if step==last and (CASE/'native_completion.json').exists():
            completion=read(CASE/'native_completion.json');assert completion['updates']==last
            files=list(CASE.glob('restart_flow*.csv'));assert len(files)==1
            nativeCSV=files[0];native=np.loadtxt(str(nativeCSV),delimiter=',',skiprows=1,usecols=range(10));native=native[np.argsort(native[:,0])]
            assert np.array_equal(native[:,0],np.arange(335726))
            nativeBits=np.array_equal(native[:,1:4].copy().view(np.uint64),a[:,1:4].copy().view(np.uint64)) and np.array_equal(native[:,4:10].copy().view(np.uint64),a[:,5:11].copy().view(np.uint64))
            assert nativeBits,'NATIVE_FINAL_SIX_STATE_BITS'
    summary=read(CASE/'acceptance_summary.json') if (CASE/'acceptance_summary.json').exists() else None
    if summary:
        assert summary['limited']==len(events) and summary['checked']==last*335726 and summary['forcedFinalNativeOutput']
    rejected=rows(CASE/'rejections.csv');assert len(rejected)==sum(int(r['lastFlowRejectionCount']) for r in steps)
    report=dict(passed=True,productionCallsiteExecuted=bool(events),realCompletedUpdates=last,limitedCandidates=len(events),
        limitedCandidatesInCompleteUpdates=sum(len(v) for s,v in byStep.items() if s<=last),
        limitedCandidatesInIncompleteUpdate=sum(len(v) for s,v in byStep.items() if s>last),uniqueLimitedNodes=len(ids),
        rawLinearIncrementRecorded=True,independentHalvingVerified=True,acceptedAlphaRange=[min(alphas),max(alphas)] if alphas else None,
        reductionsRange=[min(reductions),max(reductions)] if reductions else None,actualSavedStateComponentsBitChecked=savedChecks,
        initialVsV236BitEqual=initialEqual,preGuardEntireFramesBitEqualToV236=unchangedFrames,nativeFinalCSVWritten=bool(nativeCSV),nativeFinalSixStatesCoordinatesBitEqual=nativeBits,
        nativeFinalSHA256=sha(nativeCSV) if nativeCSV else None,nativeSummary=summary,recoveryEvents=len(rejected),
        errorMetrics=errors,nativeCompletion=completion,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'limited_state_conditions.json',dict(events=conditions,
        entropyInterpretation='local state diagnostic only; cooled walls/transport require a complete entropy flux budget; not a proof of entropy inequality failure'))
    save(out/'report.json',report);print(json.dumps(report))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v237 production audit closed')
    if len(sys.argv)>1 and sys.argv[1]=='task':audit()
    else:
        run.run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runPhaseBranchV237.py',ROOT/'scripts/runPhaseAcceptanceV237.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            ROOT/'scripts/singleDomainFrame.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',
            ROOT/'data/thermo/phase_aware_dry_air_v1.json',run.GEOMETRY,BASE/'process_rocket/report.json']+[p for p in CASE.iterdir() if p.is_file()]
        raise SystemExit(0 if run.run.guard.supervise('production_audit',[sys.executable,str(Path(__file__).resolve()),'task'],60,1024,locks,BASE,'audit') else 1)
