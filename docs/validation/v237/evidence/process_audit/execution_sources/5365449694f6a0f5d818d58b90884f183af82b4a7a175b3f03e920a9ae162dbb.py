"""模块职责：用独立分段Cv反算核对真实事件的更新接受、原生primitive及能量坐标；不推进流场。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runPhaseAcceptanceV237 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy
BASE=run.BASE;DOC=run.DOC;CLOSED=False


def rows(p):
    with p.open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);air=PhaseAwareViscousAir();errors={}
    controls=rows(BASE/'control01/accepted.csv');events=rows(ROOT/'outputs/v236_validation/rocket01/rejections.csv')
    assert len(controls)==len(events)==466 and read(BASE/'control01/report.json')['nativeCNSAcceptedWithoutFallback']
    def physical(u,offset=300000.):
        if not np.isfinite(u).all() or u[0]<=0:raise ValueError('density or finite')
        e=u[4]/u[0]-.5*np.sum((u[1:4]/u[0])**2)-offset
        t=inverse_energy(air,np.array([e]));return air.evaluate(u[0]*air.R*t,t)
    def valid(u,offset=300000.):
        try:return physical(u,offset)
        except ValueError:return None
    def check(name,a,b,tol=2e-9):
        a=np.asarray(a);b=np.asarray(b);error=float(np.max(abs(a-b)/np.maximum(abs(b),1)))
        errors[name]=max(errors.get(name,0),error);assert error<tol,(name,error)
    reductionCounts={};temperatures=[];alphas=[]
    for number,(r,e) in enumerate(zip(controls,events),1):
        assert int(r['event'])==number and int(r['globalID'])==int(e['globalID']) and r['update']==e['update']
        old=np.array([float(r['old%d'%j]) for j in range(5)]);delta=np.array([float(r['delta%d'%j]) for j in range(5)]);accepted=np.array([float(r['accepted%d'%j]) for j in range(5)])
        sourceOld=np.array([float(e['old%d'%j]) for j in range(5)]);trial=np.array([float(e['trial%d'%j]) for j in range(5)])
        assert np.array_equal(old.view(np.uint64),sourceOld.view(np.uint64)) and np.array_equal(delta.view(np.uint64),(trial-old).view(np.uint64))
        assert valid(old) is not None and valid(trial) is None
        expectedAlpha=1.;reductions=0
        while valid(old+expectedAlpha*delta) is None:
            expectedAlpha*=.5;reductions+=1;assert expectedAlpha>=1e-10 and reductions<=34
        assert expectedAlpha==float(r['alpha'])==float(r['shiftedAlpha']) and reductions==int(r['reductions'])
        check('sameFiveStateDirection',accepted,old+expectedAlpha*delta)
        q=physical(accepted);check('nativePressure',float(r['p']),q['pressurePa']);check('nativeTemperature',float(r['T']),q['temperatureK'])
        physicalOld=old.copy();physicalOld[4]-=300000*old[0];physicalDelta=delta.copy();physicalDelta[4]-=300000*delta[0]
        ref=valid(physicalOld+expectedAlpha*physicalDelta,0.);assert ref is not None
        check('referencePressure',ref['pressurePa'],q['pressurePa']);check('referenceTemperature',ref['temperatureK'],q['temperatureK'])
        assert float(r['legalAlpha'])==1 and float(r['zeroAlpha'])==1 and float(r['shiftedMaxError'])<2e-12
        reductionCounts[reductions]=reductionCounts.get(reductions,0)+1;temperatures.append(float(r['T']));alphas.append(expectedAlpha)
    rejected=rows(BASE/'control01/rejected.csv');assert len(rejected)==7
    expected={'invalid_old':'PHASE_STEP_ILLEGAL_OLD','nonfinite_delta':'PHASE_STEP_NONFINITE_INCREMENT',
      'negative_upper':'PHASE_STEP_UPPER_OR_REFERENCE','above_one_upper':'PHASE_STEP_UPPER_OR_REFERENCE','nan_upper':'PHASE_STEP_UPPER_OR_REFERENCE',
      'zero_progress':'PHASE_STEP_NO_ADMISSIBLE_PROGRESS','no_admissible_progress':'PHASE_STEP_NO_ADMISSIBLE_PROGRESS'}
    assert {r['case']:r['error'] for r in rejected}==expected
    save(out/'report.json',dict(passed=True,realEvents=466,independentAnalyticCvInverse=True,nativeCNSAcceptedWithoutFallback=True,
        acceptedAlphaRange=[min(alphas),max(alphas)],reductionsHistogram=reductionCounts,acceptedTemperatureRangeK=[min(temperatures),max(temperatures)],
        errorMetrics=errors,legalAndZeroUpdatesPreserved=932,structuralRejections=len(rejected),
        energyReferenceGuardEquivalentForSameUpper=True,legacyRelativeRelaxationReferenceInvariant=False,
        replayIncrement='recorded trial minus old; raw linear solve increment and original relaxation not saved in v236',
        productionCallsiteExecuted=False,flowUpdates=0,newQualifiedAreaM2=0,physicalQualification=False))
    print(json.dumps(dict(passed=True,events=466,alpha=[min(alphas),max(alphas)],maxError=max(errors.values()))))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v237 acceptance audit closed')
    if len(sys.argv)>1 and sys.argv[1]=='task':audit()
    else:
        run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runPhaseAcceptanceV237.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
          ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/su2PhaseAirFields.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json',
          ROOT/'outputs/v236_validation/rocket01/rejections.csv']+list((BASE/'control01').glob('*.csv'))+[BASE/'control01/report.json']
        raise SystemExit(0 if run.guard.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],60,1024,locks,BASE,'audit') else 1)
