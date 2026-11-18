"""模块职责：80位独立分段积分/二分参考，核对真实收到的target及成对微扰精度。"""
import csv,json,re,time,bisect
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=ROOT/'outputs/v260_validation';CASE=BASE/'case01';CLOSED=False
TABLE=ROOT/'outputs/v220_validation/source/SU2_CFD/include/fluid/PhaseAirTable.hpp'


def audit():
    started=time.monotonic();out=BASE/'audit01';out.mkdir(exist_ok=False)
    native=read(CASE/'native_completion.json');process=read(BASE/'process_query/report.json')
    assert process['passed'] and process['inputsStable'] and not process['liveGroupAfterCleanup'] and native['normalReturn']
    assert native['primitiveInputsBitwiseMatched'] and native['oldTracedMappingBitwiseOriginal22Fields'] and native['longDoubleDigits']==64
    raw=np.array([[float(x) for x in line.strip('{} ,').split(',')] for line in TABLE.read_text(encoding='utf8').splitlines() if re.match(r'^\{[0-9]',line)])
    assert raw.shape==(3882,6)
    def dec(x):return Decimal.from_float(float(x))
    records=list(csv.DictReader((CASE/'states.csv').open()));controls=list(csv.DictReader((CASE/'controls.csv').open()));hs=list(csv.DictReader((CASE/'hs.csv').open()))
    rows={};air=PhaseAwareViscousAir();errors=[];references={}
    with localcontext() as ctx:
        ctx.prec=80;table=[[dec(x) for x in row] for row in raw];temperatures=[r[0] for r in table];R=dec(287.0491267875601);shift=Decimal(300000)
        def value(t,kind):
            i=min(len(table)-2,bisect.bisect_right(temperatures,t)-1);a,b=table[i:i+2];dt=t-a[0]
            return a[2]+a[1]*dt+(b[1]-a[1])/(b[0]-a[0])*dt*dt/2-(R*t if kind=='E' else 0)+shift
        def reference(target,kind):
            key=(float(target).hex(),kind)
            if key in references:return references[key]
            target=dec(target);lo,hi=temperatures[0],temperatures[-1]
            for k in range(160):
                mid=(lo+hi)/2
                if value(mid,kind)<target:lo=mid
                else:hi=mid
            t=(lo+hi)/2;res=abs(value(t,kind)-target)
            result=(t,res,res<Decimal('1e-20'));references[key]=result;return result
        for row in records:
            sample=int(row['sample']);mode=row['mode'];rows[(sample,mode)]=row
            target=float(row['receivedE']);tref,res,unique=reference(target,'E');assert unique,('ACTUAL_REPLAY_REFERENCE_NOT_UNIQUE',sample,mode,str(res))
            t=float(row['T']);p=float(row['p']);g=air.evaluate(p,t);tr=air.interpolated_transport(p,t)
            checks=[abs(float(g['densityKgm3'])/float(row['rho'])-1),abs((float(g['internalEnergyJkg'])+300000-float(row['e']))/max(1,abs(float(row['e'])))),
                abs((float(g['enthalpyJkg'])+300000-float(row['h']))/max(1,abs(float(row['h'])))),abs(float(tr['dynamicViscosityPas'])/float(row['mu'])-1),abs(float(tr['conductivityWmK'])/float(row['k'])-1)]
            assert max(checks)<2e-10
            err=dec(t)-tref;errors.append(dict(sample=sample,mode=mode,targetE=target,referenceT=str(tref),errorK=float(err),referenceResidualJkg=str(res),physicalCheckMax=max(checks),method=row['method']))
        controlsReport=[]
        for row in controls:
            expect=row['expectedValid']=='1';assert (row['oldAccepted']=='1')==expect
            entry=dict(row);entry['acceptancePassed']=(row['newAccepted']=='1')==expect
            if row['newAccepted']=='1':
                target=float(row['target']);ref,res,unique=reference(target,row['kind']);error=abs(dec(float(row['newT']))-ref)
                entry.update(referenceT=str(ref),referenceResidualJkg=str(res),uniqueReference=unique,errorK=float(error));assert error<=Decimal('2e-10')
                if row['label']=='node':assert float(row['newT'])==float(row['sourceT'])
            controlsReport.append(entry)
        for row in hs:
            if row['mode'].endswith('H'):
                ref,res,unique=reference(float(row['targetH']),'H');assert unique and abs(dec(float(row['T']))-ref)<=Decimal('2e-10')
        identities=read(CASE/'input_identity.json')['samples'];byNode={};errorMap={(r['sample'],r['mode']):r for r in errors}
        for item in identities:byNode.setdefault(item['globalID'],{})[item['sourceFrame']]=item['sample']
        pairs=[]
        for node,stages in byNode.items():
            for frame in ['sigma','half']:
                i,j=stages['initial'],stages[frame];ri=Decimal(errorMap[(i,'old')]['referenceT']);rj=Decimal(errorMap[(j,'old')]['referenceT']);signal=rj-ri
                item=dict(globalID=node,frame=frame,referenceDeltaT=str(signal))
                for mode in ['old','new']:
                    a,b=rows[(i,mode)],rows[(j,mode)];observed=dec(float(b['T']))-dec(float(a['T']));error=observed-signal
                    item[mode+'ErrorK']=float(error);item[mode+'RelativeError']=float(abs(error/signal)) if signal else None
                pairs.append(item)
    for sample in range(native['samples']):
        for key in ['receivedRho','receivedE']:assert float(rows[(sample,'old')][key]).hex()==float(rows[(sample,'new')][key]).hex()
    maxima={mode:max(abs(e['errorK']) for e in errors if e['mode']==mode) for mode in ['old','new']}
    rms={mode:float(np.sqrt(np.mean([p[mode+'ErrorK']**2 for p in pairs]))) for mode in ['old','new']}
    precisionPassed=max(maxima.values())<=2e-10;gain=bool(maxima['new']<=.5*maxima['old'] and rms['new']<=.5*rms['old'])
    acceptance=all(x['acceptancePassed'] for x in controlsReport)
    save(out/'state_errors.json',errors);save(out/'pair_errors.json',pairs);save(out/'controls.json',controlsReport)
    report=dict(passed=True,precisionPassed=precisionPassed,domainAcceptancePassed=acceptance,measurableGainPassed=gain,candidateAccepted=precisionPassed and acceptance and gain,
        nativeSummary=native,reference='80 decimal digits; global piecewise-integral bisection, 160 iterations; independent of rationalized root formula',
        actualReplaySamples=native['samples'],pairedTemperatureDifferences=len(pairs),maximumTemperatureErrorK=maxima,pairedDeltaErrorRMSK=rms,
        maximumErrorNewToOldRatio=maxima['new']/maxima['old'],pairedRMSErrorNewToOldRatio=rms['new']/rms['old'],
        allActualReferenceRootsUnique=True,ambiguousKnotReferenceControls=[x['id'] for x in controlsReport if x.get('uniqueReference') is False],
        domainFailures=[x for x in controlsReport if not x['acceptancePassed']],entropyTemperatureBitMismatches=native['entropyTemperatureBitMismatches'],
        oldTraceBitwiseOriginal=True,newMethodCounts={method:sum(e['method']==method and e['mode']=='new' for e in errors) for method in set(e['method'] for e in errors if e['mode']=='new')},
        CFDExecuted=False,physicalQualification=False,newQualifiedAreaM2=0,uniqueCauseOfV258FailureProved=False,
        inputSHA256={p.name:sha(p) for p in CASE.glob('*') if p.is_file()},elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report);print(json.dumps({k:report[k] for k in ['passed','precisionPassed','domainAcceptancePassed','measurableGainPassed','maximumTemperatureErrorK','pairedDeltaErrorRMSK','elapsedSeconds']}))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v260 fixed caloric audit closed; follow NEXT_DECISION')
    audit()
