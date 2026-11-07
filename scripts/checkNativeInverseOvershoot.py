"""一次已知域内真解的原生反算区分，不重跑箭体或小CFD。"""
from pathlib import Path
import json,shutil,time,re
from prepareOpenFoamControl import ROOT,sha,save


def main():
    from runOpenFoamProcess import run
    start=time.monotonic();out=ROOT/'outputs/v203_validation/scalar01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_entry.py')
    old=ROOT/'outputs/v202_validation/native01';prior=json.loads((old/'report.json').read_text());binary=old/'probe/rahNistProbe'
    assert prior['passed'] and sha(binary)==prior['queryBinarySHA256']
    sample=out/'samples.txt';sample.write_text('1197.0279007383137 100.01 180\n')
    review=json.loads((ROOT/'outputs/v203_validation/source01/report.json').read_text());assert review['passed']
    row=next(r for r in review['cases'] if r['trueTemperatureK']==100.01 and r['seedK']==180)
    assert row['result']=='trial-out-of-range'
    locks=[Path(__file__),binary,old/'mixture.dict',sample,ROOT/'docs/validation/v203/SCALAR_PLAN.md']
    for r in review['files']:
        p=Path('/opt/openfoam13')/r['path'];assert sha(p)==r['sha256'];locks.append(p)
    save(out/'admission.json',dict(seconds=10,automaticContinuation=0,CFDRuns=0,inputSHA256={str(p):sha(p) for p in locks}))
    process=run([binary,old/'mixture.dict',sample],out,out/'native',10,locks,16,2048)
    text=(out/'native/stdout.txt').read_text();m=re.search(r'([\deE.+-]+) out of range 100 to 2500',text)
    observed=float(m.group(1)) if m else None
    passed=process['nativeExitCode']==1 and observed is not None and abs(observed-row['history'][0]['trialTemperatureK'])<1e-10
    passed=passed and not process['liveGroupAfterCleanup'] and all(process['inputStable'].values()) and process['logExportPassed']
    result=dict(passed=passed,knownTrueTemperatureK=100.01,seedK=180,targetEnergyJkg=row['targetEnergyJkg'],
        predictedTrialTemperatureK=row['history'][0]['trialTemperatureK'],nativeOutOfRangeQueryK=observed,
        falseRejectionOfInRangeEnergy=passed,exactRocketRootCauseEstablished=False,process=process,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','nativeOutOfRangeQueryK','falseRejectionOfInRangeEnergy','elapsedSeconds']}))


if __name__=='__main__':main()
