"""只读审核v198真实停止与资源；不重跑CFD补日志，不伪造缺失末场。"""
from pathlib import Path
import json,re,time
from prepareOpenFoamControl import ROOT,BASE,sha,save

def main():
    start=time.monotonic();out=BASE/'failure_audit01';out.mkdir(exist_ok=False)
    runs=[]
    for n in ['control_run01','control_run02']:
        directory=BASE/n;r=json.loads((directory/'report.json').read_text());p=r['process']
        assert p['nativeExitCode']==-15 and p['stopReason']=='scientific-guard'
        assert not p['liveGroupAfterCleanup'] and all(p['inputStable'].values())
        assert r['completed600'] is False
        for name in ['runOpenFoamControl.py','runOpenFoamProcess.py','openFoamFields.py','openFoamAirThermo.py']:
            key='/mnt/e/RocketAeroHeat/scripts/'+name;assert sha(directory/name)==p['inputSHA256'][key]
        text=(directory/'native/stdout.txt').read_text(encoding='utf-8')
        queries=[];step=None;completed=[]
        for line in text.splitlines():
            m=re.match(r'Time = (\d+)s',line)
            if m:step=int(m.group(1))
            if line.startswith('ExecutionTime ='):completed.append(step)
            m=re.search(r'out of temperature range 200 -> 2500;\s+T = ([^\s]+)',line)
            if m:queries.append(dict(step=step,requestedTemperatureK=float(m.group(1))))
        case=BASE/('development01' if n=='control_run01' else 'development02')
        saved=[x.name for x in case.iterdir() if x.is_dir() and x.name.isdigit() and x.name!='0']
        assert not saved
        runs.append(dict(name=n,processPassed=p['passed'],nativeExitCode=p['nativeExitCode'],
            guardReason=p['guardReason'],supervisorSeconds=p['elapsedSeconds'],entrySeconds=r['elapsedSeconds'],
            peakGroupRSSBytes=p['peakObservedGroupRSSBytes'],lockedInputs=len(p['inputSHA256']),
            allInputsStable=True,groupEmpty=True,lastAnnouncedStep=step,lastCompletedLogStep=max(completed,default=None),
            warningQueries=queries,savedDevelopmentFrames=saved,normalCompleted600=False))
    assert runs[0]['lastCompletedLogStep'] is None and runs[1]['warningQueries'][0]['step']==39
    assert runs[1]['warningQueries'][0]['requestedTemperatureK']==181.36899944325523
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    checks={p:sha(ROOT/p)==h for p,h in protected.items()};assert len(checks)==29 and all(checks.values())
    result=dict(auditPassed=True,candidateAccepted=False,profile='of13-nasa7-tetra-lts-600-v1',runs=runs,
        originalResidualThresholdPassed=False,diagonalSolverZeroIsNotEquationConvergence=True,
        warningMeans='native inverse-temperature trial was outside range; actual transient full state/location/cause unknown',
        finalFieldAvailable=False,wallPressureValidated=False,newQualifiedHeatAreaM2=0,
        originalRocketFaces=19752,originalRocketAreaM2=1.9646292565715553,CFDOnRealRocket=0,
        protectedFiles=checks,elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    save(out/'report.json',result)
    print(json.dumps(dict(auditPassed=True,candidateAccepted=False,runs=runs,protectedFiles=len(checks),elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
