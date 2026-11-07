"""原箭体温域失败只读审核；保留原生退出和监督事后观察的区别。"""
from pathlib import Path
import json,re,shutil,time
from prepareOpenFoamControl import ROOT,sha,save
from auditOpenFoamDevelopment import log_progress
from readOpenFoamNistRocketState import verify_identity,mapped_geometry


def main():
    start=time.monotonic();base=ROOT/'outputs/v202_validation';out=base/'rocket_failure01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    for name in ['readOpenFoamNistRocketState.py','readOpenFoamRocketState.py','auditOpenFoamDevelopment.py']:
        shutil.copyfile(ROOT/'scripts'/name,out/name)
    run=json.loads((base/'rocket_run01/report.json').read_text());admission=json.loads((base/'rocket_run01/admission.json').read_text())
    process=run['process'];assert process['nativeExitCode']==1 and process['stopReason']=='scientific-guard-after-exit'
    assert not run['passed'] and not run['saved'] and all(process['inputStable'].values()) and not process['liveGroupAfterCleanup']
    assert process['logExportPassed'];verify_identity(base/'rocket01',run,admission)
    g,m,ids,owner,surface=mapped_geometry(base/'rocket01')
    log=(base/'rocket_run01/native/stdout.txt').read_text();progress=log_progress(log)
    queries=[float(x) for x in re.findall(r'([\deE.+-]+) out of range 100 to 2500',log)]
    assert len(queries)==1 and 99<queries[0]<100 and progress['announcedSteps']==progress['completedLogSteps']+1
    assert 'FOAM FATAL ERROR' in log and 'FOAM exiting' in log
    times=[]
    for p in (base/'rocket01').iterdir():
        if not p.is_dir():continue
        try:t=float(p.name)
        except ValueError:continue
        if t>0:times.append(t)
    assert not times
    co_error=abs(run['initialMaxCo']/admission['initialMaxCoExpected']-1);assert co_error<5e-6
    result=dict(passed=True,candidatePassed=False,actualSolverDispatchVerified=run['actualSolverDispatchVerified'],
        nativeExitCode=1,nativeFatalPrecededSupervisorObservation=True,supervisorTerminatedRunningSolver=False,
        nativeTemperatureTableRequestsK=queries,completeAbnormalFieldAvailable=False,fullFieldMinimumTemperatureK=None,
        offendingCellUnknown=True,uniqueRootCauseEstablished=False,progress=progress,completeSavedDevelopmentFields=0,
        initialCoRelativeError=co_error,originalWallFaces=len(ids),originalAreaM2=float(g['area'][ids].sum()),
        originalWallExact=True,modelAndConditionIdentityPassed=True,inputCount=len(process['inputStable']),inputsStable=True,
        observationIORetries=process['observationIORetries'],process=process,
        newQualifiedHeatAreaM2=0,physicalQualification=False,numericallyConverged=False,productionReadiness=False,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','candidatePassed','progress','nativeTemperatureTableRequestsK','initialCoRelativeError','elapsedSeconds']}))


if __name__=='__main__':main()
