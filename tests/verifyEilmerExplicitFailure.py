"""真实完整域初始化失败的只读终态审核；失败不得导出为成功结果或再启动。"""
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import sha,export_loads
from auditEilmerCompletion import audit_completion
from inspectEilmerRun import inspect
from runPreparedEilmerCase import execute


def verify(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'outputs/v184_validation';case=base/'real01';rp=base/'real01_run/report.json'
    r=json.loads(rp.read_text());proc=r['commands'][0];parent=json.loads((base/'parent.json').read_text())
    assert not r['passed'] and proc['nativeExitCode']==1 and proc['stopReason'] is None
    assert proc['inputsStable'] and not proc['liveGroupAfterCleanup'] and parent['status']=='terminal'
    assert r['outputBinding']=='at-supervisor-completion-audit-failed'
    assert all(sha(path)==digest for path,digest in r['outputSHA256'].items())
    log=base/'real01_run/process00/stdout.txt';text=log.read_text()
    assert 'OutOfMemoryError' in text and 'UFluidBlock.init_grid_and_flow_arrays' in text
    assert 'Number-type: real' in text and 'Step=' not in text and 'FINAL-STEP' not in text
    snapshots=sorted(p.name for p in (case/'lmrsim/snapshots').iterdir() if p.is_dir())
    history=case/'lmrsim/transient-residuals.dat'
    history_rows=[line for line in history.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')] if history.exists() else []
    assert snapshots==['0000'] and not history_rows
    assert not list((case/'lmrsim/loads').glob('*/*-wall.dat'))
    source=ROOT/'outputs/v178_validation/real01'
    assert sha(case/'lmrsim/snapshots/0000/fluid-0000')==sha(source/'lmrsim/snapshots/0000/fluid-0000')
    status=inspect(case,rp,out/'inspection.json');assert not status['passed'] and not status['completedWallDataAvailable']
    rejected=[]
    def reject(label,fn):
        try:fn()
        except ValueError as e:rejected.append(dict(label=label,reason=str(e)))
        else:raise AssertionError(label)
    reject('native-completion',lambda:audit_completion(case,log))
    reject('failed-run-export',lambda:export_loads(case,case/'lmrsim/loads/0000/blk-0000-bndry-0-wall.dat',out/'must-not.csv',rp))
    reject('closed-prepare',lambda:execute(case,'prepare',out/'must-not-prepare',120,1024))
    reject('closed-run',lambda:execute(case,'run',out/'must-not-run',660,1024,base/'real01_prepare/report.json'))
    assert not any(out.glob('must-not*'))
    result=dict(passed=True,failureConfirmed='native-OutOfMemoryError-during-init_grid_and_flow_arrays',
                nativeExitCode=proc['nativeExitCode'],nativeSupervisionSeconds=proc['elapsedSeconds'],
                entrySeconds=r['elapsedSeconds'],parentSeconds=parent['elapsedSeconds'],
                nativePeakGroupRSSBytes=proc['peakObservedGroupRSSBytes'],addressSpaceLimitMiB=proc['addressSpaceLimitMiB'],
                groupRSSLimitMiB=proc['groupRSSLimitMiB'],peakVirtualSizeKnown=False,allocationLineKnown=False,
                inputCount=len(proc['inputSHA256']),inputsStable=True,initialFieldUnchanged=True,savedSnapshots=snapshots,
                newNativeRuns=0,newIterations=0,newDevelopedFields=0,newWallLoads=0,historyFileHasHeader=history.exists(),historyDataRows=0,
                screen='not-evaluable-no-iterations',rejected=rejected,physicalQualification=False,newQualifiedAreaM2=0,
                sourceSHA256={name:sha(ROOT/name) for name in ['scripts/readEilmerNative.py','scripts/runPreparedEilmerCase.py','scripts/auditEilmerExplicitCompletion.py','tests/verifyEilmerExplicitFailure.py']},
                elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);verify(p.parse_args().out)
