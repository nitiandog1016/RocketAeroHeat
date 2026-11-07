"""完整域新存储初始化失败的只读终态检查；实际容量、来源和失败拒绝独立核对。"""
from pathlib import Path
import argparse
import json
import re
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runEilmerProcess import sha
from readEilmerNative import export_loads,audit_run_scope
from auditEilmerCompletion import audit_completion
from eilmerExplicitAdmission import audit_admission
from eilmerDenseRuntime import audit_dense_execution
from inspectEilmerRun import inspect
from runPreparedEilmerCase import execute
from prepareEilmerDevelopmentCase import prepare


def verify(out):
    start=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_verifier.py')
    result=dict(passed=False,newNativeRuns=0,physicalQualification=False,newQualifiedAreaM2=0)
    try:
        base=ROOT/'outputs/v186_validation';case=base/'real01';rp=base/'real01_run/report.json'
        r=json.loads(rp.read_text());proc=r['commands'][0];m=json.loads((case/'case_manifest.json').read_text())
        parent=json.loads((base/'run_parent.json').read_text())
        assert not r['passed'] and proc['nativeExitCode']==1 and proc['stopReason'] is None
        assert proc['inputsStable'] and not proc['liveGroupAfterCleanup'] and parent['status']=='terminal' and parent['exitCode']==1
        assert r['outputBinding']=='at-supervisor-completion-audit-failed'
        assert all(sha(p)==h for p,h in r['outputSHA256'].items())
        assert audit_dense_execution(m['condition'],r)
        log=base/'real01_run/process00/stdout.txt';text=log.read_text()
        assert 'OutOfMemoryError' in text and 'UFluidBlock.init_grid_and_flow_arrays' in text
        assert 'Step=' not in text and 'FINAL-STEP' not in text
        assert sorted(p.name for p in (case/'lmrsim/snapshots').iterdir() if p.is_dir())==['0000']
        assert not list((case/'lmrsim/loads').glob('*/*-wall.dat'))
        history=case/'lmrsim/transient-residuals.dat'
        assert not [s for s in history.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
        # 容量依据已绑定的原几何边界计数，不为失败重复解析大网格。
        a,_=audit_admission(m,case)
        geom=json.loads(Path(a['evidence']['sourceGeometry']['path']).read_text())
        matches=re.findall(r'^RAH-STORAGE cells=(\d+) ghosts=(\d+) neq=(\d+) nftl=(\d+) scalarBytes=(\d+) stageLength=(\d+) payloadBytes=(\d+)$',text,re.M)
        assert len(matches)==1
        cells,ghosts,neq,levels,width,length,payload=map(int,matches[0])
        assert cells==geom['cells']==705490 and ghosts==sum(geom['boundaryCounts'].values())==26324
        assert (neq,levels,width)==(7,4,8) and length==(cells+ghosts)*neq and payload==2*levels*length*width==327852672
        source=ROOT/'outputs/v178_validation/real01'
        assert sha(case/'lmrsim/snapshots/0000/fluid-0000')==sha(source/'lmrsim/snapshots/0000/fluid-0000')
        status=inspect(case,rp,out/'inspection.json');assert not status['passed'] and not status['completedWallDataAvailable']
        rejected=[]
        def reject(label,fn):
            try:fn()
            except ValueError as error:rejected.append(dict(label=label,reason=str(error)))
            else:raise AssertionError(label)
        reject('native-completion',lambda:audit_completion(case,log))
        reject('failed-export',lambda:export_loads(case,case/'lmrsim/loads/0000/blk-0000-bndry-0-wall.dat',out/'must-not.csv',rp))
        reject('closed-run',lambda:execute(case,'run',out/'must-not-run',660,1024,base/'real01_prepare/report.json'))
        reject('closed-prepare',lambda:prepare(source,ROOT/'outputs/v178_validation/real01_prepare/report.json',
            ROOT/'outputs/v178_validation/real01_geometry.json',out/'must-not-case',out/'must-not-prepare',m['condition'],120,1024))
        reject('closed-scope',lambda:audit_run_scope(m,case))
        assert not list(out.glob('must-not*'))
        result.update(passed=True,failureConfirmed='native-OOM-in-grid-and-flow-initialization-after-dense-array-allocation',
            nativeSupervisionSeconds=proc['elapsedSeconds'],entrySeconds=r['elapsedSeconds'],parentSeconds=parent['elapsedSeconds'],
            nativeExitCode=1,nativePeakGroupRSSBytes=proc['peakObservedGroupRSSBytes'],
            inputCount=len(proc['inputSHA256']),inputsStable=True,allRunOutputSHA256Preserved=True,
            initialFieldUnchanged=True,historyDataRows=0,newIterations=0,newDevelopedFields=0,newWallLoads=0,
            observedStorage=dict(cells=cells,ghosts=ghosts,equations=neq,timeLevels=levels,arrayPayloadBytes=payload,
                                 previousArrayPayloadBytes=payload*levels,arrayPayloadReductionBytes=payload*(levels-1)),
            failedAllocationLineKnown=False,peakVirtualSizeKnown=False,screen='not-evaluable-no-iterations',rejected=rejected,
            sourceSHA256={name:sha(ROOT/name) for name in ['scripts/readEilmerNative.py','scripts/eilmerDenseRuntime.py',
                'scripts/runPreparedEilmerCase.py','scripts/auditEilmerExplicitCompletion.py','tests/verifyEilmerDenseFailure.py']})
    except Exception as error:
        result['error']=repr(error);raise
    finally:
        result['elapsedSeconds']=time.monotonic()-start
        (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);verify(p.parse_args().out)
