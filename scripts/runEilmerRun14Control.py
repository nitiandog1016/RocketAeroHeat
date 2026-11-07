"""具名Run14控制的原生准备/唯一有预算运行；复用Eilmer进程监督。"""
from pathlib import Path
import argparse
import json
import shutil
import time
from runEilmerProcess import run,sha,COMMIT
from auditEilmerCompletion import audit_completion


def execute(case,stage,record,seconds,prepared=None,audit=None):
    started=time.monotonic();case=Path(case).resolve();record=Path(record).resolve()
    assert not record.exists();m=json.loads((case/'case_manifest.json').read_text())
    assert m['schema']=='RAH.EilmerRun14Control.v1' and m['sourceCommit']==COMMIT
    locks=[case/'case_manifest.json',Path(__file__),Path(__file__).with_name('auditEilmerCompletion.py')]
    for name,digest in m['filesSHA256'].items():
        assert sha(case/name)==digest;locks.append(case/name)
    if stage=='prepare':
        assert not (case/'lmrsim').exists()
        commands=[['prep-gas','gas.inp','nitrogen.gas'],['lmr','prep-grids','--job=grid.lua'],['lmr','prep-sim','--job=job.lua']]
    else:
        assert prepared and audit
        previous=json.loads(Path(prepared).read_text());check=json.loads(Path(audit).read_text())
        assert previous['passed'] and previous['stage']=='prepare' and previous['case']==str(case)
        assert check['passed'] and check['caseManifestSHA256']==sha(case/'case_manifest.json')
        for path,digest in previous['nativeInputSHA256'].items():assert sha(path)==digest;locks.append(Path(path))
        commands=[['lmr','run','--max-wall-clock='+str(int(seconds)-5)]]
    record.mkdir(parents=True);shutil.copyfile(__file__,record/'executed_entry.py')
    result=dict(stage=stage,case=str(case),caseManifestSHA256=sha(case/'case_manifest.json'),
                passed=False,physicalQualification=False,commands=[])
    try:
        for i,cmd in enumerate(commands):
            remaining=seconds-(time.monotonic()-started);assert remaining>0
            receipt=run(cmd,case,record/('process%02d'%i),remaining,256,locks)
            result['commands'].append(receipt)
            if not receipt['passed']:break
        else:
            if stage=='run':result['nativeCompletion']=audit_completion(case,record/'process00/stdout.txt')
            result['passed']=True
        key='nativeInputSHA256' if stage=='prepare' else 'outputSHA256'
        result[key]={str(p):sha(p) for p in case.rglob('*') if p.is_file()}
    except Exception as error:result['error']=repr(error)
    finally:
        result['elapsedSeconds']=time.monotonic()-started
        (record/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--record',required=True)
    p.add_argument('--stage',choices=['prepare','run'],required=True);p.add_argument('--seconds',type=float,required=True)
    p.add_argument('--prepared');p.add_argument('--audit');a=p.parse_args()
    r=execute(a.case,a.stage,a.record,a.seconds,a.prepared,a.audit)
    print(json.dumps({k:r.get(k) for k in ('passed','stage','elapsedSeconds','error')}));raise SystemExit(0 if r['passed'] else 1)
