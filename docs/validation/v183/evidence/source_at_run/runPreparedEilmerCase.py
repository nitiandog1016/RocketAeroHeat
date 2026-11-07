"""固定Eilmer输入的准备/有限运行入口；由MATLAB或WSL调用，原网格始终只读。"""
from pathlib import Path
import argparse
import json
import shutil
import time
from runEilmerProcess import run, sha, COMMIT
from readEilmerNative import audit_prepared,audit_numerics,audit_run_scope
from auditEilmerCompletion import audit_completion


def execute(case, stage, record, seconds, output_mib, prepared_receipt=None):
    started=time.monotonic()
    case, record = Path(case).resolve(), Path(record).resolve()
    assert not record.exists(), '已有回执不覆盖'
    manifest_path=case/'case_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    assert manifest['schema']=='RAH.EilmerPreparedCase.v1' and manifest['sourceCommit']==COMMIT
    audit_run_scope(manifest)
    mesh=Path(manifest['mesh']['nativeSourcePath'])
    assert sha(mesh)==manifest['mesh']['sourceSHA256'], '原网格身份发生变化'
    locks=[manifest_path,mesh,Path(__file__),Path(__file__).with_name('readEilmerNative.py'),Path(__file__).with_name('auditEilmerCompletion.py')]
    if 'explicit' in manifest['condition']:
        locks += [Path(__file__).with_name(n) for n in ['eilmerExplicitProfile.py','auditEilmerExplicitCompletion.py']]
    for name,digest in manifest['filesSHA256'].items():
        assert sha(case/name)==digest, '准备输入发生变化: '+name
        locks.append(case/name)
    if stage=='prepare':
        assert not (case/'lmrsim').exists(), '拒绝覆盖已有原生目录'
        commands=[['prep-gas','gas.inp','rah-air.gas'],
                  ['lmr','prep-grids','--job=grid.lua'],['lmr','prep-sim','--job=job.lua']]
    else:
        assert prepared_receipt is not None, '运行必须指定已通过的准确准备回执'
        previous=json.loads(Path(prepared_receipt).read_text())
        assert previous['stage']=='prepare' and previous['passed'] and previous['case']==str(case)
        for path,digest in previous['nativeInputSHA256'].items():
            assert sha(path)==digest, '原生输入来源不符: '+path
            locks.append(Path(path))
        # exit0的准备过程不能替代几何/热壁合同；负体积或错CFL不得进入求解。
        if 'development' in manifest['condition'] or 'explicit' in manifest['condition']:
            assert previous.get('uniformInitialFieldBitwiseIdentical') and previous.get('nativeGridReused')
            reuse=manifest['reusedPreparation'];geometry=Path(reuse['geometryAudit'])
            assert sha(geometry)==reuse['geometryAuditSHA256']
            audit=json.loads(geometry.read_text());assert audit['passed'] and audit['nativeGridSHA256']==sha(case/'lmrsim/grid/grid-0000')
            if 'explicit' in manifest['condition']:
                from eilmerExplicitProfile import audit_explicit_config
                audit_explicit_config(json.loads((case/'lmrsim/config').read_text()),manifest['condition'],
                                      json.loads((case/'lmrsim/control').read_text()))
            else:
                audit_numerics(json.loads((case/'lmrsim/nk-config').read_text()),manifest['condition'])
            locks.append(geometry)
        else:audit_prepared(case)
        # 大场终态保存和输入核验需要余量；原生预算与整个阶段预算分别记录。
        reserve=30 if manifest['mesh']['cells']>100000 else 5
        assert seconds>reserve
        native_seconds=int(seconds-(time.monotonic()-started)-reserve)
        assert native_seconds>0,'来源核验后没有剩余原生预算'
        if manifest['condition'].get('development',{}).get('profile')=='rah-nk-real-v1' or 'explicit' in manifest['condition']:
            # 直接监督固定实数二进制，取得真实求解退出码，绕开system()包装状态。
            commands=[['lmr-run','--max-cpus=2','--max-wall-clock='+str(native_seconds)]]
        else:
            commands=[['lmr','run','--max-wall-clock='+str(native_seconds)]]
    record.mkdir(parents=True)
    shutil.copyfile(__file__,record/'executed_entry.py')
    result=dict(stage=stage,case=str(case),sourceCommit=COMMIT,caseManifestSHA256=sha(manifest_path),
                physicsQualified=False,newQualifiedAreaM2=0,passed=False,commands=[])
    try:
        for i,command in enumerate(commands):
            remaining=seconds-(time.monotonic()-started)
            assert remaining>0, '本阶段总预算耗尽'
            report=run(command,case,record/('process%02d'%i),remaining,output_mib,locks)
            result['commands'].append(report)
            if not report['passed']:
                break
        else:
            if stage=='run':
                result['nativeCompletion']=audit_completion(case,record/'process00/stdout.txt')
            result['passed']=True
            if stage=='prepare':
                # 原生配置、网格、初场和热物性绑定；运行不从其他引擎或旧场接续。
                result['nativeInputSHA256']={str(p):sha(p) for p in case.rglob('*') if p.is_file()}
            else:
                result['outputSHA256']={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()}
                result['outputBinding']='at-supervisor-completion'
    except Exception as error:
        result['error']=repr(error)
    finally:
        result['elapsedSeconds']=time.monotonic()-started
        (record/'report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',required=True);p.add_argument('--stage',choices=['prepare','run'],required=True)
    p.add_argument('--record',required=True);p.add_argument('--seconds',type=float,required=True)
    p.add_argument('--output-mib',type=int,default=128);p.add_argument('--prepared-receipt')
    a=p.parse_args()
    r=execute(a.case,a.stage,a.record,a.seconds,a.output_mib,a.prepared_receipt)
    print(json.dumps({k:r.get(k) for k in ('stage','passed','elapsedSeconds','error')}))
    raise SystemExit(0 if r['passed'] else 1)
