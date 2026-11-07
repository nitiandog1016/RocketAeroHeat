# -*- coding: utf-8 -*-
"""具名Euler小控制一次性有限运行；热流未知，拒绝旧RANS/整箭准入。"""
from pathlib import Path
import argparse
import json
import shutil
import time
from runEilmerProcess import run, sha
from readHybridEulerNative import audit_inputs
from eilmerWorkspaceRuntime import audit_workspace_runtime


def execute(admission, record):
    start=time.monotonic();admission=Path(admission);record=Path(record)
    plan=json.loads(admission.read_text());case=Path(plan['case'])
    assert plan['schema']=='RAH.HybridEulerControlAdmission.v1'
    assert plan['budget']==dict(nativeSeconds=75,entrySeconds=90,memoryMiB=4096,outputMiB=128,
                                steps=400,cpuCount=2,automaticContinuation=0)
    for p,digest in plan['evidenceSHA256'].items():assert sha(p)==digest,'准入来源变化: '+p
    assert not record.exists() and str(record.resolve())==plan['record']
    audit,grid,mapping=audit_inputs(case,plan['preparation'])
    assert audit['cells']==2400 and audit['ghostCells']==2720
    root=Path(__file__).resolve().parents[1];runtime=root/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json'
    package,locks=audit_workspace_runtime(dict(explicit=dict(profile='rah-explicit-workspace-v1',
                                        runtimePackage=str(runtime),runtimePackageSHA256=plan['runtimeSHA256'])))
    with (case/'hybrid-run-attempt.json').open('x') as f:
        json.dump(dict(admissionSHA256=sha(admission),automaticContinuation=0),f)
    record.mkdir(parents=True);shutil.copyfile(__file__,record/'executed_entry.py')
    sources=[Path(__file__).with_name(n) for n in ['readHybridEulerNative.py','eilmerAirThermo.py',
              'readEilmerNative.py','hybridExternalState.py']]
    for p in sources:shutil.copyfile(p,record/p.name)
    locked=list(json.loads(Path(plan['preparation']).read_text())['nativeInputSHA256'])
    locks += locked+list(plan['evidenceSHA256'])+sources+[admission,Path(__file__),case/'hybrid-run-attempt.json']
    r=run([package['binary'],'--max-cpus=2','--max-wall-clock=75'],case,record/'native',
          90-(time.monotonic()-start),128,locks,4096,4096)
    # 原生保存是否完整由独立读取器审核，不以exit0宣告收敛。
    result=dict(schema='RAH.HybridEulerControlRun.v1',processPassed=r['passed'],
        elapsedSeconds=time.monotonic()-start,admissionSHA256=sha(admission),
        outputSHA256={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()},
        physicalQualification=False,convergencePassed=None)
    (record/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(processPassed=r['passed'],nativeExitCode=r['nativeExitCode'],elapsedSeconds=result['elapsedSeconds'])))
    return r['passed']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--admission',required=True);p.add_argument('--record',required=True)
    a=p.parse_args();raise SystemExit(0 if execute(a.admission,a.record) else 1)
