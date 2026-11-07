"""只读恢复显式接口审核误拒；原生失败、预算停止及原始失败回执不得改写。"""
from pathlib import Path
import argparse,copy,json,shutil,time
from runEilmerProcess import sha
from auditEilmerCompletion import audit_completion
from readEilmerNative import need


def reaudit(case, original, out):
    start=time.monotonic();case=Path(case).resolve();original=Path(original).resolve();out=Path(out).resolve()
    need(not out.exists(),'恢复回执不覆盖')
    r=json.loads(original.read_text());m=json.loads((case/'case_manifest.json').read_text())
    need(r['stage']=='run' and r['case']==str(case) and 'explicit' in m['condition'],'恢复来源错误')
    need(r['passed'] is False and r.get('error')=="ValueError('wrong-explicit-dispatch')",'仅限v183直接入口提示误拒')
    need(r['caseManifestSHA256']==sha(case/'case_manifest.json'),'恢复工况变化')
    need(len(r['commands'])==1,'恢复原生任务数不符')
    command=r['commands'][0]
    need(command['passed'] and command['argv'][0]=='lmr-run' and command['nativeExitCode']==0 and
         command['stopReason'] is None and command['inputsStable'] and not command['liveGroupAfterCleanup'],
         '不能恢复原生异常/预算退出')
    for name,digest in command['inputSHA256'].items():
        path=Path(name)
        if path.is_relative_to(case):need(sha(path)==digest,'恢复时原生输入改变')
    # 原错误发生在输出SHA封存之前；明确绑定当前只读审核时刻，不能冒称退出时逐字证据。
    before={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()}
    if r.get('outputSHA256'):need(before==r['outputSHA256'],'原已封存输出变化')
    need(sum(p.stat().st_size for p in case.rglob('*') if p.is_file())==command['finalCaseBytes'],'恢复时目录大小改变')
    log=original.parent/'process00/stdout.txt';completion=audit_completion(case,log)
    need(before=={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()},'只读审核期间输出变化')
    out.mkdir(parents=True);(out/'process00').mkdir();shutil.copyfile(log,out/'process00/stdout.txt')
    for name in ['reauditEilmerExplicitRun.py','auditEilmerExplicitCompletion.py','auditEilmerCompletion.py']:
        shutil.copyfile(Path(__file__).with_name(name),out/name)
    recovered=copy.deepcopy(r)
    recovered.update(passed=True,nativeCompletion=completion,outputSHA256=before,
                     outputBinding='post-run-read-only-recovery',originalRunReceiptSHA256=sha(original),
                     originalRunReceipt=str(original),originalError=r['error'],
                     recoveryElapsedSeconds=time.monotonic()-start,newNativeRuns=0,
                     recoverySourceSHA256={p.name:sha(p) for p in out.glob('*.py')})
    recovered.pop('error',None)
    (out/'report.json').write_text(json.dumps(recovered,indent=2)+'\n')
    return recovered


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--original',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();r=reaudit(a.case,a.original,a.out)
    print(json.dumps({k:r[k] for k in ['passed','nativeCompletion','recoveryElapsedSeconds']}))
