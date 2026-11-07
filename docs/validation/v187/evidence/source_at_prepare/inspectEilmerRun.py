"""只读生成求解状态，明确区分启动器exit0、原生失败和没有壁面结果。"""
from pathlib import Path
import argparse,hashlib,json,time
from auditEilmerCompletion import audit_completion,completion_log_for_receipt


def inspect(case,receipt,out):
    start=time.monotonic();case=Path(case).resolve();receipt=Path(receipt);out=Path(out)
    assert not out.exists()
    r=json.loads(receipt.read_text());manifest=case/'case_manifest.json'
    assert r['case']==str(case) and r['stage']=='run'
    assert r['caseManifestSHA256']==hashlib.sha256(manifest.read_bytes()).hexdigest()
    result=dict(passed=False,status='failed',solverCompleted=False,physicsQualified=False,newQualifiedAreaM2=0,
        launcherReportedPassed=r['passed'],launcherExitCode=r['commands'][0]['nativeExitCode'],
        processGroupPeakRSSBytes=r['commands'][0].get('peakObservedGroupRSSBytes'),
        liveProcessGroupAfterCleanup=r['commands'][0].get('liveGroupAfterCleanup'),
        sourceRunReceiptSHA256=hashlib.sha256(receipt.read_bytes()).hexdigest(),
        completedWallDataAvailable=False)
    try:
        log=completion_log_for_receipt(case,receipt,r)
        result['nativeCompletion']=audit_completion(case,log)
        m=json.loads(manifest.read_text())
        if m['condition'].get('explicit',{}).get('profile')=='rah-explicit-workspace-v1':
            from eilmerWorkspaceRuntime import audit_workspace_execution
            audit_workspace_execution(m['condition'],r)
        if m['condition'].get('explicit',{}).get('profile') in ('rah-explicit-dense-v1','rah-explicit-dense-development-v1'):
            from eilmerDenseRuntime import audit_dense_execution
            audit_dense_execution(m['condition'],r)
        assert r['passed'] and all(c['passed'] for c in r['commands']),'原监督拒绝'
        result.update(passed=True,status='completed-research',solverCompleted=True,completedWallDataAvailable=True)
    except (ValueError,AssertionError,FileNotFoundError) as e:result['failureReason']=str(e)
    result['elapsedSeconds']=time.monotonic()-start
    out.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--receipt',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(inspect(a.case,a.receipt,a.out)))
