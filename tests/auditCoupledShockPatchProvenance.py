"""核对每次真实MATLAB进程的计算源码版本、显式输入及视图所用资源。"""

from pathlib import Path
import argparse
import json
import subprocess
from auditThermalPressureSurfaceFlow import sha

ROOT=Path(__file__).resolve().parents[1]
JOBS={'coupledPatchSmoke100':0,'coupledPatchFull100':1,'coupledPatchScope100':1,
      'coupledPatchNormalized100':1,'coupledPatchDiagnostics100':0,'coupledPatchView100':0,'coupledPatchViewCamera100':0}
WEB={'web/geometry.html','web/geometry.css','web/geometry.js','web/rocket-sensitivity.css',
     'web/shock-layer-study.html','web/shock-layer-study.css','web/shock-layer-study.js','web/shock-layer-geometry.js'}


def audit(folder,target):
    """模块用途：实际计算闭包与无关并发修改分开；失败终态和原门槛不得改写。"""
    assert not target.exists();archives={sha(p):p for p in (ROOT/'docs/validation/v100/raw').rglob('*.raw')}
    versions={};processes=[]
    for job,expected_exit in JOBS.items():
        path=folder/(job+'_process.json');d=json.loads(path.read_text(encoding='utf-8-sig'))
        assert d['status']=='terminal' and d['exitCode']==expected_exit
        keys={p for p in d['sourceSHA256'] if p.endswith('.m')}|set(d['explicitInputs'])
        if 'View' in job:
            # 导出器没有调用tests/*.m；并发外部工作新增/改写的测试不属于导出计算闭包。
            # 原全快照及其SHA保留，但不能声称未归档的无关测试版本仍有原字节。
            keys={p for p in keys if not p.startswith('tests/')}|set(d['explicitInputs'])|WEB
        assert all(d['sourceStableAfter'][p] for p in keys)
        for relative in keys:
            digest=d['sourceSHA256'][relative];key=(relative,digest)
            if key in versions:continue
            current=ROOT/relative
            if current.is_file() and sha(current)==digest:resolved=current.relative_to(ROOT).as_posix()
            elif digest in archives:resolved=archives[digest].relative_to(ROOT).as_posix()
            else:
                old=subprocess.run(['git','show','27c6ccf4ebaa1d8ac434a58ac6e044450cb60fae:'+relative],cwd=ROOT,capture_output=True)
                assert old.returncode==0 and __import__('hashlib').sha256(old.stdout).hexdigest()==digest,(job,relative,digest)
                resolved='git:27c6ccf4ebaa1d8ac434a58ac6e044450cb60fae:'+relative
            versions[key]=dict(path=relative,sha256=digest,exactVersionAt=resolved)
        unrelated=[p for p,v in d['sourceStableAfter'].items() if not v and p not in keys]
        processes.append(dict(job=job,status='terminal',exitCode=expected_exit,engineProcessIDs=d['engineProcessIDs'],
            calculationAndExplicitInputsStable=True,wholeSnapshotStable=all(d['sourceStableAfter'].values()),
            unrelatedConcurrentChanges=unrelated,sourceRecordSHA256=sha(path),
            versionScope='all non-test MATLAB sources, actual web resources and explicit inputs' if 'View' in job else 'all MATLAB sources and explicit inputs'))
    report=dict(status='passed-exact-calculation-source-epochs-and-terminal-process-audit',processes=processes,
                sourceVersions=list(versions.values()),versionCount=len(versions),oldFailuresPreserved=True,
                scopeNote='首轮只读来源审计要求导出期间全部无关测试的历史原文，因并发新增测试无旧副本而拒绝；改为实际导出闭包，保留原全快照，不声称无关测试原文已归档。',
                completePhysicalValidation=False,newProductionSolvedAreaM2=0)
    target.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],processes=len(processes),exactVersions=len(versions),
                          backgroundEditEpochs=sum(not p['wholeSnapshotStable'] for p in processes))))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder',type=Path);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();audit(a.folder,a.report)
