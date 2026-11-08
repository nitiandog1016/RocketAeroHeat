"""模块职责：核对准确执行源、保护并发/原资产，关闭本轮入口并归档小型可复核证据。"""
from pathlib import Path
import ast,hashlib,json,re,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v218_validation';DOC=ROOT/'docs/validation/v218'
NAMES=['inspectIndependentReferenceV218.py','extractIndependentReferenceSourcesV218.py',
       'assessIndependentReferenceV218.py','runIndependentReferenceV218.py']

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def save(p,data):p.write_bytes((json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n').encode('utf8'))

def main():
    start=time.monotonic();out=BASE/'final02';out.mkdir(exist_ok=False)
    # 模块一：区分调用入口、任务和被导入的监督模块，逐项与运行回执SHA绑定。
    records={}
    for stage in ('reconstruct','source','assess'):
        directory=BASE/('offline_'+stage);r=json.loads((directory/'report.json').read_text(encoding='utf8'))
        assert r['passed'] and r['inputsStable'] and r['nativeExitCode']==0 and not r['liveGroupAfterCleanup']
        assert r['stopReason'] is None
        caller=BASE/('entry_'+stage)/'executed_entry.py'
        assert sha(caller) in r['inputSHA256'].values()
        if stage=='source':
            task=BASE/'entry_source/executed_task.py';supervisor=directory/'executed_runner.py'
            assert sha(supervisor)==sha(ROOT/'scripts/runEilmerProcess.py')
        else:
            task=directory/'executed_task.py';supervisor=directory/'executed_entry.py'
            assert sha(supervisor)==sha(ROOT/'scripts/runCompleteRegionDecisionV216.py')
        assert sha(task) in r['inputSHA256'].values()
        records[stage]=dict(seconds=r['elapsedSeconds'],sampledRSSBytes=r['peakObservedGroupRSSBytes'],
            caller=caller.relative_to(ROOT).as_posix(),callerSHA256=sha(caller),
            task=task.relative_to(ROOT).as_posix(),taskSHA256=sha(task),
            supervisor=supervisor.relative_to(ROOT).as_posix(),supervisorSHA256=sha(supervisor))
    total=sum(x['seconds'] for x in records.values());assert total<120
    assessment=json.loads((BASE/'assessment01/report.json').read_text(encoding='utf8'))
    assert assessment['passed'] and not assessment['fullNumericalAcceptance'] and not assessment['productionAdmitted']
    for p,h in assessment['inputSHA256'].items():assert sha(ROOT/p)==h,p
    src=json.loads((BASE/'source01/report.json').read_text(encoding='utf8'));assert len(src['sources'])==8
    for row in src['sources']:assert sha(ROOT/row['path'])==row['sha256']
    assert sum(x['bytes'] for x in src['sources'])<16*1024**2
    contract=json.loads((DOC/'route_contract.json').read_text(encoding='utf8'))
    assert contract['goalStatus']=='active' and contract['decision']=='closed-to-production'
    assert all(contract[k] is False for k in ('nativeQueriesAdmitted','buildsAdmitted','meshGenerationAdmitted',
               'CFDAdmitted','fullHighPrecisionReferenceAdmitted','fullNumericalAcceptance'))
    # 模块二：核对全部受保护文件，以及既有母网格和七域资产；不删除或重生数据。
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf8'));assert len(protected)==29
    for p,h in protected.items():assert sha(ROOT/p.replace('\\','/'))==h,p
    prior=json.loads((ROOT/'outputs/v216_validation/final01/report.json').read_text(encoding='utf8'));assets={}
    for p,h in prior['sourceSHA256'].items():
        if p.endswith(('.npz','.mat')):
            assert sha(ROOT/p.replace('\\','/'))==h,p;assets[p]=h
    # 模块三：保留执行时源码，然后关闭CLI并确认拒绝不产生新目录。
    closed=BASE/'closed_entries';closed.mkdir(exist_ok=False)
    for name in NAMES:
        p=ROOT/'scripts'/name;shutil.copyfile(p,closed/name);text=p.read_text(encoding='utf8')
        match=re.search(r"^if __name__\s*==\s*['\"]__main__['\"]\s*:",text,re.M);assert match
        p.write_bytes((text[:match.start()]+"if __name__ == '__main__':\n    raise ValueError('v218评估已完成、冻结生产分支关闭；请按NEXT_DECISION登记新任务。')\n").encode('utf8'))
        ast.parse(p.read_text(encoding='utf8'))
    before={p.name for p in BASE.iterdir()}
    refused=subprocess.run([sys.executable,str(ROOT/'scripts/runIndependentReferenceV218.py'),'assess'],
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=5)
    assert refused.returncode!=0 and b'ValueError' in refused.stderr
    assert before=={p.name for p in BASE.iterdir()}
    (out/'closed_entry_stderr.txt').write_bytes(refused.stderr)
    for name in NAMES+[Path(__file__).name]:ast.parse((ROOT/'scripts'/name).read_text(encoding='utf8'))
    save(out/'report.json',dict(passed=True,meaningOfPassed='Integrity/closure/archive checks only',
        tasks=records,supervisedSeconds=total,sourceOutputBytes=sum(x['bytes'] for x in src['sources']),
        sourceOutputNotCoveredByWorkDirectoryGuard=True,sourceOutputSizeIndependentlyAccepted=True,
        protectedFilesUnchanged=29,originalAssetSHA256=assets,CLIClosed=True,CLIRefusalWithoutNewWorkPassed=True,
        fullNumericalAcceptance=False,CFDAdmitted=False,newQualifiedAreaM2=0,goalComplete=False,
        bookkeepingSeconds=time.monotonic()-start))
    # 模块四：小证据公开；上游完整源码仍留本地，只公开带出处的必要摘录和许可。
    dest=DOC/'evidence';dest.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);item=dict(source=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
        if rel.parts[:2]==('source01','source'):
            local.append(item);continue
        assert p.stat().st_size<1024**2
        q=dest/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==item['sha256'];item['path']=q.relative_to(ROOT).as_posix();rows.append(item)
    license=ROOT/'docs/validation/v209/evidence/GDTK_COPYING';q=dest/'GDTK_COPYING';shutil.copyfile(license,q)
    rows.append(dict(source=license.relative_to(ROOT).as_posix(),path=q.relative_to(ROOT).as_posix(),sha256=sha(q),bytes=q.stat().st_size))
    sources=['AGENTS.md','README.md','.gitattributes','docs/PROJECT_GOALS.md',
        'docs/HYBRID_INDEPENDENT_REFERENCE_DECISION_V218.md']+['scripts/'+name for name in NAMES]
    sources += [Path(__file__).resolve().relative_to(ROOT).as_posix()]
    sources += [p.relative_to(ROOT).as_posix() for p in sorted(DOC.glob('*.md'))]
    sources += [(DOC/'route_contract.json').relative_to(ROOT).as_posix()]
    assert sum(x['bytes'] for x in rows)<2*1024**2
    save(DOC/'evidence_manifest.json',dict(schema='RAH.V218Evidence.v1',files=rows,localAssets=local,
        sourceFiles=sources,sourceSHA256={p:sha(ROOT/p) for p in sources},count=len(rows),bytes=sum(x['bytes'] for x in rows),
        fullNumericalAcceptance=False,CFDAdmitted=False,protectedFilesUnchanged=29,goalComplete=False))
    print(json.dumps(dict(passed=True,evidenceFiles=len(rows),evidenceBytes=sum(x['bytes'] for x in rows),
        supervisedSeconds=total,bookkeepingSeconds=time.monotonic()-start,protectedFilesUnchanged=29)))

if __name__=='__main__':main()
