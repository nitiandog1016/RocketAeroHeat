"""仅收集v183小证据，检查并发保护、原失败与MATLAB源码稳定，不删除独有文件。"""
from pathlib import Path
import hashlib,json,shutil,time
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'outputs/v183_validation';DEST=ROOT/'docs/validation/v183/evidence'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((BASE/name).read_text('utf-8-sig'))


def collect():
    start=time.monotonic();assert not DEST.exists()
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks={p:sha(ROOT/Path(p.replace('\\','/')))==h for p,h in protected.items()}
    assert len(checks)==29 and all(checks.values())
    prep=read('small01_prepare/report.json');raw=read('small01_run/report.json');recovered=read('small01_recovery/report.json')
    verify=read('verify01/report.json');matlab=read('explicit01_process.json');parent=read('matlab_parent.json')
    assert prep['passed'] and raw['passed'] is False and raw['commands'][0]['passed']
    assert recovered['passed'] and recovered['newNativeRuns']==0 and verify['passed']
    assert matlab['status']=='terminal' and matlab['exitCode']==0 and all(matlab['sourceStableAfter'].values())
    assert parent['status']=='terminal' and parent['exitCode']==0 and parent['elapsedSeconds']<=90
    entries=[];large=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.suffix=='.pyc':continue
        record=dict(source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
        if record['bytes']>1024**2:large.append(record);continue
        target=DEST/p.relative_to(BASE)
        if target.suffix=='.log':target=target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        assert sha(target)==record['sha256'];entries.append(dict(path=target.relative_to(ROOT).as_posix(),**record))
    sources=['scripts/prepareEilmerRocketCase.py','scripts/prepareEilmerDevelopmentCase.py',
             'scripts/eilmerExplicitProfile.py','scripts/runPreparedEilmerCase.py','scripts/auditEilmerCompletion.py',
             'scripts/auditEilmerExplicitCompletion.py','scripts/reauditEilmerExplicitRun.py','scripts/readEilmerNative.py',
             'tests/verifyEilmerExplicitRun.py','tests/testEilmerExplicitEntry.m']
    m=dict(files=entries,filesCount=len(entries),totalSmallBytes=sum(x['bytes'] for x in entries),largeLocalOnly=large,
           finalSourcesSHA256={p:sha(ROOT/p) for p in sources},protectedFiles=checks,uniqueFilesDeleted=0,newQualifiedAreaM2=0,
           note='Exact prepare/run sources and original failed receipt retained. Recovery is post-run read-only binding, not original exit-time output binding. One six-step native run, no full domain.')
    (DEST.parent/'evidence_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    b=dict(preparationSeconds=prep['elapsedSeconds'],preparationBudgetSeconds=60,
           runEntrySeconds=raw['elapsedSeconds'],runEntryBudgetSeconds=30,nativeSupervisionSeconds=raw['commands'][0]['elapsedSeconds'],
           recoverySeconds=recovered['recoveryElapsedSeconds'],independentVerificationSeconds=verify['elapsedSeconds'],
           sourceContractSeconds=read('source_contract/report.json')['elapsedSeconds'],readOnlyBudgetSeconds=120,
           matlabParentSeconds=parent['elapsedSeconds'],matlabParentLimitSeconds=90,
           matlabInsideTestSeconds=read('matlab01/report.json')['elapsedSeconds'],matlabSourcesAndInputsStable=len(matlab['sourceStableAfter']),
           nativeRuns=1,nativeUpdates=6,fullDomainRuns=0,engineBuilds=0,
           nativePeakGroupRSSBytes=raw['commands'][0]['peakObservedGroupRSSBytes'],nativeCaseBytes=raw['commands'][0]['finalCaseBytes'],
           outputDirectoryLogicalBytes=sum(x['bytes'] for x in entries+large),collectionSeconds=time.monotonic()-start,
           originalAuditFailureRetained=True,sourceAtRunDiffersFromFinalAudit=True,
           uniqueFilesDeleted=0,physicalQualification=False,newQualifiedAreaM2=0,
           timingScope='Program time excludes editing/reading/manual coordination; only this version output directory counted, logical bytes not allocated disk usage.')
    (DEST.parent/'budget_receipt.json').write_text(json.dumps(b,indent=2)+'\n')
    print(json.dumps(dict(files=len(entries),bytes=m['totalSmallBytes'],large=len(large),protected=len(checks),matlabLocked=len(matlab['sourceStableAfter']))))


if __name__=='__main__':collect()
