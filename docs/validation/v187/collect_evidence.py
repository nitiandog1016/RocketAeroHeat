"""汇集本轮实际构建/相容/失败，小证据逐字保存，大数据仅登记。"""
from pathlib import Path
import hashlib,json,shutil,time

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'outputs/v187_validation'
DEST=ROOT/'docs/validation/v187/evidence'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()


def read(name):
    return json.loads((BASE/name).read_text('utf-8-sig'))


def collect():
    """模块一：检查真实终态、准确阶段源和保护文件，不重新计算。"""
    start=time.monotonic();assert not DEST.exists()
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks={p:sha(ROOT/p.replace('\\','/'))==h for p,h in protected.items()}
    assert len(checks)==29 and all(checks.values())
    build=read('build_result.json');small=read('small01_run/report.json');verify=read('verify01/report.json')
    prep=read('full/real01_prepare/report.json');run=read('full/real01_run/report.json');proc=run['commands'][0]
    fail=read('full/failure01/report.json');closed=read('full/closed_admission01/report.json')
    assert build['passed'] and build['sharedFilesStable'] and build['process']['nativeExitCode']==0
    assert small['passed'] and verify['passed'] and verify['comparedScalars']==658 and verify['maxScaledError']==0
    assert all(verify['fieldAndWallBytesIdentical'].values()) and len(verify['rejected'])==18
    assert sha(BASE/'verify01/executed_verifier.py')==verify['sourceSHA256']['tests/verifyEilmerWorkspaceRuntime.py']
    assert prep['passed'] and prep['uniformInitialFieldBitwiseIdentical'] and prep['nativeGridReused']
    assert not run['passed'] and proc['nativeExitCode']==1 and proc['stopReason'] is None
    assert proc['inputsStable'] and len(proc['inputSHA256'])==46 and not proc['liveGroupAfterCleanup']
    assert fail['passed'] and fail['newIterations']==0 and fail['historyDataRows']==0 and len(fail['rejected'])==5
    assert closed['passed'] and len(closed['rejected'])==16
    for report in [fail,closed]:
        assert all(sha(ROOT/n)==h for n,h in report['sourceSHA256'].items())
    parents={name:read(name+'.json') for name in ['scope_parent','scope02_parent','scope03_parent','build_parent',
        'prepare_parent','run_parent','verify_parent','full/admission_parent','full/prepare_parent',
        'full/run_parent','full/closed_admission_parent','full/verify_parent']}
    for name,r in parents.items():
        assert r['status']=='terminal' and r['exitCode']==(1 if name in ['scope_parent','scope02_parent','full/run_parent'] else 0),name
    for receipt,folder in [(small,BASE/'source_at_run'),(run,BASE/'full/source_at_run')]:
        for p in folder.iterdir():
            actual=receipt['commands'][0]['inputSHA256'].get(str(ROOT/'scripts'/p.name))
            if actual is not None:assert sha(p)==actual,p
    storage=read('full/storage_receipt.json')
    assert storage['uniqueFilesDeleted']==0 and storage['allRunOutputSHA256Preserved']
    for record in storage['records']:
        assert (ROOT/record['path']).samefile(ROOT/record['source'])
    # 模块二：保留当前轮所有小证据，三份大网格/初场只登记SHA。
    entries=[];large=[];observed={}
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.suffix=='.pyc':continue
        item=dict(source=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path))
        observed[str(path)]=item['sha256']
        if item['bytes']>1024**2:large.append(item);continue
        target=DEST/path.relative_to(BASE)
        if target.suffix=='.log':target=target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        assert sha(target)==item['sha256']
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),**item))
    for receipt in [small,run]:assert all(observed[p]==h for p,h in receipt['outputSHA256'].items())
    package=ROOT/'tools/runtime/eilmer-explicit-workspace-v1'
    runtime=json.loads((package/'runtime.json').read_text())
    assert all(sha(p)==h for p,h in runtime['filesSHA256'].items())
    for name in ['runtime.json','buildinfo.json']:
        path=package/name;target=DEST/'runtime_metadata'/name;target.parent.mkdir(exist_ok=True)
        shutil.copyfile(path,target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),source=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path)))
    for name in ['lmr-run','lmr-run.o']:
        path=package/name
        assert path.is_file()
        large.append(dict(source=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path)))
    names=['scripts/buildEilmerWorkspaceRuntime.py','scripts/eilmerWorkspaceRuntime.py','scripts/eilmerExplicitProfile.py',
        'scripts/eilmerExplicitAdmission.py','scripts/readEilmerNative.py','scripts/prepareEilmerDevelopmentCase.py',
        'scripts/runPreparedEilmerCase.py','scripts/auditEilmerExplicitCompletion.py','scripts/inspectEilmerRun.py',
        'tests/verifyEilmerWorkspaceRuntime.py','tests/testEilmerWorkspaceAdmission.py','tests/verifyEilmerWorkspaceFailure.py',
        'patches/eilmer/explicit-workspace-checked-v1.patch']
    manifest=dict(files=entries,filesCount=len(entries),totalSmallBytes=sum(r['bytes'] for r in entries),largeLocalOnly=large,
        finalSourcesSHA256={n:sha(ROOT/n) for n in names},protectedFiles=checks,allRunOutputSHA256Preserved=True,
        physicalQualification=False,newQualifiedAreaM2=0,uniqueFilesDeleted=0,originalFailuresRetained=True)
    (DEST.parent/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # 模块三：分阶段程序/父时间分开，不虚构一个包住所有工作的总时间。
    budget=dict(parentReceipts={name:dict(elapsedSeconds=r['elapsedSeconds'],exitCode=r['exitCode'],limitSeconds=r['parentLimitSeconds']) for name,r in parents.items()},
        layoutCompilerInvocations=3,layoutFailures=2,layoutLinkedExecutables=0,
        buildSupervisionSeconds=build['process']['elapsedSeconds'],buildEntrySeconds=build['elapsedSeconds'],engineBuilds=1,
        compiledDSourceCount=build['compiledDSourceCount'],buildInputs=len(build['process']['inputSHA256']),
        smallNativeSupervisionSeconds=small['commands'][0]['elapsedSeconds'],smallEntrySeconds=small['elapsedSeconds'],
        smallVerificationSeconds=verify['elapsedSeconds'],smallNativeRuns=1,
        fullPreparationSeconds=prep['elapsedSeconds'],fullNativeSupervisionSeconds=proc['elapsedSeconds'],fullRunEntrySeconds=run['elapsedSeconds'],
        fullVerificationSeconds=fail['elapsedSeconds'],fullDomainAttempts=1,fullNativeUpdates=0,
        fullNativePeakGroupRSSBytes=proc['peakObservedGroupRSSBytes'],addressSpaceLimitMiB=8192,groupRSSLimitMiB=8192,
        fullCapacityPartialLowerBoundBytes=read('capacity02/report.json')['partialNecessaryLowerBoundBytes'],
        nominalRemovedWorkspacePayloadBytes=read('capacity02/report.json')['knownRemovalPayloadBytes'],
        fullMemoryFitVerified=False,matlabRuns=0,guiRuns=0,newQualifiedAreaM2=0,
        aliasBytesDeduplicated=storage['aliasBytesDeduplicated'],uniqueFilesDeleted=0,
        postprocessingSharedParent=False,collectionSeconds=time.monotonic()-start,
        timingScope='Nested process/entry/parent times are not additive; source reading, editing and coordination excluded. No shared parent for all postprocessing.')
    (DEST.parent/'budget_receipt.json').write_text(json.dumps(budget,indent=2)+'\n')
    print(json.dumps(dict(files=len(entries),smallBytes=manifest['totalSmallBytes'],largeFiles=len(large),protected=29,collectionSeconds=budget['collectionSeconds'])))


if __name__=='__main__':collect()
