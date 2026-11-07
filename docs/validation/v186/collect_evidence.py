"""收集v186准确失败与工作区依据；仅遍历本轮目录，大场只登记，不删独有数据。"""
from pathlib import Path
import hashlib,json,shutil,time

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'outputs/v186_validation'
DEST=ROOT/'docs/validation/v186/evidence'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):h.update(block)
    return h.hexdigest()


def read(name):
    return json.loads((BASE/name).read_text('utf-8-sig'))


def collect():
    """模块一：核对原生终态、关闭后检查和保护文件，不重新求解或解析大网格。"""
    start=time.monotonic();assert not DEST.exists()
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks={p:sha(ROOT/p.replace('\\','/'))==h for p,h in protected.items()}
    assert len(checks)==29 and all(checks.values())
    prep=read('real01_prepare/report.json');run=read('real01_run/report.json');proc=run['commands'][0]
    failure=read('failure01/report.json');admission=read('admission01/report.json');closed=read('closed_admission01/report.json')
    workspace=read('workspace01/report.json');storage=read('storage_receipt.json')
    parents={k:read(k+'_parent.json') for k in ['prepare','run','verify']}
    assert all(p['status']=='terminal' and p['exitCode']==(1 if k=='run' else 0) for k,p in parents.items())
    assert prep['passed'] and prep['nativeGridReused'] and prep['uniformInitialFieldBitwiseIdentical']
    assert run['passed'] is False and proc['nativeExitCode']==1 and proc['stopReason'] is None
    assert proc['inputsStable'] and not proc['liveGroupAfterCleanup'] and len(proc['inputSHA256'])==44
    assert failure['passed'] and len(failure['rejected'])==5 and failure['historyDataRows']==0
    assert admission['passed'] and len(admission['rejected'])==14 and closed['passed'] and len(closed['rejected'])==15
    for receipt in [failure,closed]:
        assert all(sha(ROOT/n)==h for n,h in receipt['sourceSHA256'].items())
    assert workspace['status']=='unbuilt-workspace-only-draft' and sha(ROOT/workspace['patch'])==workspace['patchSHA256']
    assert all(sha(p)==h for p,h in workspace['sourceSHA256'].items())
    assert all(sha(p)==h for p,h in workspace['runtimeSourceSHA256'].items())
    for record in storage['records']:
        target,source=ROOT/record['path'],ROOT/record['source']
        assert target.samefile(source) and sha(target)==record['sha256']
    assert storage['uniqueFilesDeleted']==0
    entries=[];large=[];observed={}
    # 模块二：小证据逐字保存；硬链接的逻辑大小不算新增磁盘占用。
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.suffix=='.pyc':continue
        record=dict(source=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,sha256=sha(path))
        observed[str(path)]=record['sha256']
        if record['bytes']>1024**2:large.append(record);continue
        target=DEST/path.relative_to(BASE)
        if target.suffix=='.log':target=target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        assert sha(target)==record['sha256']
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),**record))
    assert all(observed[p]==h for p,h in run['outputSHA256'].items())
    sources=['scripts/eilmerExplicitAdmission.py','scripts/eilmerExplicitProfile.py','scripts/eilmerDenseRuntime.py',
        'scripts/prepareEilmerDevelopmentCase.py','scripts/runPreparedEilmerCase.py','scripts/readEilmerNative.py',
        'scripts/auditEilmerExplicitCompletion.py','scripts/inspectEilmerRun.py',
        'tests/testEilmerDenseAdmission.py','tests/verifyEilmerDenseFailure.py','patches/eilmer/explicit-workspace-v1.patch']
    manifest=dict(files=entries,filesCount=len(entries),totalSmallBytes=sum(x['bytes'] for x in entries),largeLocalOnly=large,
        finalSourcesSHA256={p:sha(ROOT/p) for p in sources},protectedFiles=checks,allRunOutputSHA256Preserved=True,
        originalFailureRetained=True,uniqueFilesDeleted=0,newQualifiedAreaM2=0,
        inheritedRuntimeEvidence='docs/validation/v185/evidence/runtime_metadata/runtime.json')
    (DEST.parent/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # 模块三：准备/运行/验证各有父回执，源码分析与整理另计，不虚构共同父时间。
    budget=dict(preparationSeconds=prep['elapsedSeconds'],preparationBudgetSeconds=120,
        preparationParentSeconds=parents['prepare']['elapsedSeconds'],preparationParentLimitSeconds=150,
        nativeSupervisionSeconds=proc['elapsedSeconds'],nativeInternalWallClockLimitSeconds=600,
        runEntrySeconds=run['elapsedSeconds'],runEntryBudgetSeconds=660,
        runParentSeconds=parents['run']['elapsedSeconds'],runParentLimitSeconds=720,
        admissionSeconds=admission['elapsedSeconds'],closedAdmissionSeconds=closed['elapsedSeconds'],
        failureVerificationSeconds=failure['elapsedSeconds'],verificationParentSeconds=parents['verify']['elapsedSeconds'],
        workspaceSourceAuditSeconds=workspace['elapsedSeconds'],storageDeduplicationSeconds=storage['elapsedSeconds'],
        postprocessingBudgetSeconds=180,postprocessingHasSingleEnclosingParent=False,
        nativeRuns=1,fullDomainAttempts=1,nativeUpdates=0,newDevelopedFields=0,newWallLoads=0,
        engineBuilds=0,smallNativeRuns=0,matlabRuns=0,guiRuns=0,
        nativePeakGroupRSSBytes=proc['peakObservedGroupRSSBytes'],addressSpaceLimitMiB=proc['addressSpaceLimitMiB'],
        groupRSSLimitMiB=proc['groupRSSLimitMiB'],peakVirtualSizeKnown=False,
        arrayPayloadBytes=327852672,arrayPayloadReductionBytes=983558016,workspacePatchBuilt=False,
        aliasBytesDeduplicated=storage['aliasBytesDeduplicated'],
        outputDirectoryLogicalBytes=sum(x['bytes'] for x in entries+large),collectionSeconds=time.monotonic()-start,
        uniqueFilesDeleted=0,physicalQualification=False,newQualifiedAreaM2=0,
        timingScope='Nested times are not additive; source reading/editing/coordination excluded. Independent preparation/run/verifier parents recorded, no shared parent for all postprocessing. Logical bytes include hardlinks, not newly allocated disk.')
    (DEST.parent/'budget_receipt.json').write_text(json.dumps(budget,indent=2)+'\n')
    print(json.dumps(dict(files=len(entries),smallBytes=manifest['totalSmallBytes'],largeFiles=len(large),protected=29,collectionSeconds=budget['collectionSeconds'])))


if __name__=='__main__':collect()
