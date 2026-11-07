"""收集v184小证据，保留真实失败及准确执行源；大场只登记，不删除独有文件。"""
from pathlib import Path
import hashlib
import json
import shutil
import time

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'outputs/v184_validation'
DEST = ROOT / 'docs/validation/v184/evidence'


def sha(path):
    """分块计算原始文件哈希，避免把大场整体载入内存。"""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read(name):
    return json.loads((BASE / name).read_text('utf-8-sig'))


def collect():
    """模块一：只读核对已有终态，不启动求解、构建或再次执行旧测试。"""
    start = time.monotonic()
    assert not DEST.exists(), '已有交付证据不覆盖'
    protected = json.loads((ROOT / 'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks = {p: sha(ROOT / p.replace('\\', '/')) == h for p, h in protected.items()}
    assert len(checks) == 29 and all(checks.values()), '并发保护文件发生变化'
    prep = read('real01_prepare/report.json')
    run = read('real01_run/report.json')
    proc = run['commands'][0]
    parent = read('parent.json')
    failure = read('failure02/report.json')
    first_failure = read('failure01/failure.json')
    closed = read('closed_admission01/report.json')
    admission = read('admission01/report.json')
    memory = read('memory01/report.json')
    storage = read('storage_receipt.json')
    assert prep['passed'] and prep['nativeGridReused'] and prep['uniformInitialFieldBitwiseIdentical']
    assert run['passed'] is False and proc['nativeExitCode'] == 1 and proc['stopReason'] is None
    assert proc['inputsStable'] and not proc['liveGroupAfterCleanup']
    assert parent['status'] == 'terminal' and parent['exitCode'] == 1 and parent['elapsedSeconds'] < 720
    assert run['outputBinding'] == 'at-supervisor-completion-audit-failed'
    assert failure['passed'] and len(failure['rejected']) == 4 and failure['historyDataRows'] == 0
    assert failure['newIterations'] == failure['newDevelopedFields'] == failure['newWallLoads'] == 0
    assert first_failure['passed'] is False and (BASE / 'failure01/executed_verifier.py').is_file()
    assert admission['passed'] and len(admission['rejected']) == 10
    assert closed['passed'] and len(closed['rejected']) == 11 and closed['legacyCompletionPassed']
    for receipt in (failure, closed):
        assert all(sha(ROOT / p) == h for p, h in receipt['sourceSHA256'].items())
    assert sha(ROOT / memory['patch']) == memory['patchSHA256']
    assert memory['status'] == 'unbuilt-allocation-only-draft' and memory['engineBuilds'] == 0
    assert storage['uniqueFilesDeleted'] == 0 and storage['allRunOutputSHA256Preserved']
    for record in storage['records']:
        path, source = ROOT / record['path'], ROOT / record['source']
        assert path.samefile(source) and sha(path) == record['sha256']

    # 模块二：仅遍历本轮目录；小证据逐字复制，大场只登记原路径与身份。
    entries, large, observed = [], [], {}
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.suffix == '.pyc':
            continue
        digest = sha(path)
        observed[str(path)] = digest
        record = dict(source=path.relative_to(ROOT).as_posix(), bytes=path.stat().st_size, sha256=digest)
        if record['bytes'] > 1024 ** 2:
            large.append(record)
            continue
        target = DEST / path.relative_to(BASE)
        if target.suffix == '.log':
            target = target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        assert sha(target) == digest
        entries.append(dict(path=target.relative_to(ROOT).as_posix(), **record))
    assert all(observed[p] == h for p, h in run['outputSHA256'].items()), '运行退出后输出身份变化'
    sources = [
        'scripts/eilmerExplicitProfile.py', 'scripts/eilmerExplicitAdmission.py',
        'scripts/prepareEilmerDevelopmentCase.py', 'scripts/runPreparedEilmerCase.py',
        'scripts/readEilmerNative.py', 'scripts/auditEilmerExplicitCompletion.py',
        'tests/testEilmerExplicitAdmission.py', 'tests/verifyEilmerExplicitFailure.py',
        'patches/eilmer/explicit-dense-storage-v1.patch',
    ]
    manifest = dict(
        files=entries, filesCount=len(entries), totalSmallBytes=sum(x['bytes'] for x in entries),
        largeLocalOnly=large, finalSourcesSHA256={p: sha(ROOT / p) for p in sources},
        protectedFiles=checks, allRunOutputSHA256Preserved=True,
        originalFailureRetained=True, uniqueFilesDeleted=0, newQualifiedAreaM2=0,
        note='One full-domain native attempt failed during initialization. Exact executed sources and the first read-only verifier failure retained. The storage patch is unbuilt; no new developed field or physical qualification.')
    (DEST.parent / 'evidence_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

    # 模块三：分别记入口、原生监督和父时间，不重复相加或当作开发总墙钟。
    budget = dict(
        preparationSeconds=prep['elapsedSeconds'], preparationBudgetSeconds=120,
        preparationPlannedParentLimitSeconds=150, preparationSeparateParentReceiptAvailable=False,
        nativeSupervisionSeconds=proc['elapsedSeconds'], nativeInternalWallClockLimitSeconds=600,
        runEntrySeconds=run['elapsedSeconds'], runEntryBudgetSeconds=660,
        runParentSeconds=parent['elapsedSeconds'], runParentLimitSeconds=parent['parentLimitSeconds'],
        admissionSeconds=admission['elapsedSeconds'], closedAdmissionSeconds=closed['elapsedSeconds'],
        failedVerifierToolSecondsIncludingStartup=first_failure['toolWallSecondsIncludingStartup'],
        recoveredVerifierSeconds=failure['elapsedSeconds'], memorySourceAuditSeconds=memory['elapsedSeconds'],
        storageDeduplicationSeconds=storage['elapsedSeconds'], postprocessingBudgetSeconds=180,
        postprocessingPlannedParentLimitSeconds=210, postprocessingEnclosingParentReceiptAvailable=False,
        nativeRuns=1, fullDomainAttempts=1, nativeUpdates=0, newDevelopedFields=0, newWallLoads=0,
        engineBuilds=0, matlabRuns=0, guiRuns=0,
        nativePeakGroupRSSBytes=proc['peakObservedGroupRSSBytes'],
        addressSpaceLimitMiB=proc['addressSpaceLimitMiB'], groupRSSLimitMiB=proc['groupRSSLimitMiB'],
        peakVirtualSizeKnown=False, failedAllocationLineKnown=False,
        nativeCaseLogicalBytes=proc['finalCaseBytes'],
        outputDirectoryLogicalBytes=sum(x['bytes'] for x in entries + large),
        largeLocalOnlyLogicalBytes=sum(x['bytes'] for x in large),
        collectionSeconds=time.monotonic() - start,
        aliasBytesDeduplicated=storage['aliasBytesDeduplicated'], uniqueFilesDeleted=0,
        proposedArrayPayloadSavingBytes=memory['potentialPayloadSavingBytes'], storagePatchBuilt=False,
        physicalQualification=False, newQualifiedAreaM2=0,
        timingScope='Nested times are not additive. Program measurements exclude editing, source reading and coordination; no single recorded parent covers all postprocessing. Logical file sizes include hardlinks and are not allocated space or newly consumed disk.')
    (DEST.parent / 'budget_receipt.json').write_text(json.dumps(budget, indent=2) + '\n')
    print(json.dumps(dict(files=len(entries), smallBytes=manifest['totalSmallBytes'], largeFiles=len(large),
                          protectedFiles=len(checks), collectionSeconds=budget['collectionSeconds'])))


if __name__ == '__main__':
    collect()
