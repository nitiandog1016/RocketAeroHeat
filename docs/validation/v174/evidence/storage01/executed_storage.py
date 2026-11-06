"""只将v174三个终态任务的六个逐字相同最终别名替为同目录硬链接。"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v174_validation'
EVIDENCE=BASE/'storage01'
TARGET=ROOT/'docs/validation/v174/storage_receipt.json'
TASKS=[('small02',3,'review_small02'),('small03',3,'audit_small03'),('real01',100,'audit_real01')]


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def save(path,value):
    temporary=path.with_name(path.name+'.pending')
    with temporary.open('w',encoding='utf8') as stream:
        json.dump(value,stream,ensure_ascii=False,indent=2,allow_nan=False)
        stream.write('\n');stream.flush();os.fsync(stream.fileno())
    os.replace(str(temporary),str(path))


def checked_path(path,exists=True):
    """变更前检查词法绝对路径及每级祖先，拒绝symlink/junction/reparse。"""
    path=Path(path);assert path.is_absolute()
    assert path==Path(os.path.abspath(str(path))) and path.resolve()==path
    path.relative_to(BASE)
    candidate=path if exists else path.parent
    while True:
        value=os.lstat(str(candidate))
        assert not candidate.is_symlink() and not (getattr(value,'st_file_attributes',0)&0x400), '拒绝reparse路径'
        if candidate==ROOT:break
        assert candidate!=candidate.parent;candidate=candidate.parent
    if exists:assert path.exists()
    return path


def info(path,digest=None):
    value=path.stat()
    return dict(path=str(path),bytes=value.st_size,sha256=digest if digest else sha(path),
        device=value.st_dev,fileID=value.st_ino,hardLinkCount=value.st_nlink,
        mtimeNs=value.st_mtime_ns,ctimeNs=value.st_ctime_ns)


def compare_and_hash(source,alias):
    """一次分块读取同时直接比较和计算两原始SHA，不解析场或猜近似重复。"""
    first=hashlib.sha256();second=hashlib.sha256()
    with source.open('rb') as left,alias.open('rb') as right:
        while True:
            a=left.read(1024**2);b=right.read(1024**2)
            assert a==b,'拒绝不同字节的别名：'+str(alias)
            if not a:break
            first.update(a);second.update(b)
    assert first.hexdigest()==second.hexdigest()
    return first.hexdigest()


def same_original_identity(path,before):
    current=info(path,before['sha256'])
    return all(current[key]==before[key] for key in ('bytes','device','fileID','mtimeNs','ctimeNs'))


def main():
    start=time.monotonic();assert BASE==Path('E:/RocketAeroHeat/outputs/v174_validation')
    checked_path(BASE);assert not EVIDENCE.exists() and not TARGET.exists()
    EVIDENCE.mkdir();(EVIDENCE/'executed_storage.py').write_bytes(Path(__file__).read_bytes())
    result=dict(passed=False,scope='only v174 small02 local3 / small03 local3 / real01 local100 final aliases',
        absoluteBoundary=str(BASE),tasks=[],pairs=[],uniqueFilesDeleted=0,runtimeModified=False,
        otherVersionsModified=False,OpenCodeModified=False,historyDirectoriesScanned=False,
        largeUniqueResultsStillLocalOnly=True,largeResultsBackedUpToGitHub=False,otherEvidenceArchived=False,
        safetyPatternSource='scripts/deduplicateSU2V173Copies.py',scriptSHA256=sha(Path(__file__)))
    before={};receipt_locks={};plan=[]
    try:
        for label,step,review_name in TASKS:
            folder=checked_path(BASE/label);assert folder.parent==BASE
            receipt_path=folder/'engine_process.json';receipt=read(receipt_path)
            assert receipt['status']=='terminal' and all(receipt['inputsStable'].values())
            proof_path=BASE/review_name/'report.json';proof=read(proof_path)
            assert proof['reviewCompleted'] and all(proof['inputsStable'].values())
            if label=='small02':
                assert proof['researchReviewPassed'] and not proof['normalTaskPassed']
                assert proof['originalSupervisorErrorPreserved'] and not proof['originalReceiptModified']
                assert proof['originalMissingFirstFrameStillUnresolved'] and receipt.get('supervisorError')
                assert proof['completeSavedLocalIterations']==[2,3]
            else:
                assert proof['interfacePassed'] and proof['savedFrameInterfacePassed']
                assert not receipt.get('supervisorError')
                if label=='small03':assert proof['normalRunCompleted'] and receipt['exitCode']==0
                else:
                    independent_path=BASE/'independent_real01/report.json';independent=read(independent_path)
                    assert independent['passed'] and all(independent['inputsStableMetadata'].values())
                    assert receipt['modelDomainExceeded'] and receipt['domainGuardTerminationRequested']
                    assert not receipt['cancelled'] and not proof['normalRunCompleted']
                    receipt_locks[str(independent_path)]=sha(independent_path)
            final=receipt['checkpoints']['lastComplete'];assert final['localIteration']==step
            assert step in [row['localIteration'] for row in proof['frames']]
            result['tasks'].append(dict(folder=str(folder),nativeExitCode=receipt['exitCode'],
                originalSupervisorError=receipt.get('supervisorError'),
                originalFailureClassificationUnchanged=True,lastCompleteLocalIteration=step,
                modelDomainExceeded=receipt.get('modelDomainExceeded',False),
                readOnlyQualification=review_name,physicalQualified=False))
            for path in (receipt_path,proof_path,BASE/('native_'+label+'_process.json')):
                receipt_locks[str(path)]=sha(path)
            for stem,ext,key in [('restart_flow','.csv','volume'),('surface_flow','.vtk','surface')]:
                source=checked_path(folder/(stem+'_%06d'%step+ext));alias=checked_path(folder/(stem+ext))
                assert source.parent==alias.parent==folder and source!=alias
                assert source.is_file() and alias.is_file()
                assert source.stat().st_dev==alias.stat().st_dev==folder.stat().st_dev
                digest=compare_and_hash(source,alias)
                assert digest==receipt['outputSHA256'][source.name]==receipt['outputSHA256'][alias.name]
                assert source.name==final[key] and digest==final[key+'SHA256']
                pending=checked_path(alias.with_name(alias.name+'.v174_dedup_link'),False)
                assert not pending.exists() and not pending.is_symlink()
                row=dict(source=str(source),alias=str(alias),sourceBefore=info(source,digest),
                    aliasBefore=info(alias,digest),alreadySameFile=os.path.samefile(source,alias),
                    directByteComparisonPassed=True,temporaryLink=str(pending))
                plan.append(row);before[str(source)]=digest;before[str(alias)]=digest
        assert len(plan)==6 and len(before)==12
        # 必须落盘完整preflight，之后才允许替换任何原有别名。
        preflight=dict(absoluteBoundary=str(BASE),allPathsWithinBoundary=True,noSymlinkOrReparse=True,
            allPairsSameDirectoryAndVolume=True,tasks=result['tasks'],pairs=plan,
            originalSHA256=before,receiptSHA256=receipt_locks,scriptSHA256=result['scriptSHA256'])
        save(EVIDENCE/'preflight.json',preflight)
        free_before=shutil.disk_usage(ROOT).free
        for row in plan:
            assert time.monotonic()-start<25,'内部操作截止，保留已完成回执，不扩预算'
            source=checked_path(row['source']);alias=checked_path(row['alias'])
            pending=checked_path(row['temporaryLink'],False)
            assert same_original_identity(source,row['sourceBefore']) and same_original_identity(alias,row['aliasBefore'])
            assert not pending.exists()
            if not row['alreadySameFile']:
                os.link(str(source),str(pending))
                checked_path(pending)
                assert os.path.samefile(source,pending) and pending.stat().st_dev==alias.stat().st_dev
                # 此处原alias仍在；只有已验证的同内容硬链接已准备好才原子替换。
                os.replace(str(pending),str(alias))
            assert os.path.samefile(source,alias) and not pending.exists()
            digest=compare_and_hash(source,alias)
            assert digest==row['sourceBefore']['sha256']==row['aliasBefore']['sha256']
            row.update(sourceAfter=info(source,digest),aliasAfter=info(alias,digest),
                pathPreserved=True,nowSameFile=True,temporaryLinkRemaining=False)
            result['pairs'].append(row)
            save(EVIDENCE/'operation_journal.json',dict(completedPairs=result['pairs'],originalSHA256=before))
        free_after=shutil.disk_usage(ROOT).free
        stable={path:sha(Path(path))==digest for path,digest in receipt_locks.items()}
        preserved={path:Path(path).is_file() and sha(Path(path))==digest for path,digest in before.items()}
        assert all(stable.values()) and all(preserved.values())
        result.update(passed=True,allMutationPathsCheckedWithinBoundary=True,noSymlinkOrReparse=True,
            allPairsSameDirectoryAndVolume=True,checkedFilePaths=12,originalSHA256=before,
            pathsAndHashesPreserved=all(preserved.values()),pathPreservation=preserved,
            duplicateLogicalBytesReleased=sum(row['aliasBefore']['bytes'] for row in plan if not row['alreadySameFile']),
            duplicateSingleLinkBytesBefore=sum(row['aliasBefore']['bytes'] for row in plan
                if not row['alreadySameFile'] and row['aliasBefore']['hardLinkCount']==1),
            observedFreeBytesBefore=free_before,observedFreeBytesAfter=free_after,
            observedFreeBytesDelta=free_after-free_before,
            freeSpaceMeasurementIncludesConcurrentSystemActivity=True,physicalDiskUsageClaimed=False,
            receiptSHA256=receipt_locks,receiptsStable=stable)
    except Exception as error:
        result['error']=type(error).__name__+': '+str(error);raise
    finally:
        result['elapsedSeconds']=time.monotonic()-start
        save(EVIDENCE/'report.json',result);save(TARGET,result)
        print(json.dumps({k:result.get(k) for k in ('passed','duplicateLogicalBytesReleased','observedFreeBytesDelta',
            'checkedFilePaths','elapsedSeconds','error')},ensure_ascii=False))


if __name__=='__main__':main()
