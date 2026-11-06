"""模块用途：仅对v173两个终态任务的同字节完整帧别名做硬链接去重，保留原路径和全部独有结果。"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

ROOT=Path(__file__).resolve().parents[1]
BASE=(ROOT/'outputs/v173_validation').resolve()
TARGET=ROOT/'docs/validation/v173/storage_receipt.json'
EVIDENCE=BASE/'storage01'


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024**2),b''):digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text('utf8'))


def save(path,value):
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')


def file_info(path):
    info=path.stat()
    return dict(path=str(path),bytes=info.st_size,sha256=sha(path),device=info.st_dev,
        fileID=info.st_ino,hardLinkCount=info.st_nlink,mtimeNs=info.st_mtime_ns)


def equal_bytes(source,alias):
    """模块一：SHA之外再分块直接比较，不解析或重算任何场数据。"""
    with source.open('rb') as first,alias.open('rb') as second:
        while True:
            a=first.read(1024**2);b=second.read(1024**2)
            if a!=b:return False
            if not a:return True


def hardlink_identical(source,alias):
    """模块二：复用v167/v148的先建硬链接再原子替换顺序，操作前限定绝对路径。"""
    source=source.resolve();alias=alias.resolve()
    source.relative_to(BASE);alias.relative_to(BASE)
    assert source.parent==alias.parent and source.parent.name in ('real01','small01')
    assert not source.is_symlink() and not alias.is_symlink() and source!=alias
    assert source.stat().st_dev==alias.stat().st_dev
    assert sha(source)==sha(alias) and equal_bytes(source,alias),'拒绝不同内容的别名'
    row=dict(source=str(source),alias=str(alias),sourceBefore=file_info(source),aliasBefore=file_info(alias),
        alreadySameFile=os.path.samefile(source,alias),directByteComparisonPassed=True)
    pending=alias.with_name(alias.name+'.v173_dedup_link')
    pending.resolve().relative_to(BASE)
    assert pending.parent==alias.parent and not pending.exists()
    if not row['alreadySameFile']:
        os.link(str(source),str(pending))
        os.replace(str(pending),str(alias))
    assert os.path.samefile(source,alias) and not pending.exists()
    row.update(sourceAfter=file_info(source),aliasAfter=file_info(alias),pathPreserved=True,
        nowSameFile=True,temporaryLink=str(pending),temporaryLinkRemaining=False)
    assert row['sourceAfter']['sha256']==row['aliasAfter']['sha256']==row['sourceBefore']['sha256']==row['aliasBefore']['sha256']
    return row


def main():
    """模块三：先确认两个终态及全部四对内容，再执行；只记录去重结果，不归档其他证据。"""
    started=time.monotonic()
    assert BASE==Path('E:/RocketAeroHeat/outputs/v173_validation').resolve()
    assert not TARGET.exists() and not EVIDENCE.exists()
    EVIDENCE.mkdir();(EVIDENCE/'executed_storage.py').write_bytes(Path(__file__).read_bytes())
    plan=[];receipt_locks={};before={};task_identity=[]
    for label,step,review in [('real01',0,'verify_real01'),('small01',7,'verify_small02')]:
        folder=(BASE/label).resolve();folder.relative_to(BASE)
        assert folder.parent==BASE and not folder.is_symlink()
        receipt_path=folder/'engine_process.json';receipt=read(receipt_path)
        assert receipt['status']=='terminal' and all(receipt['inputsStable'].values())
        review_path=BASE/review/'report.json';review_report=read(review_path)
        assert review_report['reviewCompleted'] and review_report['passed']
        inventory_path=BASE/review/'read_only_checkpoint_inventory.json';inventory=read(inventory_path)
        if 'checkpoints' in receipt:assert inventory==receipt['checkpoints']
        assert inventory['lastComplete']['localIteration']==step
        selected=next(row for row in inventory['complete'] if row['localIteration']==step)
        task_identity.append(dict(folder=str(folder),exitCode=receipt['exitCode'],timedOut=receipt.get('timedOut',False),
            lastCompleteLocalIteration=step,receiptNativeCheckpointsPresent='checkpoints' in receipt))
        for path in (receipt_path,review_path,inventory_path,BASE/('native_'+label+'_process.json')):
            receipt_locks[str(path)]=sha(path)
        for stem,extension,key in [('restart_flow','.csv','volume'),('surface_flow','.vtk','surface')]:
            source=(folder/(stem+'_%06d'%step+extension)).resolve()
            alias=(folder/(stem+extension)).resolve()
            for path in (source,alias):
                path.relative_to(BASE)
                assert path.parent==folder and not path.is_symlink()
                before[str(path)]=sha(path)
                assert before[str(path)]==receipt['outputSHA256'][path.name]
            assert source.name==selected[key] and before[str(source)]==selected[key+'SHA256']
            assert before[str(source)]==before[str(alias)] and equal_bytes(source,alias)
            plan.append(dict(source=str(source),alias=str(alias),bytes=alias.stat().st_size,sha256=before[str(alias)]))
    assert len(plan)==4 and len(before)==8
    save(EVIDENCE/'preflight.json',dict(absoluteBoundary=str(BASE),allPathsWithinBoundary=True,
        tasks=task_identity,pairs=plan,originalSHA256=before,receiptSHA256=receipt_locks))
    free_before=shutil.disk_usage(ROOT).free
    pairs=[hardlink_identical(Path(row['source']),Path(row['alias'])) for row in plan]
    free_after=shutil.disk_usage(ROOT).free
    preserved={path:Path(path).is_file() and sha(Path(path))==digest for path,digest in before.items()}
    stable={path:sha(Path(path))==digest for path,digest in receipt_locks.items()}
    assert all(preserved.values()) and all(stable.values())
    released=sum(row['aliasBefore']['bytes'] for row in pairs if not row['alreadySameFile'])
    single=sum(row['aliasBefore']['bytes'] for row in pairs if not row['alreadySameFile'] and row['aliasBefore']['hardLinkCount']==1)
    result=dict(passed=True,scope='v173 real01 local0 / small01 local7 named complete frames and final aliases only',
        absoluteBoundary=str(BASE),allMutationPathsCheckedWithinBoundary=True,tasks=task_identity,pairs=pairs,
        checkedFilePaths=8,originalSHA256=before,pathsAndHashesPreserved=all(preserved.values()),
        uniqueFilesDeleted=0,duplicateLogicalBytesReleased=released,duplicateSingleLinkBytesBefore=single,
        observedFreeBytesBefore=free_before,observedFreeBytesAfter=free_after,observedFreeBytesDelta=free_after-free_before,
        freeSpaceMeasurementIncludesConcurrentSystemActivity=True,physicalDiskUsageClaimed=False,
        receiptSHA256=receipt_locks,receiptsStable=stable,runtimeModified=False,otherVersionsModified=False,
        OpenCodeModified=False,historyDirectoriesScanned=False,largeUniqueResultsStillLocalOnly=True,
        largeResultsBackedUpToGitHub=False,otherEvidenceArchived=False,scriptSHA256=sha(Path(__file__)),
        safetyPatternSource=['scripts/deduplicateSU2V167Copies.py','scripts/deduplicateSU2V148Copies.py'],
        elapsedSeconds=time.monotonic()-started)
    save(TARGET,result)
    print(json.dumps({key:result[key] for key in ('passed','duplicateLogicalBytesReleased','observedFreeBytesDelta','checkedFilePaths','elapsedSeconds')}))


if __name__=='__main__':main()
