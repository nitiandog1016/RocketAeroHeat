"""仅将v172终态中逐字相同的local0最终别名原子替换为硬链接，保留路径与SHA。"""

from pathlib import Path
import hashlib
import json
import os
import shutil


ROOT = Path(__file__).resolve().parents[2]
OUT = (ROOT/'outputs/v172_validation').resolve()
FOLDER = (OUT/'real01').resolve()
TARGET = ROOT/'docs/validation/v172/storage_receipt.json'
PREFLIGHT = OUT/'storage_preflight.json'
FOLDER.relative_to(OUT)
OUT.relative_to(ROOT.resolve())
assert not TARGET.exists() and not PREFLIGHT.exists(), '保留既有回执，不覆盖或重复去重'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024**2),b''):
            h.update(block)
    return h.hexdigest()


def identity(path):
    st = path.stat()
    return dict(path=str(path),bytes=st.st_size,sha256=sha(path),device=st.st_dev,
                fileID=st.st_ino,hardLinkCount=st.st_nlink,mtimeNs=st.st_mtime_ns)


def save_new(path,value):
    with path.open('x',encoding='utf-8') as stream:
        json.dump(value,stream,ensure_ascii=False,indent=2)
        stream.write('\n')


# 原生异常终态仍可去重已完整保存的相同别名；不触碰未保存异常状态。
engine_path = FOLDER/'engine_process.json'
parent_path = OUT/'native_real01_process.json'
engine = json.loads(engine_path.read_text(encoding='utf-8'))
parent = json.loads(parent_path.read_text(encoding='utf-8'))
assert engine['status']==parent['status']=='terminal'
assert engine['exitCode']==parent['exitCode']==134 and not engine['timedOut']
receipt_hashes = {str(p):sha(p) for p in (engine_path,parent_path)}
pairs = [('restart_flow_000000.csv','restart_flow.csv'),('surface_flow_000000.vtk','surface_flow.vtk')]
records=[]
for source_name,alias_name in pairs:
    source=(FOLDER/source_name).resolve()
    alias=(FOLDER/alias_name).resolve()
    assert source.parent==alias.parent==FOLDER and source!=alias
    assert not source.is_symlink() and not alias.is_symlink()
    first,second=identity(source),identity(alias)
    assert first['bytes']==second['bytes']
    assert first['sha256']==second['sha256']==engine['outputSHA256'][source_name]==engine['outputSHA256'][alias_name]
    temp=FOLDER/(alias.name+'.v172_dedup_link')
    assert not temp.exists()
    records.append(dict(source=str(source),alias=str(alias),bytes=first['bytes'],sha256=first['sha256'],
                        alreadySameFile=os.path.samefile(source,alias),sourceBefore=first,aliasBefore=second,
                        temporaryLink=str(temp)))
save_new(PREFLIGHT,dict(status='verified-before-atomic-link-replacement',pairs=records,
                       receiptSHA256=receipt_hashes,noUniqueDeletionAuthorized=True))

free_before=shutil.disk_usage(ROOT).free
for record in records:
    source=Path(record['source']);alias=Path(record['alias']);temp=Path(record['temporaryLink'])
    assert sha(source)==sha(alias)==record['sha256']
    if not record['alreadySameFile']:
        os.link(source,temp)
        assert os.path.samefile(source,temp) and sha(temp)==record['sha256']
        # 临时链接先核验；仅原子替换已经逐字证明重复的别名，不移动或删除独有场。
        os.replace(temp,alias)
    assert os.path.samefile(source,alias) and sha(source)==sha(alias)==record['sha256']
    record.update(pathPreserved=source.is_file() and alias.is_file(),nowSameFile=True,
                  sourceAfter=identity(source),aliasAfter=identity(alias),temporaryLinkRemaining=temp.exists())
free_after=shutil.disk_usage(ROOT).free
stable={name:sha(Path(name))==value for name,value in receipt_hashes.items()}
assert all(stable.values())
report=dict(passed=True,sourceTaskExitCode=134,sourceTaskTimedOut=False,lastCompleteLocalIteration=0,
            pairs=records,uniqueFilesDeleted=0,pathsAndHashesPreserved=True,
            duplicateLogicalBytesReleased=sum(r['bytes'] for r in records if not r['alreadySameFile']),
            duplicateSingleLinkBytesBefore=sum(r['bytes'] for r in records if not r['alreadySameFile'] and r['aliasBefore']['hardLinkCount']==1),
            observedFreeBytesBefore=free_before,observedFreeBytesAfter=free_after,
            observedFreeBytesDelta=free_after-free_before,
            freeSpaceMeasurementIncludesConcurrentSystemActivity=True,
            receiptSHA256=receipt_hashes,receiptsStable=stable,
            checkedFilePaths=4,allFourContentHashesPreserved=True,
            scope='v172 real01 two complete local0 files and their final aliases only',
            historyDirectoriesScanned=False,largeUniqueResultsStillLocalOnly=True)
save_new(TARGET,report)
print(json.dumps({key:value for key,value in report.items() if key not in {'pairs','receiptSHA256','receiptsStable'}},ensure_ascii=False,indent=2))
