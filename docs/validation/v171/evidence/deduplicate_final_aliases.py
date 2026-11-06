"""模块用途：仅把已终态且逐字相同的两个末帧别名改为硬链接，路径和内容保持。"""
from pathlib import Path
import hashlib
import json
import os
import shutil

root=Path(__file__).resolve().parents[2]
folder=(root/'outputs/v171_validation/real01').resolve()
folder.relative_to((root/'outputs/v171_validation').resolve())
receipt=json.loads((folder/'engine_process.json').read_text())
assert receipt['status']=='terminal'
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024**2),b''):h.update(block)
    return h.hexdigest()

pairs=[('restart_flow_001500.csv','restart_flow.csv'),('surface_flow_001500.vtk','surface_flow.vtk')]
records=[];free_before=shutil.disk_usage(root).free
for source_name,alias_name in pairs:
    source=(folder/source_name).resolve();alias=(folder/alias_name).resolve()
    assert source.parent==folder and alias.parent==folder and source!=alias
    digest=sha(source);assert digest==sha(alias)==receipt['outputSHA256'][alias_name]
    assert digest==receipt['outputSHA256'][source_name]
    before=alias.stat();same=os.path.samefile(source,alias)
    temp=folder/(alias.name+'.v171_dedup_link')
    assert not temp.exists()
    if not same:
        os.link(source,temp)
        assert os.path.samefile(source,temp) and sha(temp)==digest
        os.replace(temp,alias)
    assert sha(alias)==digest and os.path.samefile(source,alias)
    records.append(dict(source=str(source),alias=str(alias),sha256=digest,bytes=before.st_size,
        alreadySameFile=same,pathPreserved=True,nowSameFile=True))
free_after=shutil.disk_usage(root).free
report=dict(passed=True,pairs=records,uniqueFilesDeleted=0,pathsAndHashesPreserved=True,
    duplicateLogicalBytesReleased=sum(r['bytes'] for r in records if not r['alreadySameFile']),
    observedFreeBytesDelta=free_after-free_before,
    freeSpaceMeasurementIncludesConcurrentSystemActivity=True,largeUniqueResultsStillLocalOnly=True)
target=root/'docs/validation/v171/storage_receipt.json'
with target.open('x',encoding='utf8') as f:json.dump(report,f,indent=2);f.write('\n')
print(report)
