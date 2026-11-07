"""终态后仅将两份逐字相同的初场别名转为硬链接，保留原路径/SHA，不删除独有数据。"""
from pathlib import Path
import hashlib,json,os,time
ROOT=Path(__file__).resolve().parents[3].resolve();base=Path(__file__).resolve().parent
assert json.loads((base/'run_parent.json').read_text())['status']=='terminal'
r=json.loads((base/'real01_run/report.json').read_text());assert not r['commands'][0]['liveGroupAfterCleanup']
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024**2),b''):h.update(block)
    return h.hexdigest()
start=time.monotonic();records=[]
for rel in ['lmrsim/snapshots/0000/fluid-0000','lmrsim/snapshots/0000/grid-0000']:
    source=(ROOT/'outputs/v178_validation/real01'/rel).resolve();target=(base/'real01'/rel).resolve()
    source.relative_to(ROOT);target.relative_to((base/'real01').resolve())
    expected=r['outputSHA256'][str(target)]
    assert sha(source)==sha(target)==expected
    with source.open('rb') as a,target.open('rb') as b:
        while True:
            x,y=a.read(1024**2),b.read(1024**2);assert x==y
            if not x:break
    size=target.stat().st_size;already=os.path.samefile(source,target)
    if not already:
        temp=target.with_name(target.name+'.rah-identical-link');assert not temp.exists()
        temp.resolve().relative_to((base/'real01').resolve())
        os.link(source,temp);os.replace(temp,target)
    assert os.path.samefile(source,target) and sha(target)==expected
    records.append(dict(path=str(target.relative_to(ROOT)),source=str(source.relative_to(ROOT)),sha256=expected,
                        bytes=size,alreadySameFile=already,aliasBytesDeduplicated=0 if already else size))
assert all(sha(Path(p))==h for p,h in r['outputSHA256'].items())
report=dict(records=records,aliasBytesDeduplicated=sum(x['aliasBytesDeduplicated'] for x in records),
            allRunOutputSHA256Preserved=True,uniqueFilesDeleted=0,elapsedSeconds=time.monotonic()-start,
            note='Only byte-identical aliases replaced by hardlinks; no filesystem allocated-free-space claim.')
(base/'storage_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
