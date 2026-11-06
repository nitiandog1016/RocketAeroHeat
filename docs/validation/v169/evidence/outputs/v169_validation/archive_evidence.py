"""只归档本轮小证据；不复制运行包/原PDF/历史大场，不删除任何源文件。"""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/v169_validation'
DEST=ROOT/'docs/validation/v169/evidence'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
selected=[Path(__file__).resolve(),OUT/'correct_recovery_pick.py',OUT/'inlet_contract.json']
selected+=list(OUT.glob('*_process.json'))+list(OUT.glob('*.log'))
dirs=['axis01','axis02','native_source','native02_source','native03_source','recovery01','recovery02']
dirs+=['precheck%02d'%i for i in range(1,4)]
dirs+=['audit3d%02d'%i for i in range(1,6)]
dirs+=['contracts%02d'%i for i in range(1,5)]
for name in dirs:
    selected += [p for p in (OUT/name).rglob('*') if p.is_file() and p.suffix in
                 ('.py','.cpp','.hpp','.inl','.json','.csv','.md','.log')]
aliases=[]
for name in ['small3d01','small3d02','small3d03']:
    folder=OUT/name
    for p in folder.iterdir():
        if not p.is_file():continue
        if p.name in ('restart_flow.csv','surface_flow.vtk'):
            target=folder/(p.stem+'_000002'+p.suffix)
            assert sha(p)==sha(target)
            aliases.append(dict(source=p.relative_to(ROOT).as_posix(),sameAs=target.relative_to(ROOT).as_posix(),sha256=sha(p)))
        elif p.suffix in ('.su2','.cfg','.json','.csv','.vtk','.log'):
            selected.append(p)
selected=sorted(set(selected))
total=sum(p.stat().st_size for p in selected)
assert total<16*1024**2 and not DEST.exists()
entries=[]
for source in selected:
    target=DEST/source.relative_to(ROOT)
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(str(source),str(target))
    assert sha(target)==sha(source)
    entries.append(dict(source=source.relative_to(ROOT).as_posix(),path=target.relative_to(ROOT).as_posix(),
                        bytes=target.stat().st_size,sha256=sha(target)))
manifest=dict(files=entries,totalFiles=len(entries),totalBytes=total,
    unarchivedByteIdenticalFinalAliases=aliases,sourceDeleted=False,largeFieldsCopied=False,
    CFDExecutableCopied=False,experimentalQualification=False)
(ROOT/'docs/validation/v169/evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')
print(json.dumps(dict(files=len(entries),bytes=total,aliases=len(aliases))))
