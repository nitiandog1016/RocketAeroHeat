"""仅收集v179小证据；大场/原文保留E盘，保护并发修改，不删除任何源数据。"""
from pathlib import Path
import hashlib
import json
import shutil

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'outputs/v179_validation'
DEST=ROOT/'docs/validation/v179/evidence'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def collect():
    assert not DEST.exists()
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks={p:sha(ROOT/p)==v for p,v in protected.items()}
    assert all(checks.values()),'并发保护文件变化，不自动接管'
    entries=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.suffix=='.pyc':continue
        # 原论文和原文截图只留本地；科学结果小图随证据交付。
        if p.stat().st_size>1024**2 or p.name in ('aiaa2004_0529.pdf','aiaa2004_0529.txt','geometry-page.png'):
            local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
            continue
        target=DEST/p.relative_to(BASE)
        if target.suffix=='.log':target=target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(target)))
    report=dict(files=entries,filesCount=len(entries),totalSmallBytes=sum(x['bytes'] for x in entries),
        largeOrSourceLocalOnly=local,protectedFiles=checks,uniqueFilesDeleted=0,newQualifiedAreaM2=0,
        note='Reference tables reused by SHA from v178; no whole-project copy, no deletion, no physical qualification.')
    (DEST.parent/'evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(files=len(entries),bytes=report['totalSmallBytes'],localOnly=len(local),protected=len(checks))))


if __name__=='__main__':collect()
