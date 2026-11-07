"""仅收集v181小输入/失败/准确源码；不复制工程、删除独有数据或触碰并发文件。"""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[3];BASE=ROOT/'outputs/v181_validation';DEST=ROOT/'docs/validation/v181/evidence'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def collect():
    assert not DEST.exists()
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text('utf-8-sig'))
    checks={p:sha(ROOT/Path(p.replace('\\','/')))==h for p,h in protected.items()};assert all(checks.values())
    entries=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.suffix=='.pyc':continue
        if p.stat().st_size>1024**2:
            local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)));continue
        target=DEST/p.relative_to(BASE)
        if target.suffix=='.log':target=target.with_suffix('.log.raw')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        entries.append(dict(path=target.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(target)))
    names=['scripts/auditEilmerCompletion.py','scripts/prepareEilmerDevelopmentCase.py','scripts/prepareEilmerRocketCase.py',
           'scripts/readEilmerNative.py','scripts/runPreparedEilmerCase.py','tests/testEilmerRealAdmission.py',
           'tests/testEilmerRealProfile.py','tests/verifyEilmerRealRun.py']
    r=dict(files=entries,filesCount=len(entries),totalSmallBytes=sum(x['bytes'] for x in entries),largeLocalOnly=local,
           finalSourcesSHA256={p:sha(ROOT/p) for p in names},protectedFiles=checks,uniqueFilesDeleted=0,newQualifiedAreaM2=0,
           note='Actual native comparison failed; scope/identity checks are separate. Original run sources retained. No full-domain CFD or build.')
    (DEST.parent/'evidence_manifest.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(files=len(entries),bytes=r['totalSmallBytes'],protected=len(checks),largeLocalOnly=len(local))))
if __name__=='__main__':collect()
