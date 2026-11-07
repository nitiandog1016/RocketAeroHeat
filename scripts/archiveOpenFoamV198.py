"""公开交付小证据；网格/二进制留E盘，按原字节SHA绑定而非复制工程包。"""
from pathlib import Path
import hashlib,json,shutil,subprocess,time

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v198_validation'
DEST=ROOT/'docs/validation/v198/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False)
    suffixes={'.json','.txt','.log','.py','.C','.H','.md','.dict','.list','.asc'}
    special={'controlDict','fvSchemes','fvSolution','physicalProperties','files','options'}
    rows=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        if p.suffix not in suffixes and p.name not in special:continue
        if '/Make/linux' in p.as_posix():continue
        rel=p.relative_to(BASE);q=DEST/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==sha(p)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    # 原生运行元数据虽在忽略目录，单独保存明确小清单，二进制不上传。
    runtime=ROOT/'tools/runtime/openfoam/13-20260624/runtime.json'
    (DEST/'runtime.json').write_bytes(runtime.read_bytes())
    rows.append(dict(path=(DEST/'runtime.json').relative_to(ROOT).as_posix(),source=runtime.relative_to(ROOT).as_posix(),bytes=runtime.stat().st_size,sha256=sha(runtime)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    report=dict(schema='RAH.V198Evidence.v1',files=rows,count=len(rows),bytes=sum(x['bytes'] for x in rows),
        originalFilesUnmodified=True,protectedFiles=29,largeAssetsUploaded=False,uniqueDeleted=0,
        elapsedSeconds=time.monotonic()-start)
    manifest=ROOT/'docs/validation/v198/evidence_manifest.json'
    manifest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['count','bytes','protectedFiles','elapsedSeconds']}))

if __name__=='__main__':main()
