"""归档v209小证据、明确未实现项；完整第三方样例和网格留在E盘。"""
from pathlib import Path
import hashlib,json,shutil,time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v209_validation'
DEST=ROOT/'docs/validation/v209/evidence'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[];local=[]
    # 模块一：仅发布审核摘录和执行/失败回执，不重复打包完整上游源码。
    for folder in ['source01','source02','source03','source04','audit01']:
        for p in sorted((BASE/folder).rglob('*')):
            if not p.is_file():continue
            rel=p.relative_to(BASE)
            selected=(folder=='audit01' or (p.parent==BASE/folder and p.suffix in ['.py','.json']))
            if selected:
                q=DEST/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
                assert sha(q)==sha(p)
                rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
            else:local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p),uploaded=False))
    license_source=BASE/'source02/source/COPYING'
    q=DEST/'GDTK_COPYING';shutil.copyfile(license_source,q)
    rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=license_source.relative_to(ROOT).as_posix(),bytes=q.stat().st_size,sha256=sha(q)))
    # 模块二：保存交付时的小源码快照，历史执行版独立保存。
    sources=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md',
        'docs/HYBRID_ANISOTROPIC_ROUTE_DECISION_V209.md','scripts/reviewAnisotropicViscousEntry.py',
        'scripts/auditAnisotropicViscousEntry.py','scripts/archiveAnisotropicRouteV209.py']
    sources += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v209').iterdir() if p.is_file() and p.suffix in ['.md','.json']]
    for name in sources:
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.source' if name=='.gitattributes' else name)
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==sha(p)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    # 模块三：检查准确状态及并发保护，不把源码摘录当作运行结果。
    r=json.loads((BASE/'audit01/report.json').read_text())
    c=json.loads((ROOT/'docs/validation/v209/route_contract.json').read_text())
    assert r['passed'] and not r['existingCombinationAccepted']
    assert not c['executionAuthorized'] and not c['physicalQualification']
    assert not c['nextImplementation']['nativeHelperBuildAuthorized']
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    assert sum(row['bytes'] for row in rows)<8*1024**2
    manifest=dict(schema='RAH.V209Evidence.v1',files=rows,count=len(rows),bytes=sum(r['bytes'] for r in rows),
        localAssets=local,sourceFiles=sources,protectedFiles=29,protectedUnchanged=True,
        uniqueDeleted=0,builds=0,CFDRuns=0,physicalQualification=False,goalComplete=False,
        elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v209/evidence_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:manifest[k] for k in ['count','bytes','protectedFiles','elapsedSeconds']}))

if __name__=='__main__':main()
