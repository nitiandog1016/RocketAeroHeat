"""按原字节归档v210小证据；二进制、对象文件与完整上游树留在E盘。"""
from pathlib import Path
import hashlib,json,shutil,time,ast
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v210_validation';DEST=ROOT/'docs/validation/v210/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        selected=p.suffix in ['.json','.txt','.py','.d','.lua','.tsv'] and p.stat().st_size<=1024**2
        if selected:
            q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
            rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
        else:local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p),uploaded=False))
    names=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_PHASE_AWARE_EILMER_V210.md',
        'scripts/prepareEilmerPhaseAwareAir.py','scripts/prepareEilmerQueryRecovery.py','scripts/runEilmerPhaseAirHelper.py',
        'scripts/verifyEilmerPhaseAirHelper.py','scripts/finalizeEilmerPhaseAirV210.py','scripts/archiveEilmerPhaseAirV210.py']
    names += [p.relative_to(ROOT).as_posix() for p in (ROOT/'native/eilmerPhaseAwareAir').rglob('*') if p.is_file()]
    names += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v210').glob('*.md')]
    for name in names:
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.source' if name=='.gitattributes' else name)
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    for name in names:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf-8'))
    for name in ['scripts/runEilmerPhaseAirHelper.py','scripts/prepareEilmerPhaseAwareAir.py','scripts/prepareEilmerQueryRecovery.py']:
        tree=ast.parse((ROOT/name).read_text(encoding='utf-8'));assert isinstance(tree.body[-1].body[0],ast.Raise)
    assert sum(r['bytes'] for r in rows)<8*1024**2
    manifest=dict(schema='RAH.V210Evidence.v1',files=rows,count=len(rows),bytes=sum(r['bytes'] for r in rows),
        localAssets=local,sourceFiles=names,protectedFiles=29,protectedUnchanged=True,uniqueDeleted=0,
        engineBuilds=0,CFDRuns=0,physicalQualification=False,goalComplete=False,elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v210/evidence_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:manifest[k] for k in ['count','bytes','protectedFiles','elapsedSeconds']}))
if __name__=='__main__':main()
