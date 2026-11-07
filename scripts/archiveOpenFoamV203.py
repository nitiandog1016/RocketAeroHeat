"""v203小源码/实际反算与压力拒绝证据；大资产及原生链接留本地。"""
from pathlib import Path
import hashlib,json,shutil,time,os
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v203_validation';DEST=ROOT/'docs/validation/v203/evidence'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[];local=[];aliases=[]
    for directory,dirs,names in os.walk(BASE):
        for name in list(dirs):
            if name in {'lnInclude','__pycache__'}:
                dirs.remove(name)
                if name=='lnInclude':aliases.append((Path(directory)/name).relative_to(ROOT).as_posix())
        for name in sorted(names):
            p=Path(directory)/name;size=p.stat().st_size
            selected=(p.suffix in {'.json','.txt','.log','.py','.C','.H','.hpp','.md','.png','.csv'} or p.name in {'controlDict','fvSchemes','fvSolution','physicalProperties','mixture.dict','files','options'})
            selected=selected and p.name not in {'wall_mapping.json','repair.json'} and size<=1024**2;digest=sha(p)
            if selected:
                q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(q)==digest
                rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest))
            else:local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest,uploaded=False))
    names=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_BOUNDED_THERMO_V203.md',
        'scripts/reviewOpenFoamInverse.py','scripts/checkNativeInverseOvershoot.py','scripts/recoverNativeInverseObservation.py',
        'scripts/deployBoostMathHeaders.py','scripts/buildBoundedNistInverse.py','scripts/checkBoundedNistDispatch.py',
        'scripts/runBoundedNistRocket.py','scripts/auditBoundedNistRocket.py','scripts/archiveOpenFoamV203.py','tests/testBoundedNistFailure.py']
    names += [p.relative_to(ROOT).as_posix() for d in ['native/openfoamBoundedNist','native/openfoamBoundedNistProbe'] for p in (ROOT/d).rglob('*') if p.is_file()]
    names += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v203').glob('*.md')]
    for name in names:
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.source' if name=='.gitattributes' else name);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    report=dict(schema='RAH.V203Evidence.v1',count=len(rows),bytes=sum(r['bytes'] for r in rows),files=rows,
        localLargeAssets=local,nativeBuildLinkDirectoriesNotFollowed=aliases,protectedFiles=29,protectedUnchanged=True,
        largeAssetsUploaded=False,uniqueDeleted=0,sourceFiles=names,elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v203/evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(count=len(rows),bytes=report['bytes'],localBytes=sum(r['bytes'] for r in local),protected=29,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
