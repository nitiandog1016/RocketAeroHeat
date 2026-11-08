"""v204按字节归档小源码/证据；论文全文、执行库与无关并发改动留在本地。"""
from pathlib import Path
import hashlib,json,os,shutil,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v204_validation';DEST=ROOT/'docs/validation/v204/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[];local=[];skipped=[]
    for directory,dirs,names in os.walk(BASE):
        for name in list(dirs):
            if name in {'lnInclude','__pycache__'}:dirs.remove(name);skipped.append((Path(directory)/name).relative_to(ROOT).as_posix())
        for name in sorted(names):
            p=Path(directory)/name;size=p.stat().st_size;digest=sha(p)
            selected=(p.suffix in {'.json','.txt','.log','.py','.C','.cpp','.H','.hpp','.bib','.md','.csv','.png'} or p.name in {'files','options','mixture.dict','LICENSE','OPENFOAM_COPYING'}) and size<=1024**2
            selected=selected and not name.startswith(('lemmon','invalid_'))
            if selected:
                q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(q)==digest
                rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest))
            else:local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest,uploaded=False))
    files=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_PRESSURE_DOMAIN_V204.md',
        'data/thermo/standard_dry_air_cp_v1.json','scripts/reviewNistPressureDomain.py','scripts/validateNistPressureDomain.py',
        'scripts/runNistPressureReference.py','scripts/standardDryAir.py','scripts/prepareStandardDryAirDomain.py',
        'scripts/buildStandardDryAirDomain.py','scripts/auditStandardNativeDomain.py','scripts/plotStandardDryAirDomain.py',
        'scripts/finalizeStandardDryAirEvidence.py','scripts/archiveStandardDryAirV204.py','tests/testStandardDryAirDomain.py']
    files += [p.relative_to(ROOT).as_posix() for d in ['native/openfoamStandardDryAir','native/openfoamStandardDryAirProbe'] for p in (ROOT/d).rglob('*') if p.is_file()]
    files += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v204').glob('*.md')]
    for name in files:
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.source' if name=='.gitattributes' else name);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'));assert len(protected)==29
    assert all(sha(ROOT/p)==h for p,h in protected.items())
    result=dict(schema='RAH.V204Evidence.v1',count=len(rows),bytes=sum(r['bytes'] for r in rows),files=rows,localAssets=local,
        skippedLinkDirectories=skipped,protectedFiles=29,protectedUnchanged=True,sourceFiles=files,uniqueDeleted=0,
        paperFullTextUploaded=False,largeLibrariesUploaded=False,CFDRuns=0,goalComplete=False,elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v204/evidence_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['count','bytes','protectedFiles','elapsedSeconds']}))

if __name__=='__main__':main()
