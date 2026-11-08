"""v206小证据按字节归档；跳过原生链接目录，大场/运行库仅登记且保留E盘。"""
from pathlib import Path
import hashlib,json,os,shutil,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v206_validation';DEST=ROOT/'docs/validation/v206/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[];local=[];skipped=[]
    for directory,dirs,names in os.walk(BASE):
        for name in list(dirs):
            if name in {'lnInclude','__pycache__'}:dirs.remove(name);skipped.append((Path(directory)/name).relative_to(ROOT).as_posix())
        for name in sorted(names):
            p=Path(directory)/name;size=p.stat().st_size;digest=sha(p)
            selected=(p.suffix in {'.json','.txt','.log','.py','.m','.C','.cpp','.H','.hpp','.md','.csv','.png'} or p.name in {'files','options','mixture.dict','LICENSE','OPENFOAM_COPYING'}) and size<=1024**2
            if selected:
                q=DEST/p.relative_to(BASE)
                if q.name=='.gitattributes':q=q.with_name('gitattributes.snapshot.txt')
                q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
                if sha(q)!=digest:raise ValueError('归档字节改变')
                rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest))
            else:local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest,uploaded=False))
    files=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_PHASE_AWARE_AIR_V206.md',
        'data/thermo/phase_aware_dry_air_v1.json','phaseAwareAirState.m','phaseAwareAirReferenceState.m',
        'scripts/phaseAwareDryAir.py','scripts/preparePhaseAwareDryAir.py','scripts/checkPhaseAwareMatlab.py',
        'scripts/buildPhaseAwareAir.py','scripts/auditPhaseAwareNative.py','scripts/runPhaseAwareAirRocket.py',
        'scripts/auditPhaseAwareAirRocket.py','scripts/finalizePhaseAwareAirEvidence.py','scripts/archivePhaseAwareAirV206.py',
        'tests/testPhaseAwareAir.py','tests/testPhaseAwareAirState.m','tests/testPhaseAwareAirRocketFailure.py']
    files += [p.relative_to(ROOT).as_posix() for d in ['native/openfoamPhaseAwareAir','native/openfoamPhaseAwareAirProbe'] for p in (ROOT/d).rglob('*') if p.is_file()]
    files += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v206').glob('*.md')]
    for name in files:
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.snapshot.txt' if name=='.gitattributes' else name)
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    if len(protected)!=29 or not all(sha(ROOT/p.replace('\\','/'))==h for p,h in protected.items()):raise ValueError('并发保护项改变')
    result=dict(schema='RAH.V206Evidence.v1',count=len(rows),bytes=sum(r['bytes'] for r in rows),files=rows,localAssets=local,
        skippedLinkDirectories=skipped,protectedFiles=29,protectedUnchanged=True,sourceFiles=files,uniqueDeleted=0,
        fullPapersUploaded=False,largeLibrariesUploaded=False,CFDRuns=1,newQualifiedHeatAreaM2=0,goalComplete=False,
        elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v206/evidence_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['count','bytes','protectedFiles','elapsedSeconds']}))

if __name__=='__main__':main()
