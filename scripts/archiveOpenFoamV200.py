"""v200小证据公开交付；大体场/网格和完整采样映射留E盘并以SHA绑定。"""
from pathlib import Path
import hashlib,json,shutil,time

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v200_validation';DEST=ROOT/'docs/validation/v200/evidence'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[]
    suffixes={'.json','.txt','.log','.py','.C','.H','.md','.csv','.png'}
    special={'controlDict','fvSchemes','fvSolution','physicalProperties'}
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        if p.name=='wall_mapping.json':continue
        if p.suffix not in suffixes and p.name not in special:continue
        q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(q)==sha(p)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    files=['AGENTS.md','docs/PROJECT_GOALS.md','docs/HYBRID_ROCKET_OPENFOAM_READY_V200.md',
        'scripts/reviewOpenFoamRestart.py','scripts/runOpenFoamRestart.py','scripts/auditOpenFoamRestart.py',
        'scripts/prepareOpenFoamRocket.py','scripts/checkOpenFoamRocket.py','scripts/auditOpenFoamRocketGeometry.py',
        'scripts/relocateOpenFoamRocket.py','scripts/prepareOpenFoamControl.py','scripts/archiveOpenFoamV200.py',
        'tests/verifyOpenFoamRelocation.py','tests/testOpenFoamV200Contracts.py']
    files += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v200').glob('*.md')]
    for name in files:
        p=ROOT/name;q=DEST/'final_source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    report=dict(schema='RAH.V200Evidence.v1',count=len(rows),bytes=sum(r['bytes'] for r in rows),files=rows,
        protectedFiles=29,protectedUnchanged=True,largeAssetsUploaded=False,uniqueDeleted=0,
        outputLogicalBytes=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file()),elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v200/evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['count','bytes','outputLogicalBytes','elapsedSeconds']}))


if __name__=='__main__':main()
