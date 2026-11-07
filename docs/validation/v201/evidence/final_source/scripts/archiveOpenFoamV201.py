"""v201公开小证据；完整体场/网格/逐原面数据留E盘并在清单中绑定。"""
from pathlib import Path
import hashlib,json,shutil,time
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v201_validation';DEST=ROOT/'docs/validation/v201/evidence'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts:continue
        selected=(p.suffix in {'.json','.txt','.log','.py','.C','.H','.md','.png'} or p.name in {'controlDict','fvSchemes','fvSolution','physicalProperties'}) and p.name not in {'wall_mapping.json'}
        if selected:
            q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
            assert sha(p)==sha(q)
            rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
        else:local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p),uploaded=False))
    files=['AGENTS.md','docs/PROJECT_GOALS.md','docs/HYBRID_ROCKET_TEMPERATURE_LIMIT_V201.md',
        'scripts/reviewOpenFoamRocketBudget.py','scripts/runOpenFoamRocketDevelopment.py','scripts/readOpenFoamRocketState.py',
        'scripts/auditOpenFoamRocketDevelopment.py','scripts/plotOpenFoamRocketDiagnostic.py','scripts/checkRocketExpansionTemperatureLimit.py',
        'scripts/archiveOpenFoamV201.py','tests/testOpenFoamRocketAdmission.py','tests/testOpenFoamRocketState.py']
    files += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v201').glob('*.md')]
    for name in files:
        p=ROOT/name;q=DEST/'final_source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    report=dict(schema='RAH.V201Evidence.v1',count=len(rows),bytes=sum(r['bytes'] for r in rows),files=rows,localLargeAssets=local,
        protectedFiles=29,protectedUnchanged=True,largeAssetsUploaded=False,uniqueDeleted=0,
        outputLogicalBytes=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file()),elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v201/evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['count','bytes','outputLogicalBytes','elapsedSeconds']}))


if __name__=='__main__':main()
