"""v202小源码证据归档；大场/库/完整长日志本地保留，原生链接不跟随。"""
from pathlib import Path
import hashlib,json,shutil,time,os
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v202_validation';DEST=ROOT/'docs/validation/v202/evidence'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    start=time.monotonic();DEST.mkdir(exist_ok=False);rows=[];local=[];aliases=[]
    # 模块一：OpenFOAM的lnInclude由WSL生成；只记录该链接目录，不在Windows追踪它。
    for directory,dirs,names in os.walk(BASE):
        for name in list(dirs):
            if name in {'lnInclude','__pycache__'}:
                dirs.remove(name)
                if name=='lnInclude':aliases.append((Path(directory)/name).relative_to(ROOT).as_posix())
        for name in sorted(names):
            p=Path(directory)/name;size=p.stat().st_size
            selected=(p.suffix in {'.json','.txt','.log','.py','.C','.H','.md','.png','.csv'} or p.name in {'controlDict','fvSchemes','fvSolution','physicalProperties','mixture.dict','files','options'})
            selected=selected and p.name not in {'wall_mapping.json','repair.json'} and size<=1024**2
            digest=sha(p)
            if selected:
                q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(q)==digest
                rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest))
            else:
                local.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=size,sha256=digest,uploaded=False))
                if name=='stdout.txt':
                    text=p.read_text(encoding='utf-8',errors='replace');lines=text.splitlines()
                    q=DEST/p.relative_to(BASE).parent/'stdout_excerpt.txt';q.parent.mkdir(parents=True,exist_ok=True)
                    q.write_text('Full log retained locally; SHA256 '+digest+'\nLines '+str(len(lines))+'\n\nFIRST 55 LINES\n'+'\n'.join(lines[:55])+'\n\nLAST 55 LINES\n'+'\n'.join(lines[-55:])+'\n',encoding='utf-8')
                    rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),kind='explicit-excerpt',bytes=q.stat().st_size,sha256=sha(q)))
    # 模块二：交付完整源码，不把本轮并发UI/其他研究文件混入提交。
    names=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_NIST_AIR_V202.md',
        'data/thermo/nist_air_shomate_v1.json',
        'scripts/nistIntegratedAir.py','scripts/reviewOpenFoamTabulatedThermo.py','scripts/buildOpenFoamNistAir.py',
        'scripts/openFoamFields.py','scripts/openFoamConservation.py','scripts/openFoamSynchronousControl.py','scripts/hybridEulerReference.py',
        'scripts/runOpenFoamNistControl.py','scripts/auditOpenFoamNistControl.py','scripts/runOpenFoamNistRocket.py',
        'scripts/readOpenFoamNistRocketState.py','scripts/auditOpenFoamNistRocketFailure.py','scripts/archiveOpenFoamV202.py',
        'tests/verifyNistIntegratedAir.py','tests/verifyNistCallbacks.py','tests/testNistRocketContracts.py']
    names += [p.relative_to(ROOT).as_posix() for d in ['native/openfoamNistAir','native/openfoamNistProbe'] for p in (ROOT/d).rglob('*') if p.is_file()]
    names += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v202').glob('*.md')]
    for name in names:
        # .gitattributes快照改用普通文件名，避免其规则在证据子目录中再次生效。
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.source' if name=='.gitattributes' else name);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,bytes=p.stat().st_size,sha256=sha(p)))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    report=dict(schema='RAH.V202Evidence.v1',count=len(rows),bytes=sum(r['bytes'] for r in rows),files=rows,
        localLargeAssets=local,nativeBuildLinkDirectoriesNotFollowed=aliases,protectedFiles=29,protectedUnchanged=True,
        largeAssetsUploaded=False,uniqueDeleted=0,sourceFiles=names,elapsedSeconds=time.monotonic()-start)
    (ROOT/'docs/validation/v202/evidence_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(count=len(rows),bytes=report['bytes'],localBytes=sum(r['bytes'] for r in local),protected=29,elapsedSeconds=report['elapsedSeconds'])))


if __name__=='__main__':main()
