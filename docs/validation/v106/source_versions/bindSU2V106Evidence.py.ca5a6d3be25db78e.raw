"""绑定v106全部终态、准确源版本、小型网格与原面映射证据，不启动求解。"""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v106_validation'
DEST=ROOT/'docs/validation/v106'
FILES=['buildSU2FlatPlateCase.m','buildSU2RocketVolume.m','writeSU2RocketCase.m','buildSU2WallProjection.m',
       'mapSU2WallResult.m','runEmbeddedSU2Case.m','readSU2SurfaceResult.m','buildRocketTestBody.m',
       'regionalFlow/axisymmetric/buildRocketAxisymmetricDomain.m','regionalFlow/axisymmetric/buildAxisymmetricFVGeometry.m',
       'tests/auditEmbeddedSU2FlatPlate.py','tests/auditSU2V106.py','tests/studyEmbeddedSU2FlatPlate.m',
       'tests/studySU2RocketVolume.m','tests/startSU2RocketVolume.m','tests/testSU2WallProjection.m',
       'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','scripts/bindSU2V106Evidence.py']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy(source,path):
    target=DEST/path;target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists():
        assert sha(target)==sha(source),str(target)
    else:
        target.write_bytes(source.read_bytes())
    return dict(path=target.relative_to(ROOT).as_posix(),sha256=sha(target),bytes=target.stat().st_size)


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    versions={}
    for path in FILES:
        p=ROOT/path
        copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))
    for folder in [ROOT/'docs/validation/v105/source_versions',DEST/'source_versions']:
        for file in folder.glob('*.raw'):
            versions[sha(file)]=file.relative_to(ROOT).as_posix()
    jobs=[]
    for file in sorted(OUT.glob('*_process.json')):
        state=json.loads(file.read_text());assert state['status']=='terminal' and all(state['sourceStableAfter'].values())
        selected={p:h for p,h in state['sourceSHA256'].items() if p in FILES}
        assert all(h in versions for h in selected.values()),[(p,h) for p,h in selected.items() if h not in versions]
        name=file.name[:-len('_process.json')]
        jobs.append(dict(name=name,exitCode=state['exitCode'],
                         lockedProjectSourceVersions={p:versions[h] for p,h in selected.items()},
                         receipt=copy(file,Path('processes')/file.name),
                         stdout=copy(OUT/(name+'_stdout.log'),Path('processes')/(name+'_stdout.log'))))
    assert len(jobs)==8 and sum(j['exitCode']==0 for j in jobs)==5
    engines=[]
    for file in sorted(OUT.rglob('engine_process.json')):
        state=json.loads(file.read_text());assert state['status']=='terminal' and all(state['inputsStable'].values())
        for name,digest in state['inputSHA256'].items():
            assert sha(Path(name))==digest
        for name,digest in state['outputSHA256'].items():
            assert sha(file.parent/name)==digest
        relative=file.parent.relative_to(OUT)
        evidence=[]
        for name in ['case.cfg','case_manifest.json','engine_request.json','engine_process.json','engine_stdout.log','surface_flow.vtk']:
            source=file.parent/name
            if source.exists():
                evidence.append(copy(source,Path('cases')/relative/name))
        engines.append(dict(name=relative.as_posix(),exitCode=state['exitCode'],timedOut=state['timedOut'],
                            dryRun=state['request'].get('dryRun',False),elapsedSeconds=state['elapsedSeconds'],evidence=evidence))
    assert len(engines)==6 and sum(e['timedOut'] for e in engines)==1 and sum(e['dryRun'] for e in engines)==2
    # 只保存一份完整3D小型输入；同几何不同攻角配置共享此网格。
    copy(OUT/'rocket_startup02/angle02/mesh.su2',Path('inputs')/'rocket_mesh.su2')
    copy(OUT/'plate_outer01/mesh.su2',Path('inputs')/'plate_outer_mesh.su2')
    copy(OUT/'rocket_projection02/report.json',Path('mapping')/'report.json')
    audit=json.loads((DEST/'independent_audit.json').read_text(encoding='utf-8'))
    for path,digest in audit['inputSHA256'].items():
        assert sha(ROOT/path)==digest
    count=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
    assert count<100*1024*1024
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
              for p in sorted(DEST.rglob('*')) if p.is_file()]
    result=dict(status='tested-limited-domain-and-three-dimensional-interface; not-qualified-flow',
                matlabJobs=jobs,engineJobs=engines,evidence=evidence,sourceVersions=versions,
                outputBytes=count,newQualifiedRocketArea=0,physicsQualified=False,
                threeDStartupOnly=True,allOriginalFacesMapped=2448,mainWorkbenchIntegrated=False,
                fineMeshConvergencePassed=False,outerDomainScreenPassed=True,
                areaConventionDifferenceRetained=True)
    (DEST/'delivery_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(matlabJobs=len(jobs),matlabExit0=sum(j['exitCode']==0 for j in jobs),
                         engineJobs=len(engines),evidenceFiles=len(evidence),outputBytes=count,newQualifiedRocketArea=0)))


if __name__=='__main__':
    main()
