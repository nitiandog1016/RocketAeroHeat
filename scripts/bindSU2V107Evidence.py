"""绑定v107全部终态、失败、准确源码及小型重开输入；不启动求解。"""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v107_validation'
DEST=ROOT/'docs/validation/v107'
FILES=['aeroHeatWorkflowState.m','rocketAeroHeatWorkbench.m','embeddedSU2Job.m','createSU2ResearchView.m',
       'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','writeSU2RocketCase.m','buildSU2RocketVolume.m',
       'buildSU2WallProjection.m','mapSU2WallResult.m','readSU2SurfaceResult.m','buildRocketTestBody.m',
       'createAerothermalGeometryView.m','updateAerothermalGeometryView.m','web/geometry.js','web/geometry.html',
       'web/geometry.css','tests/testEmbeddedSU2Workbench.m','tests/testSU2ResearchReplay.m',
       'tests/testEmbeddedSU2Lifecycle.py','tests/testSU2AsyncBinding.m','tests/testAeroHeatWorkflowState.m',
       'tests/auditSU2V107.py','scripts/bindSU2V107Evidence.py','createCaseConfig.m','calcFreestream.m']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def copy(p,target):
    dst=DEST/target;dst.parent.mkdir(parents=True,exist_ok=True)
    if dst.exists():
        assert sha(dst)==sha(p),str(dst)
    else:
        dst.write_bytes(p.read_bytes())
    return dict(path=dst.relative_to(ROOT).as_posix(),sha256=sha(dst),bytes=dst.stat().st_size)


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    for name in FILES:
        p=ROOT/name;copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))
    versions={sha(p):p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw')}
    jobs=[]
    for file in sorted(OUT.glob('*_process.json')):
        state=json.loads(file.read_text());assert state['status']=='terminal' and all(state['sourceStableAfter'].values())
        chosen={p:h for p,h in state['sourceSHA256'].items() if p in FILES}
        assert all(h in versions for h in chosen.values()),[(p,h) for p,h in chosen.items() if h not in versions]
        name=file.name[:-len('_process.json')]
        jobs.append(dict(name=name,exitCode=state['exitCode'],sourceVersions={p:versions[h] for p,h in chosen.items()},
                         receipt=copy(file,Path('processes')/file.name),
                         stdout=copy(OUT/(name+'_stdout.log'),Path('processes')/(name+'_stdout.log'))))
    assert len(jobs)==4 and sum(j['exitCode']==0 for j in jobs)==3
    ui=json.loads((OUT/'workbench02/report.json').read_text(encoding='utf8'))
    successful=Path(ui['actualEngineFolder']).resolve()
    folders=[successful,Path(ui['cancelFolder']).resolve(),Path(ui['closeFolder']).resolve(),OUT/'binding01',
             OUT/'lifecycle01/invalid_config',OUT/'lifecycle01/time_limit']
    engines=[]
    for i,folder in enumerate(folders):
        state=json.loads((folder/'engine_process.json').read_text(encoding='utf8'))
        assert state['status']=='terminal'
        changed=[p for p,valid in state['inputsStable'].items() if not valid]
        assert (len(changed)==1 and changed[0].endswith('input_snapshot.mat')) if i==3 else not changed
        for file,digest in state['inputSHA256'].items():
            p=Path(file)
            assert sha(p)==digest or digest in versions or str(p) in changed,(file,digest)
        for name,digest in state['outputSHA256'].items():
            assert sha(folder/name)==digest,(folder,name)
        target=Path('engines')/('case%02d'%i)
        for name in ('case.cfg','case_manifest.json','engine_request.json','engine_process.json',
                     'engine_stdout.log','supervisor.log','engine_cancel.json','history.csv'):
            p=folder/name
            if p.is_file():copy(p,target/name)
        engines.append(dict(folder=folder.relative_to(ROOT).as_posix(),exitCode=state['exitCode'],
                            cancelled=state['cancelled'],timedOut=state['timedOut'],intentionalChangedInputs=changed))
    assert len(engines)==6 and sum(e['cancelled'] for e in engines)==3 and sum(e['timedOut'] for e in engines)==1
    for name in ('mapped_research.mat','input_snapshot.mat','engine_process.json','case_manifest.json','case.cfg','history.csv','surface_flow.vtk'):
        copy(successful/name,Path('replay_seed')/name)
    # 同一个网格已由v106交付；避免再保存相同的1.43MB数据。
    shared=ROOT/'docs/validation/v106/inputs/rocket_mesh.su2';assert sha(shared)==sha(successful/'mesh.su2')
    for name in ('workbench02/report.json','workbench02/capture_done.json','replay01/report.json','replay01/capture_done.json',
                 'replay01/research_snapshot.csv','replay01/research_snapshot.csv.json','lifecycle01/report.json','binding01/binding_report.json'):
        copy(OUT/name,Path('checks')/name)
    audit=json.loads((DEST/'independent_audit.json').read_text(encoding='utf8'))
    for name,digest in audit['inputSHA256'].items():assert sha(ROOT/name)==digest
    roots=[OUT]+folders[:3]
    output_bytes=sum(p.stat().st_size for folder in roots for p in folder.rglob('*') if p.is_file())
    assert output_bytes<100*1024*1024
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
              for p in sorted(DEST.rglob('*')) if p.is_file()]
    result=dict(status='tested-main-workbench-and-owned-engine-lifecycle; research-only',matlabJobs=jobs,engineJobs=engines,
                sourceVersions=versions,evidence=evidence,sharedMesh=dict(path=shared.relative_to(ROOT).as_posix(),sha256=sha(shared)),
                outputBytes=output_bytes,evidenceBytes=sum(x['bytes'] for x in evidence),physicsQualified=False,newQualifiedRocketArea=0,
                nativeMainWindowInspected=True,arbitrarySectionForEmbeddedBranchImplemented=False,
                warmStartPhysicsTested=False,originalMeshAndBoundaryImproved=False)
    (DEST/'delivery_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:result[k] for k in ('outputBytes','evidenceBytes','newQualifiedRocketArea')}))


if __name__=='__main__':
    main()
