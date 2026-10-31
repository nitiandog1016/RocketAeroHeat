"""绑定v109三维初始化研究和主工作台截面证据；不启动或延长数值运行。"""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v109_validation';DEST=ROOT/'docs/validation/v109'
FILES=['scripts/prepareSU2RocketInitialState.py','tests/auditSU2RocketInitialState.py','tests/auditSU2V109.py',
       'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','scripts/bindSU2V109Evidence.py',
       'createSU2ResearchView.m','createAerothermalGeometryView.m','updateAerothermalGeometryView.m',
       'rocketAeroHeatWorkbench.m','aeroHeatWorkflowState.m','su2WorkbenchMeshOptions.m',
       'buildSU2RocketVolume.m','buildSU2WallProjection.m','mapSU2WallResult.m','readSU2SurfaceResult.m',
       'web/geometry.js','web/geometry.html','web/geometry-su2.html','web/geometry.css','web/ale-section-core.js','web/su2-section-view.js',
       'tests/testSU2Sections.m','tests/testSU2SectionView.js','tests/replaySU2SectionUI.m','tests/reviewSU2RocketCheckpoint.m']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def save(data,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==data,str(p)
    else:p.write_bytes(data)
    return p.relative_to(ROOT).as_posix()


def copy(p,dst):return save(p.read_bytes(),DEST/dst)


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    for name in FILES:
        p=ROOT/name;copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))
    versions={sha(p):p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw')}
    jobs=[]
    for p in sorted(OUT.glob('*_process.json')):
        s=json.loads(p.read_text());assert s['status']=='terminal' and all(s['sourceStableAfter'].values())
        selected={name:h for name,h in s['sourceSHA256'].items() if name in FILES}
        for name,h in selected.items():
            if h not in versions:
                raw=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(ROOT));assert hashlib.sha256(raw).hexdigest()==h,(name,h)
                versions[h]=save(raw,DEST/'source_versions'/(Path(name).name+'.'+h[:16]+'.raw'))
        name=p.name[:-len('_process.json')];copy(p,Path('processes')/p.name);copy(OUT/(name+'_stdout.log'),Path('processes')/(name+'_stdout.log'))
        jobs.append(dict(name=name,exitCode=s['exitCode'],sources={p:versions[h] for p,h in selected.items()}))
    assert len(jobs)==3 and sum(j['exitCode']==0 for j in jobs)==2
    folder=OUT/'rocket_seed01';engine=json.loads((folder/'engine_process.json').read_text())
    assert engine['status']=='terminal' and engine['timedOut'] and not engine['cancelled'] and all(engine['inputsStable'].values())
    for p,h in engine['inputSHA256'].items():assert sha(Path(p))==h
    for name,h in engine['outputSHA256'].items():assert sha(folder/name)==h
    for name in ('case.cfg','case_manifest.json','source_case_manifest.json','engine_process.json','engine_request.json',
                 'engine_stdout.log','initialization_manifest.json','initial_audit.json','surface_flow.vtk'):
        copy(folder/name,Path('engine')/name)
    # 保存稀疏历史原行用于检查点/末段审核；不把它冒充完整3922行历史。
    lines=(folder/'history.csv').read_bytes().splitlines(keepends=True)
    wanted={0,1,2,3750,3921}|set(range(0,3922,250))|set(range(3671,3922))
    selected=[line for line in lines[1:] if int(line.split(b',')[2]) in wanted]
    save(lines[0]+b''.join(selected),DEST/'engine'/'history_selected_rows.csv')
    for sub in ('checkpoint01','sections02','final_ui01'):
        for p in (OUT/sub).glob('*'):
            if p.is_file():copy(p,Path('checks')/sub/p.name)
    audit=json.loads((DEST/'independent_audit.json').read_text())
    for p,h in audit['inputSHA256'].items():assert sha(ROOT/p)==h
    seed=ROOT/'docs/validation/v107/replay_seed';shared=[]
    for name in ('mapped_research.mat','input_snapshot.mat','engine_process.json'):
        p=seed/name;shared.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)))
    size=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file());assert size<100*1024*1024
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(DEST.rglob('*')) if p.is_file()]
    report=dict(status='tested-main-section-UI-and-improved-startup; three-dimensional-run-timed-out-unconverged',
                matlabJobs=jobs,engine=dict(exitCode=engine['exitCode'],timedOut=True,elapsedSeconds=engine['elapsedSeconds'],
                    completeHistorySHA256=sha(folder/'history.csv'),completeHistoryRows=len(lines)-1,historyPublication='selected-exact-rows-only'),
                sourceVersions=versions,evidence=evidence,sharedUIReplay=shared,outputBytes=size,evidenceBytes=sum(x['bytes'] for x in evidence),
                currentRunShownInFinalUIScreenshot=False,UIScreenshotSource='actual v107 replay; M4 alpha20 beta5',
                checkpointIteration=3750,terminalIteration=3921,residualConverged=False,newQualifiedRocketArea=0,
                largeInitialStateAndRestartBackedUpOnGitHub=False)
    (DEST/'delivery_manifest.json').write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    print(json.dumps({k:report[k] for k in ('outputBytes','evidenceBytes','newQualifiedRocketArea')}))


if __name__=='__main__':main()
