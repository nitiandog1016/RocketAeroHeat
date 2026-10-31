"""绑定v110空间诊断、审核及主UI实际终态；不执行求解或覆盖旧证据。"""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v110_validation';DEST=ROOT/'docs/validation/v110'
FILES=['auditSU2ResultConvergence.m','describeSU2Convergence.m','createSU2ResearchView.m','rocketAeroHeatWorkbench.m',
       'readSU2SurfaceResult.m','mapSU2WallResult.m','buildSU2WallProjection.m','aeroHeatWorkflowState.m',
       'createAerothermalGeometryView.m','updateAerothermalGeometryView.m','embeddedSU2Job.m',
       'web/geometry.js','web/geometry.css','web/geometry-su2.html','web/ale-section-core.js','web/su2-section-view.js',
       'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','scripts/prepareSU2CheckpointProbe.py',
       'scripts/diagnoseSU2Checkpoint.py','scripts/plotSU2CheckpointDiagnostics.py','scripts/bindSU2V110Evidence.py',
       'tests/testSU2ConvergenceAudit.m','tests/auditSU2V110.py']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def save(data,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():assert path.read_bytes()==data,str(path)
    else:path.write_bytes(data)
    return path.relative_to(ROOT).as_posix()


def copy(p,relative):return save(p.read_bytes(),DEST/relative)


def main():
    assert not (DEST/'delivery_manifest.json').exists()
    for name in FILES:
        p=ROOT/name;copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))
    # 初版空间报告仅少三个解释/位置字段；重构原源码后必须命中原报告的完整SHA。
    old=json.loads((OUT/'spatial_report.json').read_text())
    p=ROOT/'scripts/diagnoseSU2Checkpoint.py';data=p.read_bytes()
    data=b''.join(line for line in data.splitlines(keepends=True) if not any(k in line for k in
                 (b'residualScaleCaveat=',b'downstreamOfBody=',b'upstreamOfBody=')))
    h=hashlib.sha256(data).hexdigest();assert h==old['inputSHA256']['scripts/diagnoseSU2Checkpoint.py']
    save(data,DEST/'source_versions'/('diagnoseSU2Checkpoint.py.'+h[:16]+'.raw'))
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
    assert len(jobs)==5 and sum(j['exitCode']==0 for j in jobs)==1
    engines=[]
    for name in ('probe01','probe02'):
        folder=OUT/name;s=json.loads((folder/'engine_process.json').read_text())
        assert s['status']=='terminal' and all(s['inputsStable'].values())
        for path,h in s['inputSHA256'].items():
            if Path(path)==ROOT/'scripts/prepareSU2CheckpointProbe.py':assert h in versions
            else:assert sha(Path(path))==h,path
        for path,h in s['outputSHA256'].items():assert sha(folder/path)==h
        for file in folder.iterdir():
            if file.is_file() and file.suffix in ('.json','.cfg','.log') or file.name=='history.csv':
                copy(file,Path('engine')/name/file.name)
        engines.append(dict(name=name,exitCode=s['exitCode'],elapsedSeconds=s['elapsedSeconds'],timedOut=s['timedOut']))
    assert [x['exitCode'] for x in engines]==[1,0]
    for name in ('spatial_report.json','spatial_report_final.json'):copy(OUT/name,name)
    for folder in sorted(OUT.glob('audit0*')):
        for p in folder.iterdir():
            if p.is_file():copy(p,Path('checks')/folder.name/p.name)
    audit=json.loads((DEST/'independent_audit.json').read_text())
    for p,h in audit['inputSHA256'].items():assert sha(ROOT/p)==h
    output_size=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file());assert output_size<100*1024*1024
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(DEST.rglob('*')) if p.is_file()]
    report=dict(status='tested-spatial-diagnostic-and-main-convergence-identity; three-dimensional-physics-incomplete',
                matlabJobs=jobs,engineJobs=engines,sourceVersions=versions,evidence=evidence,
                outputBytes=output_size,evidenceBytes=sum(x['bytes'] for x in evidence),
                currentDiagnosticShownInUIScreenshot=False,UIScreenshotSource='v107 M4 alpha20 beta5 iteration1',
                checkpointSourceIteration=3750,newQualifiedRocketArea=0,physicsQualified=False,
                largeVolumeAndRestartBackedUpOnGitHub=False,noLongSolverRunThisVersion=True)
    (DEST/'delivery_manifest.json').write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
    print(json.dumps({k:report[k] for k in ('outputBytes','evidenceBytes','newQualifiedRocketArea')}))


if __name__=='__main__':main()
