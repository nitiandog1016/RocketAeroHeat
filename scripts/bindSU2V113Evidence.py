"""绑定v113真实网格/引擎/UI任务与原源版本；大体网格和重启动留本地。"""
from pathlib import Path
import hashlib,json,subprocess,os
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v113_validation';DEST=ROOT/'docs/validation/v113'
FILES=['buildSU2HybridVolume.m','readSU2HybridVolume.m','buildSU2WallProjection.m','writeSU2RocketCase.m',
       'rocketAeroHeatWorkbench.m','su2WorkbenchMeshOptions.m','aeroHeatWorkflowState.m','embeddedSU2Job.m',
       'runEmbeddedSU2Case.m','buildRocketTestBody.m','mapSU2WallResult.m','readSU2SurfaceResult.m',
       'auditSU2ResultConvergence.m','identifySU2NativeCheckpoint.m','createSU2ResearchView.m',
       'createAerothermalGeometryView.m','updateAerothermalGeometryView.m','createCaseConfig.m',
       'estimateSU2CheckpointStorage.m','readSU2TimeoutCheckpoint.m','describeSU2Convergence.m',
       'scripts/installGmshRuntime.py','scripts/meshSU2HybridVolume.py','scripts/runEmbeddedSU2.py',
       'scripts/runMatlabResearch.py','scripts/bindSU2V113Evidence.py','tests/testSU2HybridVolume.m',
       'tests/testSU2HybridWorkbench.m','tests/testSU2HybridRejection.m','tests/auditSU2V113.py','tests/auditSU2V110.py',
       'web/geometry.js','web/geometry.css','web/geometry-su2.html','web/ale-section-core.js','web/su2-section-view.js',
       'docs/validation/v113/PLAN.md']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def save(b,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==b,str(p)
    else:p.write_bytes(b)
    return p.relative_to(ROOT).as_posix()
def copy(p,d):return save(p.read_bytes(),DEST/d)
def epoch(p):return copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))

def main():
    assert not (DEST/'delivery_manifest.json').exists()
    for n in FILES:epoch(ROOT/n)
    versions={sha(p):p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw')}
    previous={sha(p):p for p in (ROOT/'docs/validation/v112/source_versions').glob('*.raw')}
    jobs=[]
    for p in sorted(OUT.glob('*_process.json')):
        j=read(p);assert j['status']=='terminal' and all(j['sourceStableAfter'].values())
        sources={n:h for n,h in j['sourceSHA256'].items() if n in FILES}
        for n,h in sources.items():
            if h not in versions:
                if h in previous:versions[h]=epoch(previous[h])
                else:
                    b=subprocess.check_output(['git','show','HEAD:'+n],cwd=str(ROOT))
                    if hashlib.sha256(b).hexdigest()!=h:b=b.replace(b'\n',b'\r\n')
                    assert hashlib.sha256(b).hexdigest()==h,(n,h)
                    versions[h]=save(b,DEST/'source_versions'/(Path(n).name+'.'+h[:16]+'.raw'))
        for n in j['explicitInputs']:assert sha(ROOT/n)==j['sourceSHA256'][n]
        copy(p,Path('processes')/p.name);log=OUT/p.name.replace('_process.json','_stdout.log');copy(log,Path('processes')/log.name)
        jobs.append(dict(name=p.stem,exitCode=j['exitCode'],sourceStable=True,sources={n:versions[h] for n,h in sources.items()}))
    assert len(jobs)==4 and sum(j['exitCode']==0 for j in jobs)==3
    ui=read(OUT/'ui01/report.json');feature=Path(ui['engineFolder'])
    meshes=list((ROOT/'outputs/hybrid_mesh').glob('*/mesh_process.json'));assert len(meshes)==1
    meshDirs=[OUT/'actual01',OUT/'actual02',meshes[0].parent];meshRecords=[]
    for label,folder in zip(['actual_failed','actual','feature'],meshDirs):
        r=read(folder/'mesh_process.json');assert r['status']=='terminal' and r['inputStable'] and r['sourceStable']
        assert r['scriptSHA256'] in versions
        assert sha(Path(r['input']))==r['inputSHA256']
        for n,h in r['outputSHA256'].items():assert sha(folder/n)==h
        for p in folder.iterdir():
            if p.suffix in ('.json','.log'):copy(p,Path('mesh')/label/p.name)
        copy(Path(r['input']),Path('inputs')/(label+'.mat'))
        meshRecords.append(dict(name=label,exitCode=r['exitCode'],seconds=r['elapsedSeconds'],script=versions[r['scriptSHA256']]))
    engines=[]
    for label,folder in [('actual',OUT/'actual_review02/engine'),('feature',feature)]:
        r=read(folder/'engine_process.json');assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        for n,h in r['inputSHA256'].items():assert sha(Path(n))==h
        for n,h in r['outputSHA256'].items():assert sha(folder/n)==h
        for p in folder.iterdir():
            if p.suffix in ('.json','.log','.cfg') or p.name in ('history.csv','surface_flow.vtk','mapped_research.mat','input_snapshot.mat'):
                copy(p,Path('engine')/label/p.name)
        engines.append(dict(name=label,exitCode=r['exitCode'],seconds=r['elapsedSeconds'],physicsQualified=False))
    audit=read(OUT/'independent_audit.json')
    for n,h in audit['inputSHA256'].items():assert sha(ROOT/n)==h
    for name in ('independent_audit.json','actual_input_source.json','storage_receipt.json'):
        copy(OUT/name,name)
    for name in ('actual_review02/report.json','actual_review02/mapped_result.mat','guard01/checks.json',
                 'ui01/report.json','ui01/feature_geometry.mat','ui01/fields.csv','ui01/fields.csv.json',
                 'ui01/section.csv','ui01/section.csv.json','ui01/capture_done.json','references/source.json'):
        copy(OUT/name,Path('checks')/name)
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(DEST.rglob('*')) if p.is_file()]
    roots=[OUT,feature,meshes[0].parent];paths=[p for d in roots for p in d.rglob('*') if p.is_file()]+[Path(read(meshes[0])['input'])]
    unique={(p.stat().st_dev,p.stat().st_ino):p.stat().st_size for p in paths};assert sum(unique.values())<100*1024**2
    m=dict(status='tested-original-surface-hybrid-mesh-and-main-workbench; two-step-flow-unqualified',
           baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip(),
           matlabJobs=jobs,meshJobs=meshRecords,engineJobs=engines,sourceVersions=versions,evidence=evidence,
           documentationSHA256={n:sha(ROOT/n) for n in ['AGENTS.md','docs/PROJECT_GOALS.md','docs/SU2_HYBRID_ORIGINAL_SURFACE_V113.md']},
           evidenceBytes=sum(e['bytes'] for e in evidence),researchOutputInodeUniqueBytes=sum(unique.values()),
           researchOutputPathBytes=sum(p.stat().st_size for p in paths),runtimeStoredSeparately=True,
           newQualifiedArea=0,physicsQualified=False,largeMeshAndRestartBackedUpOnGitHub=False)
    save((json.dumps(m,ensure_ascii=False,indent=2)+'\n').encode('utf-8'),DEST/'delivery_manifest.json')
    print(json.dumps({k:m[k] for k in ['evidenceBytes','researchOutputInodeUniqueBytes','newQualifiedArea']}))

if __name__=='__main__':main()
