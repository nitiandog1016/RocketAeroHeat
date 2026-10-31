"""绑定经典热边界、实际任务、最终UI和面积勘误，保留失败及准确源版本。"""
from pathlib import Path
import hashlib,json,subprocess

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v114_validation';DEST=ROOT/'docs/validation/v114'
FILES=['buildSchaeferFergusonBody.m','schaeferFergusonCondition.m',
       'resolveSU2ThermalBoundary.m','attachSU2ThermalBoundary.m','writeSU2RocketCase.m',
       'rocketAeroHeatWorkbench.m','su2WorkbenchMeshOptions.m','aeroHeatWorkflowState.m',
       'buildSU2HybridVolume.m','readSU2HybridVolume.m','buildSU2WallProjection.m',
       'embeddedSU2Job.m','mapSU2WallResult.m','readSU2SurfaceResult.m',
       'auditSU2ResultConvergence.m','identifySU2NativeCheckpoint.m','createSU2ResearchView.m',
       'createAerothermalGeometryView.m','updateAerothermalGeometryView.m','createCaseConfig.m',
       'estimateSU2CheckpointStorage.m','readSU2TimeoutCheckpoint.m','describeSU2Convergence.m',
       'scripts/meshSU2HybridVolume.py','scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py',
       'scripts/bindSU2V114Evidence.py','scripts/mergeSU2V114Duplicates.ps1',
       'tests/testSchaeferThermalContract.m','tests/testSchaeferAdiabaticWorkbench.m',
       'tests/reviewSU2ThermalView.m','tests/testSU2HybridAreaMetadata.m',
       'tests/auditSU2V114.py','tests/auditSU2V110.py',
       'web/geometry.js','web/geometry.css','web/geometry-su2.html',
       'web/ale-section-core.js','web/su2-section-view.js']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def save(data,p):
    p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists():assert p.read_bytes()==data,str(p)
    else:p.write_bytes(data)
    return p.relative_to(ROOT).as_posix()
def copy(p,d):return save(p.read_bytes(),DEST/d)
def epoch(p):return copy(p,Path('source_versions')/(p.name+'.'+sha(p)[:16]+'.raw'))

def main():
    assert not (DEST/'delivery_manifest.json').exists()
    for n in FILES:epoch(ROOT/n)
    versions={sha(p):p.relative_to(ROOT).as_posix() for p in (DEST/'source_versions').glob('*.raw')}
    previous={sha(p):p for p in (ROOT/'docs/validation/v113/source_versions').glob('*.raw')}
    jobs=[]
    for p in sorted(OUT.glob('*_process.json')):
        j=read(p);assert j['status']=='terminal' and j['exitCode']==0
        assert all(j['sourceStableAfter'].values())
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
        copy(p,Path('processes')/p.name)
        log=OUT/p.name.replace('_process.json','_stdout.log');copy(log,Path('processes')/log.name)
        jobs.append(dict(name=p.stem,exitCode=0,sourceStable=True,
                         sources={n:versions[h] for n,h in sources.items()}))
    assert len(jobs)==4
    engine=Path(read(OUT/'actual01/report.json')['engineFolder'])
    mesh=ROOT/'outputs/hybrid_mesh/a75f5d87-b861-4478-b972-3095018283e6'
    mr=read(mesh/'mesh_process.json')
    assert mr['status']=='terminal' and mr['exitCode']==0 and mr['sourceStable'] and mr['inputStable']
    assert sha(Path(mr['input']))==mr['inputSHA256'] and mr['scriptSHA256'] in versions
    for n,h in mr['outputSHA256'].items():assert sha(mesh/n)==h
    for p in mesh.iterdir():
        if p.suffix in ('.json','.log'):copy(p,Path('mesh')/p.name)
    copy(Path(mr['input']),Path('inputs/classic_surface.mat'))
    er=read(engine/'engine_process.json')
    assert er['status']=='terminal' and er['exitCode']==0 and all(er['inputsStable'].values())
    for n,h in er['inputSHA256'].items():assert sha(Path(n))==h
    for n,h in er['outputSHA256'].items():assert sha(engine/n)==h
    for p in engine.iterdir():
        if p.suffix in ('.json','.cfg','.log') or p.name in (
                'history.csv','surface_flow.vtk','mapped_research.mat','input_snapshot.mat'):
            copy(p,Path('engine')/p.name)
    audit=read(OUT/'independent_audit.json')
    for n,h in audit['inputSHA256'].items():assert sha(ROOT/n)==h
    assert audit['unmatchedStationCount']==1 and not audit['experimentalErrorCalculated']
    for n in ['independent_audit.json','storage_receipt.json',
              'unit01/geometry.mat','unit01/report.json',
              'actual01/classic_result.mat','actual01/report.json','actual01/ui_checks.json',
              'actual01/fields.csv','actual01/fields.csv.json',
              'actual01/section.csv','actual01/section.csv.json','actual01/capture_done.json',
              'review01/report.json','review01/checks.json','review01/fields.csv',
              'review01/fields.csv.json','review01/section.csv','review01/section.csv.json',
              'review01/legacy_isothermal.csv','review01/legacy_isothermal.csv.json',
              'area01/report.json','area01/corrected_input/case_manifest.json']:
        copy(OUT/n,Path('checks')/n)
    storage=read(OUT/'storage_receipt.json')
    for pair in storage['pairs']:
        a=Path(pair['target']);b=Path(pair['retained'])
        assert sha(a)==sha(b)==pair['sha256']
        assert (a.stat().st_dev,a.stat().st_ino)==(b.stat().st_dev,b.stat().st_ino)
    paths=[p for d in [OUT,engine,mesh] for p in d.rglob('*') if p.is_file()]+[Path(mr['input'])]
    unique={(p.stat().st_dev,p.stat().st_ino):p.stat().st_size for p in paths}
    assert sum(unique.values())<100*1024**2
    evidence=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
              for p in sorted(DEST.rglob('*')) if p.is_file()]
    docs=['AGENTS.md','docs/PROJECT_GOALS.md','docs/SU2_CLASSIC_THERMAL_BOUNDARY_V114.md',
          'tests/data/schaefer_ferguson_1961/pressure_fig8.json']
    m=dict(status='tested-classic-thermal-workflow-and-area-metadata; flow-unconverged',
           baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip(),
           matlabJobs=jobs,meshSeconds=mr['elapsedSeconds'],engineSeconds=er['elapsedSeconds'],
           sourceVersions=versions,evidence=evidence,documentationSHA256={n:sha(ROOT/n) for n in docs},
           evidenceBytes=sum(x['bytes'] for x in evidence),researchOutputInodeUniqueBytes=sum(unique.values()),
           researchOutputPathBytes=sum(p.stat().st_size for p in paths),runtimeStoredSeparately=True,
           spatialTemperatureExecuted=False,adiabaticHeatFluxIsPrescribed=True,
           unmatchedPressureStationCount=1,newQualifiedArea=0,physicsQualified=False,
           largeMeshAndRestartBackedUpOnGitHub=False)
    save((json.dumps(m,ensure_ascii=False,indent=2)+'\n').encode('utf-8'),DEST/'delivery_manifest.json')
    print(json.dumps({k:m[k] for k in ['evidenceBytes','researchOutputInodeUniqueBytes','newQualifiedArea']}))

if __name__=='__main__':main()
