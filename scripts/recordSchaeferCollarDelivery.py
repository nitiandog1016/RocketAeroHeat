"""封存v116实际进程、输入及准确源码；失败和未合格状态不改写为成功。"""
from pathlib import Path
import hashlib,json,shutil
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v116_validation';DEST=ROOT/'docs/validation/v116'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def copy(p,q):
    q.parent.mkdir(parents=True,exist_ok=True)
    if q.exists():assert sha(p)==sha(q),str(q)
    else:shutil.copyfile(str(p),str(q))
def archive(p):
    h=sha(p);q=DEST/'source_versions'/(p.name+'.'+h[:16]+'.raw');copy(p,q)
    return dict(sha256=h,archived=q.relative_to(ROOT).as_posix())
def main():
    target=DEST/'delivery_manifest.json';assert not target.exists()
    used=['readSchaeferHeatTransferData.m','tests/prepareSchaeferCollarStudy.m','tests/reviewSchaeferCollarStudy.m',
          'scripts/transferSU2HybridState.py','scripts/runEmbeddedSU2.py','scripts/su2ThermalRuntime.py',
          'scripts/runMatlabResearch.py','scripts/meshSU2HybridVolume.py','buildSU2HybridVolume.m',
          'readSU2HybridVolume.m','buildSchaeferFergusonBody.m','writeSU2RocketCase.m',
          'buildSU2WallProjection.m','readSU2SurfaceResult.m','mapSU2WallResult.m',
          'attachSU2ThermalBoundary.m','resolveSU2ThermalBoundary.m','auditSU2ResultConvergence.m',
          'identifySU2NativeCheckpoint.m','createSU2ResearchView.m','describeSU2Convergence.m',
          'createAerothermalGeometryView.m','updateAerothermalGeometryView.m']
    used += ['web/'+n for n in ['geometry-su2.html','geometry.html','geometry.css','geometry.js','ale-section-core.js','su2-section-view.js']]
    used += ['tests/auditSchaeferCollarPreparation.py','tests/auditSU2V116.py','tests/auditSU2V110.py',
             'tests/data/schaefer_ferguson_1961/heat_fig11.json','scripts/recordSchaeferCollarDelivery.py']
    sources={rel:archive(ROOT/rel) for rel in used};jobs=[]
    for name in ['prepare01','review01','review02']:
        receipt=read(OUT/(name+'_process.json'))
        assert receipt['status']=='terminal' and all(receipt['sourceStableAfter'].values())
        bound={}
        for rel in used:
            if rel not in receipt['sourceSHA256']:continue
            h=receipt['sourceSHA256'][rel];current=ROOT/rel
            if sha(current)==h:bound[rel]=archive(current)
            else:
                versions=[p for p in DEST.glob('*.raw') if sha(p)==h]
                assert len(versions)==1,(name,rel,h)
                bound[rel]=dict(sha256=h,archived=versions[0].relative_to(ROOT).as_posix())
        jobs.append(dict(name=name,exitCode=receipt['exitCode'],sourceStable=True,sources=bound))
        for suffix in ['_process.json','.log','_stdout.log']:
            copy(OUT/(name+suffix),DEST/'processes'/(name+suffix))
    cases={}
    for name in ['control','candidate']:
        folder=OUT/'study01'/name;r=read(folder/'engine_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        for p,h in r['inputSHA256'].items():assert sha(Path(p))==h
        for p,h in r['outputSHA256'].items():assert sha(folder/p)==h
        for file in ['case.cfg','case_manifest.json','engine_request.json','engine_process.json','engine_stdout.log','history.csv','surface_flow.vtk','projection.mat']:
            copy(folder/file,DEST/'cases'/name/file)
        if (folder/'initial.transfer.json').exists():copy(folder/'initial.transfer.json',DEST/'cases'/name/'initial.transfer.json')
        cases[name]=dict(seconds=r['elapsedSeconds'],exitCode=0,physicalSolutionQualified=False)
    for file in ['gmsh.log','mesh_process.json','mesh_report.json']:
        copy(OUT/'study01/mesh32'/file,DEST/'mesh32'/file)
    mesher=read(OUT/'study01/mesh32/mesh_process.json');copy(Path(mesher['input']),DEST/'mesh32/input.mat')
    for file in ['preparation.json','preflight_audit.json']:
        copy(OUT/'study01'/file,DEST/'checks'/file)
    for file in ['initial_comparison.json','initial_wall_profiles.json','independent_audit.json','duplicate_review.json']:
        copy(OUT/file,DEST/'checks'/file)
    for part in ['review01','review02']:
        for p in (OUT/part).iterdir():
            if p.is_file() and p.suffix in ['.mat','.csv','.json']:copy(p,DEST/'checks'/part/p.name)
    manifest=dict(baseCommit='cb36c83e96a8dbc07b08d51a94132a0aac5fa8b0',
        status='tested-interface-and-paired-screening; candidate-not-adopted; thermal-data-draft',
        matlabJobs=jobs,engineJobs=cases,mesherTerminalExitCode=mesher['exitCode'],sources=sources,
        newQualifiedArea=0,physicalSolutionQualified=False,nominalStationMatchedPairs=0,
        evidence={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(DEST.rglob('*')) if p.is_file()},
        localOnly=['full old and new hybrid meshes','full restart and initial conservative fields','runtimes and historical object store'])
    target.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(manifest['evidence']),sourceFiles=len(sources),matlabJobs=len(jobs),engineJobs=len(cases))))
if __name__=='__main__':main()
