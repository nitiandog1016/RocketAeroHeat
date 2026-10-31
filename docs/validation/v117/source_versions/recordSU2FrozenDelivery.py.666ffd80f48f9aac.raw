"""封存v117构建、机制诊断、候选和原面显示的小型证据；不上传大体场或编译缓存。"""
from pathlib import Path
import json,hashlib,shutil
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v117_validation';DEST=ROOT/'docs/validation/v117'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def copy(p,q):
    q.parent.mkdir(parents=True,exist_ok=True)
    if q.exists():assert sha(p)==sha(q)
    else:shutil.copyfile(str(p),str(q))
def archive(p):
    h=sha(p);q=DEST/'source_versions'/(p.name+'.'+h[:16]+'.raw');copy(p,q)
    return dict(sha256=h,archived=q.relative_to(ROOT).as_posix())
def main():
    target=DEST/'delivery_manifest.json';assert not target.exists()
    used=['scripts/buildSU2FrozenEnergyRuntime.py','scripts/buildSU2ThermalRuntime.py','scripts/prepareSU2ThermalBuild.py',
          'scripts/prepareSU2FrozenEnergyStudy.py','scripts/su2ThermalRuntime.py','scripts/runEmbeddedSU2.py',
          'scripts/runMatlabResearch.py','scripts/recordSU2FrozenDelivery.py',
          'third_party/su2_thermal/RAHFrozenEnergyRecovery.hpp','third_party/su2_thermal/RAHFrozenEnergyDiagnostics.hpp',
          'third_party/su2_thermal/FrozenEnergyHook.cpp.inc','tests/testSU2FrozenEnergyRecovery.cpp',
          'tests/testSU2FrozenRuntime.py','tests/auditSU2FrozenDirections.py','tests/auditSU2FrozenResult.py',
          'tests/auditSU2V110.py','tests/reviewSU2FrozenEnergyStudy.m','describeSU2Runtime.m','createSU2ResearchView.m',
          'readSU2SurfaceResult.m','mapSU2WallResult.m','auditSU2ResultConvergence.m','identifySU2NativeCheckpoint.m',
          'attachSU2ThermalBoundary.m','resolveSU2ThermalBoundary.m','describeSU2Convergence.m','buildSchaeferFergusonBody.m',
          'createAerothermalGeometryView.m','updateAerothermalGeometryView.m']
    used+=['web/'+n for n in ['geometry-su2.html','geometry.html','geometry.css','geometry.js','ale-section-core.js','su2-section-view.js']]
    sources={rel:archive(ROOT/rel) for rel in used}
    # 编译输入全清单及固定源补丁保留；归档本身有公开固定来源，运行时本地不分发。
    for name in ['extract_compiler','unit_compile','unit_run','configure','compile']:
        rec=read(OUT/(name+'_process.json'));assert rec['status']=='terminal' and rec['exitCode']==0 and not rec['timedOut']
        for suffix in ['_process.json','.log']:copy(OUT/(name+suffix),DEST/'build'/(name+suffix))
    for name in ['source_tree.json','source_patch.diff','build_sources.json','runtime_installed.json','runtime_contract_tests.log']:
        copy(OUT/name,DEST/'build'/name)
    cases={}
    gate=read(OUT/'diagnostic_gate.json')
    for name in (['diagnostic','apply'] if gate['passed'] else ['diagnostic']):
        folder=OUT/name;r=read(folder/'engine_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        for p,h in r['inputSHA256'].items():assert sha(Path(p))==h
        for p,h in r['outputSHA256'].items():assert sha(folder/p)==h
        for file in ['case.cfg','case_manifest.json','engine_request.json','engine_process.json','engine_stdout.log','history.csv',
                     'rah_frozen_counts.csv','rah_frozen_directions.csv','surface_flow.vtk','projection.mat']:
            copy(folder/file,DEST/'cases'/name/file)
        cases[name]=dict(seconds=r['elapsedSeconds'],exitCode=r['exitCode'],mode=r['frozenEnergyMode'],physicsQualified=False)
    job=read(OUT/'review01_process.json');assert job['status']=='terminal' and job['exitCode']==0 and all(job['sourceStableAfter'].values())
    for rel in used:
        if rel in job['sourceSHA256']:assert sha(ROOT/rel)==job['sourceSHA256'][rel]
    for suffix in ['_process.json','.log','_stdout.log']:copy(OUT/('review01'+suffix),DEST/'processes'/('review01'+suffix))
    for file in (OUT/'review01').iterdir():
        if file.is_file() and file.suffix in ['.mat','.csv','.json']:copy(file,DEST/'checks/review01'/file.name)
    for name in ['diagnostic_gate.json','independent_audit.json','launch_gate_rejection.json','preflight_preservation.json']:
        copy(OUT/name,DEST/'checks'/name)
    result=dict(baseCommit='3f208cbfdb54fe92764b0abc8595b96851d4bba7',
        status='tested-bounded-mechanism-study; experimental-not-qualified',newQualifiedArea=0,physicsQualified=False,
        candidateAuthorized=gate['passed'],candidateRun=('apply' in cases),
        sources=sources,engineCases=cases,matlabJob=dict(exitCode=0,sourcesStable=True),
        evidence={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(DEST.rglob('*')) if p.is_file()},
        localOnly=['full meshes and initial/restart conservative fields','controlled executables','fixed build archives','shared historical results'])
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(result['evidence']),sourceFiles=len(sources),engineJobs=len(cases))))
if __name__=='__main__':main()
