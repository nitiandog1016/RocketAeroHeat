"""绑定v119实际续算、拒绝、物性/近壁及公开参考证据，不复制大型体场。"""
from pathlib import Path
import json,hashlib,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v119_validation';DEST=ROOT/'docs/validation/v119'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()

def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def copy(p,q):
    q.parent.mkdir(parents=True,exist_ok=True)
    if q.exists():assert sha(p)==sha(q)
    else:shutil.copyfile(str(p),str(q))

def archive(p):
    digest=sha(p);q=DEST/'source_versions'/(p.name+'.'+digest[:16]+'.raw');copy(p,q)
    return dict(sha256=digest,archived=q.relative_to(ROOT).as_posix())

def main():
    target=DEST/'delivery_manifest.json';assert not target.exists();sources={};jobs={}
    phases={
      'controls01':['tests/testSU2ContinuationWorkbench.m','rocketAeroHeatWorkbench.m','writeSU2CheckpointInitialState.m',
                    'readSU2HybridVolume.m','aeroHeatWorkflowState.m','estimateSU2CheckpointStorage.m','scripts/prepareSU2CheckpointContinuation.py'],
      'review01':['tests/reviewActualRocketContinuation.m','tests/auditSU2PropertyAssumption.m','airThermalProperties.m',
                  'createSU2ResearchView.m','auditSU2ResultConvergence.m','identifySU2NativeCheckpoint.m','describeSU2Convergence.m',
                  'describeSU2Runtime.m','readSU2SurfaceResult.m','mapSU2WallResult.m','attachSU2ThermalBoundary.m',
                  'createAerothermalGeometryView.m','updateAerothermalGeometryView.m',
                  'web/geometry-su2.html','web/geometry.html','web/geometry.css','web/geometry.js','web/ale-section-core.js','web/su2-section-view.js']}
    for name,used in phases.items():
        r=read(OUT/(name+'_process.json'));assert r['status']=='terminal' and r['exitCode']==0 and all(r['sourceStableAfter'].values())
        for rel in used:
            assert sha(ROOT/rel)==r['sourceSHA256'][rel],(name,rel)
            sources[rel]=archive(ROOT/rel)
        jobs[name]=dict(exitCode=0,sourceStable=True,usedSourceBindings={rel:sources[rel]['sha256'] for rel in used})
        for suffix in ['_process.json','.log','_stdout.log']:copy(OUT/(name+suffix),DEST/'processes'/(name+suffix))
        for file in (OUT/name).iterdir():
            if file.is_file() and file.suffix in ['.json','.csv','.mat']:copy(file,DEST/'checks'/name/file.name)
    for rel in ['scripts/runMatlabResearch.py','scripts/runEmbeddedSU2.py','scripts/su2ThermalRuntime.py',
                'scripts/recordActualContinuationDelivery.py','embeddedSU2Job.m','writeSU2RocketCase.m',
                'tests/testSU2ContinuationPreparation.py','tests/auditActualRocketContinuation.py',
                'tests/auditActualRocketWallSpacing.py','tests/checkPublicDNSReference.py','tests/auditSU2V110.py',
                'tools/runtime/su2-thermal-build/source/Common/src/geometry/CGeometry.cpp',
                'tools/runtime/su2-thermal-build/source/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl',
                'tools/runtime/su2-thermal-build/source/config_template.cfg',
                'tools/runtime/su2-thermal-build/source/Common/src/CConfig.cpp']:
        sources[rel]=archive(ROOT/rel)
    folder=OUT/'continuation01';r=read(folder/'engine_process.json')
    assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and not r['cancelled'] and all(r['inputsStable'].values())
    for p,h in r['inputSHA256'].items():assert sha(Path(p))==h,p
    for p,h in r['outputSHA256'].items():assert sha(folder/p)==h,p
    for name in ['case.cfg','case_manifest.json','continuation.json','engine_request.json','engine_process.json','engine_stdout.log',
                 'history.csv','surface_flow_000150.vtk','surface_flow_000299.vtk']:
        copy(folder/name,DEST/'cases/continuation01'/name)
    for name in ['preparation_checks.json','independent_audit.json','property_comparison.json','wall_spacing.json']:
        copy(OUT/name,DEST/'checks'/name)
    for p in [OUT/'matlab_adapter01/continuation.json',OUT/'matlab_adapter01/case.cfg']:
        copy(p,DEST/'checks/matlab_adapter'/p.name)
    preserved=read(ROOT/'outputs/v117_validation/preflight_preservation.json')
    for p,h in preserved.items():assert sha(ROOT/p)==h,p
    (DEST/'checks/preservation.json').write_text(json.dumps(dict(checked=len(preserved),allUnchanged=True,
        reference='outputs/v117_validation/preflight_preservation.json'),indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file()]
    old={(p.stat().st_dev,p.stat().st_ino) for p in (ROOT/'outputs/v118_validation').rglob('*') if p.is_file()}
    objects={};inventory=[]
    for p in files:
        st=p.stat();assert st.st_ino;key=(st.st_dev,st.st_ino);objects.setdefault(key,st.st_size)
        inventory.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=st.st_size,fileID=str(st.st_ino),sha256=sha(p)))
    storage=dict(paths=len(files),uniqueObjects=len(objects),pathBytes=sum(p.stat().st_size for p in files),
                 uniqueFileBytes=sum(objects.values()),newComparedWithV118Bytes=sum(v for k,v in objects.items() if k not in old),
                 actualNTFSAllocationMeasured=False,storageException='full original-domain native checkpoint research; explicit256MiB checkpoint budget',
                 uniqueResultsDeleted=False,v117BlockedCleanupRetried=False,inventory=inventory)
    (DEST/'checks/local_storage.json').write_text(json.dumps(storage,indent=2)+'\n')
    result=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip(),
                status='tested-native-continuation-and-UI; one-predeclared-gate-failed; not-physically-qualified',
                physicsQualified=False,newQualifiedArea=0,matlabJobs=jobs,engineJobs=1,engineSeconds=r['elapsedSeconds'],sources=sources,
                evidence={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(DEST.rglob('*')) if p.is_file()},
                localOnly=['full meshes and conservative restart fields','official/controlled runtime and caches','historical baseline results'])
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(evidenceFiles=len(result['evidence']),sourceVersions=len(sources),matlabJobs=len(jobs),engineJobs=1,storage=storage['uniqueFileBytes'])))

if __name__=='__main__':main()
