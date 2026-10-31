"""绑定v120固定外域近壁对照、物性能力探测与原面UI，不复制完整体场。"""
from pathlib import Path
import json,hashlib,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v120_validation';DEST=ROOT/'docs/validation/v120'
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
    h=sha(p);q=DEST/'source_versions'/(p.name+'.'+h[:16]+'.raw');copy(p,q)
    return dict(sha256=h,archived=q.relative_to(ROOT).as_posix())
def main():
    target=DEST/'delivery_manifest.json';assert not target.exists();sources={};jobs={}
    phases={'preparation01':['tests/prepareActualRocketNearWallStudy.m','readSU2HybridVolume.m','writeSU2RocketCase.m','resolveSU2ThermalBoundary.m','buildSU2WallProjection.m'],
      'review01':['tests/reviewActualRocketNearWallStudy.m','createSU2ResearchView.m','auditSU2ResultConvergence.m','identifySU2NativeCheckpoint.m',
                  'describeSU2Convergence.m','describeSU2Runtime.m','readSU2SurfaceResult.m','mapSU2WallResult.m','attachSU2ThermalBoundary.m',
                  'createAerothermalGeometryView.m','updateAerothermalGeometryView.m','web/geometry-su2.html','web/geometry.html',
                  'web/geometry.css','web/geometry.js','web/ale-section-core.js','web/su2-section-view.js']}
    for name,used in phases.items():
        r=read(OUT/(name+'_process.json'));assert r['status']=='terminal' and r['exitCode']==0 and all(r['sourceStableAfter'].values())
        for rel in used:assert sha(ROOT/rel)==r['sourceSHA256'][rel];sources[rel]=archive(ROOT/rel)
        jobs[name]=dict(exitCode=0,sourceStable=True,usedSourceBindings={rel:sources[rel]['sha256'] for rel in used})
        for suffix in ['_process.json','.log','_stdout.log']:copy(OUT/(name+suffix),DEST/'processes'/(name+suffix))
    extra=['scripts/runMatlabResearch.py','scripts/runEmbeddedSU2.py','scripts/su2ThermalRuntime.py','scripts/meshSU2HybridVolume.py',
           'scripts/redistributeSU2NormalCollar.py','scripts/transferSU2NormalColumns.py','scripts/prepareSU2CheckpointContinuation.py',
           'scripts/recordActualNearWallDelivery.py','scripts/mergeActualNearWallDuplicates.py',
           'tests/prepareActualRocketNearWallRuns.py','tests/auditActualRocketNearWallStudy.py','tests/prepareSU2CoolPropProbe.py','tests/auditSU2V110.py',
           'tests/locateNearWallHeatDifference.py',
           'tools/runtime/su2-thermal-build/source/Common/src/CConfig.cpp',
           'tools/runtime/su2-thermal-build/source/SU2_CFD/src/fluid/CCoolProp.cpp',
           'tools/runtime/su2-thermal-build/source/SU2_CFD/include/fluid/CCoolProp.hpp',
           'tools/runtime/su2-thermal-build/source/SU2_CFD/src/solvers/CEulerSolver.cpp',
           'tools/runtime/su2-thermal-build/source/config_template.cfg','tools/runtime/su2-thermal-build/source/meson_options.txt',
           'tools/runtime/su2-thermal-build/source/meson.build']
    for rel in extra:sources[rel]=archive(ROOT/rel)
    mesh=OUT/'mesh32';mr=read(mesh/'mesh_process.json');assert mr['status']=='terminal' and mr['exitCode']==0 and mr['sourceStable'] and not mr['gmshRun']
    for p,h in mr['sourceHashes'].items():assert sha(Path(p))==h
    for p,h in mr['outputSHA256'].items():assert sha(mesh/p)==h
    for p in mesh.iterdir():
        if p.suffix in ['.json','.log'] or p.name=='mesh_input.mat':copy(p,DEST/'mesh'/p.name)
    engines={}
    for name in ['control','candidate','coolprop_probe']:
        folder=OUT/name if name=='coolprop_probe' else OUT/'study01'/name;r=read(folder/'engine_process.json')
        assert r['status']=='terminal' and not r['timedOut'] and not r['cancelled'] and all(r['inputsStable'].values())
        if name!='coolprop_probe':assert r['exitCode']==0
        for p,h in r['inputSHA256'].items():assert sha(Path(p))==h
        for p,h in r['outputSHA256'].items():assert sha(folder/p)==h
        for p in folder.iterdir():
            if p.suffix in ['.cfg','.json','.log'] or p.name in ['history.csv','surface_flow_000100.vtk','surface_flow_000199.vtk']:
                copy(p,DEST/'cases'/name/p.name)
        engines[name]=dict(exitCode=r['exitCode'],elapsedSeconds=r['elapsedSeconds'],dryRun=r['request']['dryRun'],sourceStable=True)
    for p in (OUT/'review01').iterdir():
        if p.is_file() and p.suffix in ['.mat','.json','.csv']:copy(p,DEST/'checks/review01'/p.name)
    for name in ['initial_audit.json','initial_flow.transfer.json','independent_audit.json','heat_difference_location.json','capability_probe.json','storage_receipt.json']:
        copy(OUT/name,DEST/'checks'/name)
    copy(OUT/'study01/preparation.json',DEST/'checks/preparation.json')
    preserved=read(ROOT/'outputs/v117_validation/preflight_preservation.json')
    for p,h in preserved.items():assert sha(ROOT/p)==h
    (DEST/'checks/preservation.json').write_text(json.dumps(dict(checked=len(preserved),allUnchanged=True),indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file()];objects={};inventory=[]
    old={(p.stat().st_dev,p.stat().st_ino) for d in ['v118_validation','v119_validation'] for p in (ROOT/'outputs'/d).rglob('*') if p.is_file()}
    for p in files:
        st=p.stat();assert st.st_ino;key=(st.st_dev,st.st_ino);objects.setdefault(key,st.st_size)
        inventory.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=st.st_size,fileID=str(st.st_ino),sha256=sha(p)))
    storage=dict(pathBytes=sum(p.stat().st_size for p in files),uniqueFileBytes=sum(objects.values()),
                 newComparedWithV118V119Bytes=sum(v for k,v in objects.items() if k not in old),actualNTFSAllocationMeasured=False,
                 explicitCheckpointBudgetMiBPerCase=512,uniqueResultsDeleted=False,v117BlockedCleanupRetried=False,inventory=inventory)
    (DEST/'checks/local_storage.json').write_text(json.dumps(storage,indent=2)+'\n')
    audit=read(OUT/'independent_audit.json')
    result=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip(),
       status='tested-nearwall-pair-and-result-component; research-screen-only',physicsQualified=False,newQualifiedArea=0,
       screeningPassed=audit['allScreeningGatesPassed'],matlabJobs=jobs,engines=engines,gmshRun=False,sources=sources,
       evidence={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(DEST.rglob('*')) if p.is_file()},
       localOnly=['full meshes and conservative restart fields','official/controlled runtime and caches','historical baseline results'])
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(dict(evidenceFiles=len(result['evidence']),sourceVersions=len(sources),uniqueFileBytes=storage['uniqueFileBytes'])))
if __name__=='__main__':main()
