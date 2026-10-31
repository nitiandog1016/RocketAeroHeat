"""收口v121源码/许可/真实终态与小型复现输入；不复制运行时或完整历史大场。"""
from pathlib import Path
import hashlib,json,shutil,subprocess,ast
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v121_validation';DEST=ROOT/'docs/validation/v121'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def copy(p,q):
    q.parent.mkdir(parents=True,exist_ok=True)
    if q.exists():assert sha(p)==sha(q),str(q)
    else:shutil.copyfile(str(p),str(q))
def archive(p):
    h=sha(p);q=DEST/'source_versions'/(p.name+'.'+h[:16]+'.raw');copy(p,q)
    return dict(sha256=h,archived=q.relative_to(ROOT).as_posix())
def bind(p,h):
    if p.is_file() and sha(p)==h:return archive(p)
    candidates=[q for q in (DEST/'source_versions').glob(p.name+'.*.raw') if sha(q)==h]
    assert len(candidates)==1,'缺少实际使用的准确源码：'+str(p)
    return dict(sha256=h,archived=candidates[0].relative_to(ROOT).as_posix())
def main():
    target=DEST/'delivery_manifest.json';assert not target.exists();sources={}
    files=['scripts/prepareSU2CoolPropBuild.py','scripts/buildSU2CoolPropRuntime.py','scripts/finalizeSU2CoolPropBuild.py',
        'scripts/installSU2CoolPropRuntime.py','scripts/auditSU2CoolPropState.py','scripts/su2ThermalRuntime.py',
        'scripts/runEmbeddedSU2.py','scripts/runMatlabResearch.py','scripts/buildSU2ThermalRuntime.py',
        'scripts/recordCoolPropCapabilityDelivery.py','tests/prepareSU2ControlledPropertyProbe.py',
        'tests/prepareSU2CoolPropThermalProbe.py','tests/auditCoolPropNativeStates.py','tests/auditCoolPropWallFlux.py',
        'tests/testSU2CoolPropContracts.py','tests/reviewSU2CoolPropCapability.m','tests/auditEmbeddedSU2FlatPlate.py',
        'readSU2CoolPropStateAudit.m','auditSU2ResultConvergence.m','describeSU2Runtime.m','createSU2ResearchView.m',
        'readSU2SurfaceResult.m','identifySU2NativeCheckpoint.m','airThermalProperties.m','third_party/su2_thermal/RAHCoolPropState.cpp']
    for rel in files:
        p=ROOT/rel;sources[rel]=archive(p)
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'))
    originalBuild=read(OUT/'build_inputs.json')
    buildBindings={}
    for path,h in originalBuild.items():
        p=Path(path)
        if p.name in ['meson.build','generate_headers.py']:
            used=OUT/(p.name+'.used');assert sha(used)==h;buildBindings[path]=archive(used)
        elif p.suffix!='.exe':buildBindings[path]=bind(p,h)
    compilation=read(OUT/'compile_process.json');assert compilation['status']=='terminal' and compilation['exitCode']==0
    final=read(OUT/'build_finalization.json');assert final['actualCompileSourcesAndGeneratedPropertyTablesStable']
    for p in OUT.iterdir():
        if p.is_file() and p.suffix in ['.json','.log','.diff','.original','.used']:
            copy(p,DEST/'build_and_checks'/p.name)
    # 记录上游接口与数据定义，实际二进制/静态库不纳入Git。
    src=ROOT/'tools/runtime/su2-thermal-build/source';cp=src/'subprojects/CoolProp'
    upstream=['SU2_CFD/src/fluid/CCoolProp.cpp','SU2_CFD/include/fluid/CCoolProp.hpp',
        'SU2_CFD/include/solvers/CFVMFlowSolverBase.inl','SU2_CFD/src/solvers/CNSSolver.cpp']
    for rel in upstream:sources['upstream/'+rel]=archive(src/rel)
    for p in [cp/'dev/fluids/Air.json',cp/'LICENSE',cp/'CoolPropBibTeXLibrary.bib']:
        copy(p,DEST/'upstream'/p.name)
    runtime=ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-v1'
    for p in (runtime/'licenses').rglob('*'):
        if p.is_file():copy(p,DEST/'licenses'/p.relative_to(runtime/'licenses'))
    for p in runtime.iterdir():
        if p.is_file() and p.name!='runtime.json':copy(p,DEST/'runtime_sources'/p.name)
    preparations=ROOT/'tools/runtime/su2-coolprop-build'
    for p in preparations.glob('*.json'):copy(p,DEST/'dependencies'/p.name)
    engines={}
    for name in ['capability01','thermal20']:
        folder=OUT/name;r=read(folder/'engine_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut'] and all(r['inputsStable'].values())
        used={}
        for path,h in r['inputSHA256'].items():
            p=Path(path)
            if p.suffix in ['.py','.md']:used[path]=bind(p,h)
            else:assert sha(p)==h,path
        for name2,h in r['outputSHA256'].items():assert sha(folder/name2)==h
        for p in folder.iterdir():
            if p.is_file() and p.name!='mesh.su2':copy(p,DEST/'cases'/name/p.name)
        for p in (folder/'thermodynamics').iterdir():copy(p,DEST/'cases'/name/'thermodynamics'/p.name)
        engines[name]=dict(exitCode=0,elapsedSeconds=r['elapsedSeconds'],usedSources=used,sourceStableAtRun=True)
    copy(OUT/'thermal20/mesh.su2',DEST/'inputs/mesh.su2')
    for folder in ['properties','contracts']:
        for p in (OUT/folder).iterdir():
            if p.is_file():copy(p,DEST/folder/p.name)
    matlab=read(OUT/'review_coolprop01_process.json');assert matlab['status']=='terminal' and matlab['exitCode']==0 and all(matlab['sourceStableAfter'].values())
    matlabUsed={}
    for rel in [p for p in files if p.endswith('.m')]:matlabUsed[rel]=bind(ROOT/rel,matlab['sourceSHA256'][rel])
    assert not read(OUT/'matlab_review.json')['audit']['checkpointResidualPassed']
    localFiles=[p for p in OUT.rglob('*') if p.is_file()];objects={}
    for p in localFiles:
        st=p.stat();objects[(st.st_dev,st.st_ino)]=st.st_size
    storage=dict(outputPathBytes=sum(p.stat().st_size for p in localFiles),outputUniqueFileBytes=sum(objects.values()),
        runtimeFileBytes=sum(p.stat().st_size for p in runtime.rglob('*') if p.is_file()),
        newBuildFileBytes=sum(p.stat().st_size for p in preparations.rglob('*') if p.is_file()),
        actualNTFSAllocationMeasured=False,sharedOldBuildTreeIncluded=False,
        note='file lengths only; no old results/cache deletion; downloaded dependencies also appear in reused source tree',v117BlockedCleanupRetried=False)
    (DEST/'storage.json').write_text(json.dumps(storage,indent=2)+'\n')
    report=dict(status='tested-compressible-property-and-native-nonzero-wall-heat-interface',physicsQualified=False,newQualifiedArea=0,
        baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip(),
        engines=engines,matlab=dict(exitCode=0,sourceStable=True,usedSources=matlabUsed),
        sources=sources,actualBuildSources=buildBindings,storage=storage,
        originalBuildPostAuditFailurePreserved=True,rawPTIdentityFailurePreserved=True,
        newThreeDimensionalCoolPropCase=False,newMainUIOrBrowserTest=False,
        buildCacheExclusionFix='syntax checked and equivalent to tested post-build audit; full build not repeated',
        localOnly=['runtime executables and static libraries','shared build source tree','old full rocket fields and history store'],
        evidence={p.relative_to(ROOT).as_posix():sha(p) for p in DEST.rglob('*') if p.is_file()})
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(evidenceFiles=len(report['evidence']),sourceFiles=len(sources),**storage),indent=2))
if __name__=='__main__':main()
