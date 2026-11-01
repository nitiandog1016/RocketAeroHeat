"""归档v146同源启动对照及科研工作台证据；大场仅本地，保留准确执行源。"""
from pathlib import Path
import hashlib,json,shutil,subprocess
from prepareSU2SSTStartupStudy import ROOT,BASE,read,sha,write

def copy(source,target):
    target.parent.mkdir(parents=True,exist_ok=True);assert not target.exists()
    shutil.copyfile(str(source),str(target));assert sha(source)==sha(target)

def main():
    out=ROOT/'docs/validation/v146';target=out/'delivery_manifest.json';assert not target.exists()
    comparison=read(BASE/'startup_comparison.json');assert read(BASE/'export_review.json')['passed']
    variants=['control','first_order']+(['restore_second_order'] if comparison['restoreAllowed'] else [])
    processes={}
    for name in variants:
        folder=BASE/name;r=read(folder/'engine_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert read(folder/'sst_state_audit.json')['interfacePassed']
        for file,h in r['outputSHA256'].items():assert sha(folder/file)==h
        for path,h in r['inputSHA256'].items():assert sha(Path(path))==h
        processes[name]=dict(exitCode=r['exitCode'],seconds=r['elapsedSeconds'],inputsStable=True)
    # MATLAB锁定的全部源必须仍然一致，未提交的并行UI依赖只在本地保留准确副本。
    ui=read(BASE/'workbench01_process.json');assert ui['status']=='terminal' and ui['exitCode']==0 and all(ui['sourceStableAfter'].values())
    for path,h in ui['sourceSHA256'].items():assert sha(ROOT/path)==h,path
    own=['auditSU2ResultConvergence.m','readSU2VerifiedSSTStudy.m','describeSU2Runtime.m',
        'scripts/su2ThermalRuntime.py','scripts/prepareSU2SSTStartupStudy.py','scripts/compareSU2SSTStartupStudy.py',
        'scripts/reviewSU2SSTStartupExports.py','scripts/archiveSU2SSTStartupV146.py',
        'tests/testSU2SSTStartupStudy.py','tests/testSU2SSTStartupWorkbench.m']
    dependencies=['scripts/prepareSU2SSTCheckpoint.py','scripts/prepareSU2CheckpointContinuation.py',
        'scripts/prepareSU2SSTCase.py','scripts/su2SSTInterface.py','scripts/runEmbeddedSU2.py',
        'scripts/auditSU2SSTState.py','scripts/runMatlabResearch.py','readSU2SSTStateAudit.m',
        'createSU2ResearchView.m','rocketAeroHeatWorkbench.m','describeSU2Convergence.m']
    for name in own+dependencies:copy(ROOT/name,out/'source_versions'/name)
    external={}
    for name in ['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py']:
        dst=BASE/'uncommitted_ui_sources'/name;copy(ROOT/name,dst)
        blob=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(ROOT))
        external[name]=dict(workingSHA256=sha(ROOT/name),baseGitSHA256=hashlib.sha256(blob).hexdigest(),
            retainedLocalPath=str(dst.relative_to(ROOT)),includedInThisCommit=False)
    for p in list(BASE.rglob('*')):
        if not p.is_file() or 'uncommitted_ui_sources' in p.parts:continue
        if (p.suffix in ['.json','.log','.cfg'] or p.name in ['history.csv','section.csv','run_bounded.py']) and p.stat().st_size<2*1024**2:
            copy(p,out/'evidence'/p.relative_to(BASE))
    files={str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*') if p.is_file()}
    evidence={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file()}
    length=sum(x['bytes'] for x in files.values());assert length<2400*1024**2
    result=dict(passed=True,processes=processes,ownedSourcesSHA256={n:sha(ROOT/n) for n in own},
        outputFiles=files,evidenceFiles=evidence,existingUncommittedDependencies=external,
        outputPathFileLengthBytes=length,sizeMeaning='file lengths including hardlink aliases; not allocated space',
        actualNativeRuns=len(variants),newGeometryRuns=0,physicsQualified=False,newQualifiedArea=0,largeResultsUploaded=False)
    write(target,result);print(json.dumps({k:result[k] for k in ['passed','actualNativeRuns','outputPathFileLengthBytes']}))
if __name__=='__main__':main()
