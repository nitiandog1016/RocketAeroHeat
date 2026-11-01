"""归档通量小检查、完整域对照与工作台；仅复制小证据，不上传大场。"""
from pathlib import Path
import hashlib,json,shutil,subprocess
from prepareSU2SSTFluxStudy import ROOT,BASE,read,sha,write

def copy(source,target):
    target.parent.mkdir(parents=True,exist_ok=True);assert not target.exists()
    shutil.copyfile(str(source),str(target));assert sha(source)==sha(target)

def main():
    out=ROOT/'docs/validation/v148';target=out/'delivery_manifest.json';assert not target.exists()
    assert read(BASE/'flux_preflight02/report.json')['passed'] and read(BASE/'export_review.json')['passed']
    comparison=read(BASE/'flux_comparison.json');processes={}
    # 模块一：完整实际场的输入/输出逐字绑定，保存失败筛选而不改默认。
    for name in ['small_up2','control','up2']:
        folder=BASE/name;r=read(folder/'engine_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert read(folder/'sst_state_audit.json')['interfacePassed']
        for p,h in r['inputSHA256'].items():assert sha(Path(p))==h,(name,p)
        for p,h in r['outputSHA256'].items():assert sha(folder/p)==h,(name,p)
        processes[name]=dict(exitCode=0,seconds=r['elapsedSeconds'],inputsStable=True)
    ui=read(BASE/'workbench01_process.json');assert ui['status']=='terminal' and ui['exitCode']==0 and all(ui['sourceStableAfter'].values())
    for p,h in ui['sourceSHA256'].items():assert sha(ROOT/p)==h,p
    # 模块二：准确最终执行源和既存未提交UI依赖分开保存。
    own=['auditSU2ResultConvergence.m','describeSU2Runtime.m','readSU2VerifiedSSTStudy.m','scripts/su2ThermalRuntime.py',
        'scripts/verifySU2AUSMUp2Flux.py','scripts/reviewSU2AUSMUp2Flux.py','scripts/prepareSU2SSTFluxStudy.py',
        'scripts/compareSU2SSTFluxStudy.py','scripts/reviewSU2SSTFluxExports.py','scripts/archiveSU2SSTFluxV148.py',
        'tests/testSU2SSTFluxWorkbench.m']
    deps=['scripts/prepareSU2SSTCheckpoint.py','scripts/prepareSU2CheckpointContinuation.py','scripts/prepareSU2SSTCase.py',
        'scripts/su2SSTInterface.py','scripts/runEmbeddedSU2.py','scripts/auditSU2SSTState.py','scripts/runMatlabResearch.py',
        'scripts/compareSU2SSTStartupStudy.py','readSU2SSTStateAudit.m','createSU2ResearchView.m','rocketAeroHeatWorkbench.m',
        'describeSU2Convergence.m']
    for p in own+deps:copy(ROOT/p,out/'source_versions'/p)
    external={}
    for p in ['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py']:
        dst=BASE/'uncommitted_ui_sources'/p;copy(ROOT/p,dst)
        blob=subprocess.check_output(['git','show','HEAD:'+p],cwd=str(ROOT))
        external[p]=dict(workingSHA256=sha(ROOT/p),baseGitSHA256=hashlib.sha256(blob).hexdigest(),
            retainedLocalPath=str(dst.relative_to(ROOT)),includedInThisCommit=False)
    # 模块三：证据小文件按字节复制；当前源和日志均已终态。
    for p in list(BASE.rglob('*')):
        if not p.is_file() or 'uncommitted_ui_sources' in p.parts:continue
        if p.stat().st_size<2*1024**2 and (p.suffix in ['.json','.log','.cfg','.cpp','.hpp','.txt','.stdout','.stderr'] or p.name in ['history.csv','section.csv','native_flux.csv','run_bounded.py']):
            copy(p,out/'evidence'/p.relative_to(BASE))
    files={str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*') if p.is_file()}
    unique={}
    for p in BASE.rglob('*'):
        if p.is_file():s=p.stat();unique[(s.st_dev,s.st_ino)]=s.st_size
    inherited={}
    for folder in [ROOT/'outputs/v145_validation/fixture02/case',ROOT/'outputs/v146_validation/control']:
        for name in ['mesh.su2','restart_flow.csv']:
            s=(folder/name).stat();inherited[(s.st_dev,s.st_ino)]=s.st_size
    unique_new=sum(size for key,size in unique.items() if key not in inherited);length=sum(v['bytes'] for v in files.values())
    assert length<2200*1024**2 and unique_new<1500*1024**2
    evidence={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file()}
    result=dict(passed=True,processes=processes,ownedSourcesSHA256={p:sha(ROOT/p) for p in own},
        outputFiles=files,evidenceFiles=evidence,existingUncommittedDependencies=external,
        outputPathFileLengthBytes=length,uniqueCurrentStageFileLengthExcludingInheritedBytes=unique_new,
        sizeMeaning='file lengths, not allocated space; inherited hardlinks excluded only from unique-current-stage count',
        nativeFullDomainRuns=2,nativeSmallRuns=1,fragmentCompileAttempts=2,fragmentExecutions=1,
        newRuntimePackages=0,sharedBuildCacheModified=False,newGeometryRuns=0,
        candidateEligibleForFurtherStudy=comparison['candidateEligibleForFurtherStudy'],
        defaultChanged=False,physicsQualified=False,newQualifiedArea=0,largeResultsUploaded=False,noUniqueDataDeleted=True)
    write(target,result);print(json.dumps({k:result[k] for k in ['passed','nativeFullDomainRuns','candidateEligibleForFurtherStudy','outputPathFileLengthBytes','uniqueCurrentStageFileLengthExcludingInheritedBytes']}))

if __name__=='__main__':main()
