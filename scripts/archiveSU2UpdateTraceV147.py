"""归档实际原生更新、失败和工作台来源；大场/运行包仅本地，不包含他人的改动。"""
from pathlib import Path
import hashlib,json,shutil,subprocess,os
from prepareSU2UpdateTraceCase import ROOT,BASE,read,sha,write

def copy(source,target):
    target.parent.mkdir(parents=True,exist_ok=True);assert not target.exists()
    shutil.copyfile(str(source),str(target));assert sha(source)==sha(target)

def main():
    out=ROOT/'docs/validation/v147';target=out/'delivery_manifest.json';assert not target.exists()
    assert read(BASE/'small_comparison_v2.json')['passed']
    assert read(BASE/'native_update_analysis.json')['passed'] and read(BASE/'export_review.json')['passed']
    cases=['small_official','small_off','small_on','small_off_v2','small_on_v2','actual']
    processes={};historical=[]
    # 模块一：实际文件与各自回执对应；后续修正过的来源必须有逐字失败版本。
    preserved={sha(p):p for p in (out/'failure01').rglob('*') if p.is_file()}
    for name in cases:
        folder=BASE/name;r=read(folder/'engine_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
        assert read(folder/'sst_state_audit.json')['interfacePassed']
        for file,h in r['outputSHA256'].items():assert sha(folder/file)==h,(name,file)
        for path,h in r['inputSHA256'].items():
            p=Path(path)
            if sha(p)==h:continue
            assert h in preserved,(name,path)
            historical.append(dict(case=name,originalPath=path,sha256=h,retained=str(preserved[h].relative_to(ROOT))))
        processes[name]=dict(exitCode=r['exitCode'],seconds=r['elapsedSeconds'],inputsStable=True)
    # 模块二：两个MATLAB任务均已终态；准确记录旧字段拒绝及后续只读修正。
    for name,expected in [('workbench01',1),('workbench02',0)]:
        r=read(BASE/(name+'_process.json'))
        assert r['status']=='terminal' and r['exitCode']==expected and all(r['sourceStableAfter'].values())
    ui=read(BASE/'workbench02_process.json');post_test=[]
    for path,h in ui['sourceSHA256'].items():
        if sha(ROOT/path)==h:continue
        # 此导出复核器不由MATLAB执行；只在终态后改为读取成功的workbench02。
        assert path=='scripts/reviewSU2UpdateTraceExports.py'
        raw=(ROOT/path).read_bytes().replace(b'workbench02',b'workbench01')
        assert hashlib.sha256(raw).hexdigest()==h
        p=out/'source_versions/ui_locked'/path;p.parent.mkdir(parents=True,exist_ok=True);assert not p.exists();p.write_bytes(raw)
        post_test.append(dict(path=path,uiLockedSHA256=h,finalSHA256=sha(ROOT/path),reason='export reviewer destination corrected after MATLAB terminal; no solver or MATLAB source changed'))
    own=['auditSU2ResultConvergence.m','describeSU2Runtime.m','describeSU2Convergence.m','readSU2VerifiedSSTStudy.m',
        'readSU2UpdateTraceAudit.m','scripts/runEmbeddedSU2.py','scripts/su2ThermalRuntime.py',
        'scripts/prepareSU2UpdateTracePatch.py','scripts/buildSU2UpdateTraceRuntime.py','scripts/prepareSU2UpdateTraceCase.py',
        'scripts/reviewSU2UpdateTrace.py','scripts/diagnoseSU2NativeUpdate.py','scripts/reviewSU2UpdateTraceExports.py',
        'scripts/archiveSU2UpdateTraceV147.py','tests/testSU2UpdateTraceWorkbench.m',
        'third_party/su2_update_trace/RAHUpdateTrace.hpp','third_party/su2_update_trace/update_trace.diff']
    deps=['scripts/prepareSU2SSTCheckpoint.py','scripts/prepareSU2CheckpointContinuation.py','scripts/prepareSU2SSTCase.py',
        'scripts/su2SSTInterface.py','scripts/auditSU2SSTState.py','scripts/runMatlabResearch.py',
        'readSU2SSTStateAudit.m','createSU2ResearchView.m','rocketAeroHeatWorkbench.m']
    for name in own+deps:copy(ROOT/name,out/'source_versions/final'/name)
    external={}
    for name in ['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py']:
        dst=BASE/'uncommitted_ui_sources'/name;copy(ROOT/name,dst)
        blob=subprocess.check_output(['git','show','HEAD:'+name],cwd=str(ROOT))
        external[name]=dict(workingSHA256=sha(ROOT/name),baseGitSHA256=hashlib.sha256(blob).hexdigest(),
            retainedLocalPath=str(dst.relative_to(ROOT)),includedInThisCommit=False)
    # 模块三：独立构建/实际源恢复审核；不覆盖原包或把待验证清单改写成物理成功。
    builds={};runtime_lengths={}
    for num,flavor in [('01','rah-update-trace-v1'),('02','rah-update-trace-v2')]:
        folder=BASE/('build'+num);r=read(folder/'compile_process.json')
        assert r['status']=='terminal' and r['exitCode']==0 and not r['timedOut']
        assert all(read(folder/'source_restore.json').values())
        runtime=read(folder/'runtime_installed.json');dest=ROOT/'tools/runtime/su2'/('8.6.0-'+flavor)
        assert read(dest/'runtime.json')==runtime
        for item in runtime['files']:assert sha(dest/item['path'])==item['sha256']
        builds[num]=dict(seconds=r['elapsedSeconds'],exitCode=0,allObjectsRebuilt=True,sourceRestored=True)
        runtime_lengths[flavor]=sum(p.stat().st_size for p in dest.rglob('*') if p.is_file())
    # 只复制小型证据；大型CSV/MAT和二进制继续保留本地。
    for p in list(BASE.rglob('*')):
        if not p.is_file() or 'uncommitted_ui_sources' in p.parts:continue
        if p.stat().st_size<2*1024**2 and (p.suffix in ['.json','.log','.cfg','.used','.patched','.ninja','.hpp','.inl'] or p.name in ['history.csv','section.csv','run_bounded.py']):
            copy(p,out/'evidence'/p.relative_to(BASE))
    files={str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*') if p.is_file()}
    # 当前目录文件ID去重，并排除逐字共享的旧源；文件长度不冒充磁盘分配量。
    unique={}
    for p in BASE.rglob('*'):
        if p.is_file():s=p.stat();unique[(s.st_dev,s.st_ino)]=s.st_size
    inherited={}
    for folder in [ROOT/'outputs/v145_validation/fixture02/case',ROOT/'outputs/v146_validation/control']:
        for name in ['mesh.su2','restart_flow.csv']:
            s=(folder/name).stat();inherited[(s.st_dev,s.st_ino)]=s.st_size
    unique_new=sum(size for key,size in unique.items() if key not in inherited)
    length=sum(v['bytes'] for v in files.values());assert length<2200*1024**2 and unique_new<1500*1024**2
    evidence={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file()}
    result=dict(passed=True,processes=processes,builds=builds,ownedSourcesSHA256={n:sha(ROOT/n) for n in own},
        retainedHistoricalSources=historical,postTestReviewOnlyChanges=post_test,existingUncommittedDependencies=external,
        outputFiles=files,evidenceFiles=evidence,outputPathFileLengthBytes=length,
        uniqueCurrentStageFileLengthExcludingInheritedBytes=unique_new,runtimePathFileLengths=runtime_lengths,
        sizeMeaning='file lengths and file-ID deduplication; not allocated disk space; build cache reported separately',
        actualNativeRuns=1,smallNativeRuns=5,newGeometryRuns=0,physicsQualified=False,newQualifiedArea=0,
        largeResultsUploaded=False,noUniqueDataDeleted=True)
    write(target,result);print(json.dumps({k:result[k] for k in ['passed','actualNativeRuns','smallNativeRuns','outputPathFileLengthBytes','uniqueCurrentStageFileLengthExcludingInheritedBytes','runtimePathFileLengths']}))

if __name__=='__main__':main()
