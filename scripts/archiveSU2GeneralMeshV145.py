"""归档v145实际来源、失败与成功证据；只读数值结果，不重跑、不打工程包。"""
from pathlib import Path
import json,shutil,subprocess
from prepareSU2CheckpointContinuation import ROOT,read,sha

def copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True);assert not dst.exists()
    shutil.copyfile(str(src),str(dst));assert sha(src)==sha(dst)

def main():
    base=ROOT/'outputs/v145_validation';out=ROOT/'docs/validation/v145';target=out/'delivery_manifest.json';assert not target.exists()
    case=Path(read(base/'main01/engine_folder.json')['folder']);r=read(case/'engine_process.json')
    assert r['status']=='terminal' and r['exitCode']==0 and all(r['inputsStable'].values())
    assert read(base/'independent_review.json')['passed'] and read(base/'main02/report.json')['webgl']['webglError']==0
    for name,h in r['outputSHA256'].items():assert sha(case/name)==h
    for p,h in r['inputSHA256'].items():assert sha(Path(p))==h
    # 实际失败的准确源码已在修改前保存；任何找不到对应字节的输入都拒绝归档。
    prior={sha(p):p for p in out.rglob('*') if p.is_file() and '_sources' in str(p)}
    versions={};processes={}
    for name in ['fixture01','fixture02','main01','main02']:
        q=read(base/(name+'_process.json'));assert q['status']=='terminal' and all(q['sourceStableAfter'].values())
        processes[name]=dict(exitCode=q['exitCode'],sourceStable=True)
        for n,h in q['sourceSHA256'].items():
            p=ROOT/n
            if sha(p)!=h:
                assert h in prior,(name,n);versions[name+'/'+n]=str(prior[h].relative_to(ROOT))
    own=['auditSU2ResultConvergence.m','describeSU2Convergence.m','readSU2HybridVolume.m','rocketAeroHeatWorkbench.m',
        'scripts/runSU2WorkbenchPipeline.py','scripts/su2ThermalRuntime.py','scripts/prepareSU2GeneralHybrid.py',
        'scripts/prepareSU2WallDistanceState.py','scripts/prepareSU2NearWallFixture.py','scripts/reviewSU2GeneralMeshRun.py',
        'scripts/diagnoseSU2GeneralStartup.py','scripts/archiveSU2GeneralMeshV145.py',
        'tests/testSU2GeneralInitialization.m','tests/testSU2GeneralMeshWorkbench.m','tests/testSU2WallDistanceState.py']
    dependencies=['scripts/su2SSTInterface.py','scripts/prepareSU2SSTCase.py','scripts/runEmbeddedSU2.py',
        'scripts/auditSU2SSTState.py','scripts/reviewSU2SSTLocalDiffusion.py','scripts/runMatlabResearch.py',
        'embeddedSU2Job.m','readSU2SSTStateAudit.m','writeSU2RocketCase.m','createSU2ResearchView.m']
    for n in own+dependencies:copy(ROOT/n,out/'source_versions'/n)
    # 并发工作树中的既存UI修改不混入本次交付；仅保留本地准确副本和公开哈希记录。
    external={}
    for n in ['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py']:
        source=ROOT/n;dst=base/'uncommitted_ui_sources'/n;copy(source,dst)
        blob=subprocess.check_output(['git','show','HEAD:'+n],cwd=str(ROOT))
        import hashlib
        external[n]=dict(workingSHA256=sha(source),baseGitSHA256=hashlib.sha256(blob).hexdigest(),
            retainedLocalPath=str(dst.relative_to(ROOT)),includedInThisCommit=False)
    for folder,label in [(base,'phase'),(case,'actual')]:
        for p in folder.rglob('*'):
            if not p.is_file() or 'uncommitted_ui_sources' in p.parts:continue
            relative=p.relative_to(folder)
            small=p.suffix in ['.json','.log','.cfg'] or p.name in ['history.csv','section.csv']
            fixture='fixture01' in p.parts or 'fixture02' in p.parts or 'near_wall_fixture' in p.parts
            if fixture and p.suffix in ['.mat','.csv','.su2','.npz']:small=True
            if small and p.stat().st_size<=2*1024**2:copy(p,out/'evidence'/label/relative)
    files={str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=sha(p)) for folder in [base,case] for p in folder.rglob('*') if p.is_file()}
    evidence={str(p.relative_to(out)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.rglob('*') if p.is_file()}
    length=sum(v['bytes'] for v in files.values());assert length<1500*1024**2
    result=dict(passed=True,ownedSourcesSHA256={n:sha(ROOT/n) for n in own},processes=processes,
        preservedEarlierVersions=versions,existingUncommittedDependencies=external,outputFiles=files,evidenceFiles=evidence,
        outputPathFileLengthBytes=length,evidenceFileLengthBytes=sum(v['bytes'] for v in evidence.values()),
        sizeMeaning='sum of file lengths including retained hard-link paths; not allocated space or newly added unique data',
        actualNativeRuns=1,smallNativeRuns=2,newGeometryRuns=0,largeResultsUploaded=False,physicsQualified=False,newQualifiedArea=0)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['passed','outputPathFileLengthBytes','evidenceFileLengthBytes','actualNativeRuns','smallNativeRuns']}))
if __name__=='__main__':main()
