"""模块职责：归档800更新/原面展示证据与容量失败，保护旧资产，不把压力稳定当热资格。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runPhaseEulerDevelopmentV232 import ROOT,BASE,DOC,CASE,BUILD,PREV,BASELINE,HEADER,CLOSED,read,save,sha,OLD_BUILD

OWN=['native/su2PhaseAwareAir/EulerDevelopmentObservationV232.hpp','scripts/runPhaseEulerDevelopmentV232.py',
 'scripts/auditPhaseEulerDevelopmentV232.py','scripts/finalizePhaseEulerDevelopmentV232.py','tests/testPhaseEulerViewerV232.js',
 'docs/HYBRID_PHASE_EULER_DEVELOPMENT_V232.md','docs/validation/v232/PLAN.md','docs/validation/v232/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def main():
    started=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：原网格、气体、旧构建和并发候选只读核对。
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256'],r['file']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(OLD_BUILD/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    assets=read(PREV/'final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    prior_dependencies=read(ROOT/'docs/validation/v231/evidence/actual_dependencies.json')
    for n,h in prior_dependencies.items():assert sha(n)==h,n
    for rows in read(ROOT/'outputs/v229_validation/build02/members.json').values():
        for r in rows:assert sha(r['path'])==r['sha256'],r['path']
    previous=read(PREV/'final01/report.json');assert sha(PREV/'recovery_build02/primitive.obj')==previous['primitiveObjectSHA256']
    assert sha(CASE/'mesh.su2')==sha(BASELINE/'rocket01/mesh.su2') and sha(CASE/'case.cfg')==sha(BASELINE/'rocket01/case.cfg')
    records={p.parent.name[8:]:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==8
    for name,r in records.items():
        assert r['passed']==(name!='export'),name
        assert r['inputsStable'] and not r['liveGroupAfterCleanup']
        assert r['stopReason']==('output' if name=='export' else None),name
    assert records['export']['exitCode']==15 and records['run']['exitCode']==0
    audit=read(BASE/'audit01/report.json');ui=read(BASE/'workbench01/report.json');plot=read(BASE/'plot02/report.json')
    region=read(BASE/'regional01/report.json')
    assert audit['passed'] and audit['completedUpdates']==800 and audit['prefix30ByteIdentical']
    assert len(audit['savedUpdates'])==18 and len(audit['savedFrameChecks'])==16
    assert audit['last50PressureStable'] and not audit['allFiveBelowMinus8'] and not audit['scientificAcceptance']
    assert audit['rejectionLastPass']==8 and not audit['momentumSameStageVerified']
    assert ui['passed'] and ui['fields']==13 and plot['passed'] and region['originalTargetFaces']==3128
    assert ui['exportSHA256']==plot['JSONSHA256']==sha(BASE/'workbench01/rocket_euler_research.json')
    assert ui['imageSHA256']==sha(BASE/'workbench01/rocket_euler_comparison.png')
    assert plot['imageSHA256']==sha(BASE/'plot02/rocket_euler_comparison.png')
    for n,h in region['sourceSHA256'].items():assert sha(ROOT/n)==h,n
    assert sha(ROOT/'web/viewer.js')==records['export2']['inputSHA256'][str(ROOT/'web/viewer.js')]
    bindings=read(BASE/'prepare01/dependencies.json')
    for n,h in bindings.items():assert sha(n)==h,n
    assert sha(HEADER)==sha(BASE/'prepare01'/HEADER.name)
    actual={};text=(BUILD/'rocket.d').read_text(encoding='utf8').replace('\\\n',' ')
    for name in text[text.index(': ')+2:].split():
        p=Path(name);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p);assert bindings[str(p)]==h,str(p);actual[str(p)]=h
    assert records['run']['inputSHA256'][str(BUILD/'rocket.exe')]==sha(BUILD/'rocket.exe')
    assert CLOSED
    refused=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseEulerDevelopmentV232.py'),'run'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert refused.returncode!=0 and b'closed' in refused.stdout
    # 模块二：三阶段原面、日志、源码公开；18体场及原生重启文件仍留E盘。
    out.mkdir();evidence.mkdir()
    for folder in BASE.glob('process_*'):shutil.copytree(folder,evidence/folder.name)
    for stage in ('prepare01','audit01','regional01','workbench01','plot02'):shutil.copytree(BASE/stage,evidence/stage)
    target=evidence/'build01';target.mkdir()
    for p in BUILD.iterdir():
        if p.suffix in ('.json','.d'):shutil.copyfile(p,target/p.name)
    target=evidence/'rocket01';target.mkdir()
    for n in ('case.cfg','contract.json','native_completion.json','prefix_verified.json','native_rejections.csv','history.csv','physical_history.csv'):
        shutil.copyfile(CASE/n,target/n)
    local={p.relative_to(BASE).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in CASE.iterdir()
      if p.is_file() and (p.name.startswith(('state_','restart_flow')) or p.name=='mesh.su2')}
    local['build01/rocket.exe']=dict(bytes=(BUILD/'rocket.exe').stat().st_size,sha256=sha(BUILD/'rocket.exe'))
    save(evidence/'local_asset_index.json',dict(base='outputs/v232_validation',files=local,volumeAndNativeRestartPublished=False,binariesPublished=False))
    save(evidence/'actual_dependencies.json',actual)
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        target=evidence/'final_source'/(n if n!='.gitattributes' else '.gitattributes.snapshot')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    totals={k:sum(r['elapsedSeconds'] for r in records.values() if r['budgetKind']==k) for k in ('offline','build','cfd','audit','ui')}
    elapsed=time.monotonic()-started
    assert totals['offline']+elapsed<120 and totals['build']<60 and totals['cfd']<600 and totals['audit']+totals['ui']<120
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<512*1024**2
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,scientificAcceptance=False,normalCompletion=True,
      completedUpdates=800,completeSavedFrames=18,newSavedFramesDomainVerified=16,originalPrefix30ByteIdentical=True,
      nativeDomainRecoveryObserved=True,lastNativePreprocessingRejectedPoints=8,actualRejectedNodeIDsKnown=False,
      lastLogPhysicalRMS=audit['lastLogPhysicalRMS'],last50PressureStable=True,
      last50PressureAreaRelativeChange=audit['last50PressureAreaRelativeChange'],momentumComparisonCrossStageOnly=True,
      steadyNetBoundaryScaled=audit['steadyBoundaryNetScaled'],spatialFarfieldQualified=False,newHeatAreaM2=0,newQualifiedAreaM2=0,
      originalFaces=19752,originalAreaM2=1.9646292565715553,exportedFields=13,frontendComputedValuesVerified=217272,
      pressureFigureScale='logarithmic',figureManuallyViewed=True,originalPlotPreserved=True,
      protectedCount=len(protected),engineeringArtifactsUnchanged=len(engineering['files']),engineeringCandidateCount=len(candidates['files']),
      oldObjectsUnchanged=len(old),oldDependenciesUnchanged=len(deps),priorRepairedDependenciesUnchanged=len(prior_dependencies),assets=assets,
      actualDriverDependencies=len(actual),primitiveObjectUnchanged=True,processSeconds=totals,finalizationSeconds=elapsed,
      outputBytesBeforeReport=size,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      actualCFDSampledGroupRSSBytes=records['run']['peakObservedGroupRSSBytes'],allGroupsEmpty=True,
      initialUIOutputGuardDoubleCountFailurePreserved=True,finalOutputGuardCountsWholeBaseOnce=True,
      targetedReadOnlyInspectionNotFullyMetered=True,MATLABExecuted=False,sharedWebModified=False,browserInteractionTested=False,
      uniqueAssetsDeleted=False,newMeshGenerated=False,physicalRunCLIClosed=True)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v232/evidence_manifest.json']+['docs/validation/v232/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,processSeconds=totals,
      completedUpdates=800,pressureStable=True,scientificAcceptance=False,bytes=size)))

if __name__=='__main__':main()
