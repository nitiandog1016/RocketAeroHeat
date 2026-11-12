"""模块职责：区分实际失败源与未执行记录器，保护既有资产并归档原面科研证据。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runPhaseEulerRocketV230 import ROOT,BASE,DOC,PREV,OLD_BUILD,HELPER,CLOSED,read,save,sha

OWN=['native/su2PhaseAwareAir/rocketEulerRun.cpp','native/su2PhaseAwareAir/rocketEulerFailure.hpp',
 'scripts/runPhaseEulerRocketV230.py','scripts/auditPhaseEulerRocketV230.py',
 'scripts/finalizePhaseEulerRocketV230.py','tests/testPhaseEulerViewerV230.js',
 'docs/HYBRID_PHASE_EULER_ROCKET_V230.md','docs/validation/v230/PLAN.md',
 'docs/validation/v230/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def main():
    started=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：只读核对并发工程候选、旧构建对象及原始几何/物性资产。
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
    mesh=ROOT/'outputs/v194_validation/rocket01/mesh.su2'
    assert sha(mesh)==sha(BASE/'rocket01/mesh.su2')=='faff7bf901d330eaf53a4616014620281f54e58e8137d1a0e4943ab14804f0e7'
    assets[str(mesh.relative_to(ROOT)).replace('\\','/')]=sha(mesh)
    # 模块二：失败保持失败；核对实际执行源，不能用后加记录器替代原版本。
    records={p.parent.name[8:]:read(p) for p in BASE.glob('process_*/report.json')}
    failures={'prepare','audit','run'}
    assert len(records)==10
    for stage,r in records.items():
        assert r['passed']==(stage not in failures),(stage,r['passed'])
        assert r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['stopReason'] is None,stage
    assert records['run']['exitCode']==3221225501
    for stage in ('build01','build02'):
        for n,h in read(BASE/stage/'dependencies.json').items():
            p=BASE/'prepare02/executed_helper.cpp' if stage=='build01' and Path(n)==HELPER else Path(n)
            assert sha(p)==h,n
    for p in ('rocketEulerRun.cpp','rocketEulerFailure.hpp'):
        assert sha(BASE/'observer01'/p)==sha(HELPER.parent/p)
    assert records['run']['inputSHA256'][str(HELPER)]==sha(BASE/'prepare02/executed_helper.cpp')
    assert records['build_observer']['inputSHA256'][str(HELPER)]==sha(HELPER)
    assert read(BASE/'crash01/report.json')['executableSHA256']==sha(BASE/'build01/rocket.exe')
    for group in read(PREV/'build02/members.json').values():
        for r in group:assert sha(r['path'])==r['sha256'],r['path']
    audit=read(BASE/'audit02/report.json');ui=read(BASE/'workbench01/report.json')
    assert audit['passed'] and not audit['scientificAcceptance'] and audit['completedUpdates']==30
    assert audit['savedFullStateUpdates']==[0,1] and not audit['underlyingFailureCauseKnown']
    assert not audit['allFiveBelowMinus8'] and not audit['fullMomentumBalanceAccepted']
    assert not audit['originalAbsoluteSumCriterionPassed'] and audit['sumConsistentWithinRoundoff']
    assert ui['passed'] and ui['fields']==11 and ui['originalFaces']==19752
    assert ui['exportSHA256']==sha(BASE/'workbench01/rocket_euler_research.json')
    assert ui['imageSHA256']==sha(BASE/'workbench01/rocket_euler_research.png')
    assert sha(ROOT/'web/viewer.js')==records['export']['inputSHA256'][str(ROOT/'web/viewer.js')]
    assert CLOSED
    refusal=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseEulerRocketV230.py'),'run'],
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert refusal.returncode!=0 and b'closed' in refusal.stdout
    # 模块三：小证据与原面显示公开；完整体场、二进制和机器WER留本地。
    out.mkdir();evidence.mkdir()
    for stage in ('prepare01','prepare02','observer01','crash01','workbench01'):
        shutil.copytree(BASE/stage,evidence/stage)
    for folder in BASE.glob('process_*'):shutil.copytree(folder,evidence/folder.name)
    for stage in ('build01','build02'):
        target=evidence/stage;target.mkdir()
        for p in (BASE/stage).iterdir():
            if p.suffix in ('.json','.d'):shutil.copyfile(p,target/p.name)
    for stage,names in {'rocket01':['case.cfg','contract.json','history.csv','physical_history.csv'],
                        'audit02':['report.json','original_faces_000000.csv','original_faces_000001.csv']}.items():
        target=evidence/stage;target.mkdir()
        for name in names:shutil.copyfile(BASE/stage/name,target/name)
    local={str(p.relative_to(BASE)).replace('\\','/'):dict(bytes=p.stat().st_size,sha256=sha(p))
      for p in [BASE/'rocket01/state_000000.csv',BASE/'rocket01/state_000001.csv',BASE/'rocket01/mesh.su2',
                BASE/'audit02/geometry.npz',BASE/'build01/rocket.exe',BASE/'build02/rocket.exe']}
    save(evidence/'local_asset_index.json',dict(base='outputs/v230_validation',files=local,
        completeVolumeFieldsPublished=False,binariesPublished=False,fullMachineWERPublished=False))
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        target=evidence/'final_source'/(n if n!='.gitattributes' else '.gitattributes.snapshot')
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    totals={kind:sum(r['elapsedSeconds'] for r in records.values() if r['budgetKind']==kind)
            for kind in ('offline','build','cfd')}
    elapsed=time.monotonic()-started
    assert totals['offline']+elapsed<120 and totals['build']<180 and totals['cfd']<600
    output_bytes=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert output_bytes<512*1024**2
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,scientificAcceptance=False,
      actualNativeStartSolverExecuted=True,completedUpdates=30,savedFullStateUpdates=[0,1],normalCompletion=False,
      originalRunExitCode=records['run']['exitCode'],terminateEstablished=True,underlyingFailureCauseKnown=False,
      observationRecorderCompiled=True,observationRecorderExecuted=False,newCFDReplays=0,
      actualRunSourceSHA256=sha(BASE/'prepare02/executed_helper.cpp'),
      unexecutedRecorderSourceSHA256=sha(HELPER),unexecutedRecorderHeaderSHA256=sha(HELPER.parent/'rocketEulerFailure.hpp'),
      newHeatAreaM2=0,newQualifiedAreaM2=0,originalFaces=19752,originalAreaM2=1.9646292565715553,
      exportedFields=11,frontendComputedValuesVerified=177768,heatFluxAndSolidTemperatureUnknown=True,
      protectedCount=len(protected),engineeringArtifactsUnchanged=len(engineering['files']),
      engineeringCandidateCount=len(candidates['files']),oldObjectsUnchanged=len(old),oldDependenciesUnchanged=len(deps),
      assets=assets,processSeconds=totals,finalizationSeconds=elapsed,outputBytesBeforeReport=output_bytes,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      actualCFDSampledGroupRSSBytes=records['run']['peakObservedGroupRSSBytes'],allGroupsEmpty=True,
      unmeteredTargetedReadOnlyInspectionPresent=True,oldFailuresPreserved=True,MATLABExecuted=False,
      sharedWebModified=False,browserInteractionTested=False,uniqueAssetsDeleted=False,physicalRunCLIClosed=True)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v230/evidence_manifest.json']+
        ['docs/validation/v230/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,
      processSeconds=totals,originalRunFailed=True,scientificAcceptance=False,recorderExecuted=False)))

if __name__=='__main__':main()
