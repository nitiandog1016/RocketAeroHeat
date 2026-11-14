"""模块职责：归档原生入口/完整控制/真实失败，核对保护与来源，不执行新数值任务。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v235_validation';DOC=ROOT/'docs/validation/v235'
OWN=['native/su2PhaseAwareAir/SingleDomainHooksV235.hpp','native/su2PhaseAwareAir/singleDomainRunV235.cpp',
 'scripts/prepareSingleDomainV235.py','scripts/runSingleDomainV235.py','scripts/auditSingleDomainV235.py',
 'scripts/singleDomainFrame.py','scripts/finalizeSingleDomainV235.py','docs/SINGLE_DOMAIN_NATIVE_ENTRY_V235.md',
 'docs/validation/v235/PLAN.md','docs/validation/v235/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    eng=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in eng['files']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    assets=read(ROOT/'docs/validation/v222/evidence/final_report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==11
    failed={'process_stage':1,'process_off':86,'process_observe':86};historical=[];totals={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['stopReason'] is None,name
        assert r['exitCode']==failed.get(name,0) and r['passed']==(name not in failed),name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            snapshot=BASE/name/'execution_sources'/(h+Path(p).suffix);assert sha(snapshot)==h,(name,p)
            historical.append(dict(process=name,path=p,executedSHA256=h,currentSHA256=sha(p)))
    for r in read(BASE/'source01/bindings.json'):assert sha(r['source'])==r['sha256'] and sha(r['snapshot'])==r['sha256']
    for n,h in read(BASE/'stage01/dependencies.json').items():assert sha(n)==h,n
    for rs in read(BASE/'build01/members.json').values():
        for r in rs:assert sha(r['path'])==r['sha256']
    a=read(BASE/'audit01/report.json');assert a['passed'] and a['nativeProductionLoopExercised'] and a['structuralRejections']==8
    assert a['fullRocketFlowUpdates']==0 and not a['newRejectionCallbackDynamicallyExercised']
    assert not read(BASE/'audit01/real_admission_decision.json')['fullRocketAdmitted']
    closed={}
    for n in OWN:
        if n.endswith('.py'):ast.parse((ROOT/n).read_text(encoding='utf8'))
    for n in OWN[2:5]:
        p=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
        text=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in text
        closed[n]=text.strip()
    out.mkdir();evidence.mkdir()
    for p in BASE.glob('process_*'):shutil.copytree(p,evidence/p.name)
    for name in ('source01','stage01','off01','observe01','off02','observe02','audit01'):
        shutil.copytree(BASE/name,evidence/name)
    target=evidence/'build01';target.mkdir()
    for p in (BASE/'build01').iterdir():
        if p.suffix in ('.json','.d'):shutil.copyfile(p,target/p.name)
    for p in BASE.glob('control_contract*.json'):shutil.copyfile(p,evidence/p.name)
    local={p.relative_to(BASE).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in (BASE/'build01').iterdir() if p.suffix in ('.exe','.a','.obj')}
    save(evidence/'local_asset_index.json',dict(base='outputs/v235_validation',files=local,largeAssetsBackedUpToGitHub=False,
      requiredRuntime='v220 build02 thin archives and v167 exact objects; explicit future budget required'))
    for n in OWN:
        target=evidence/'final_source'/(n if n!='.gitattributes' else '.gitattributes.snapshot');target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,target)
    elapsed=time.monotonic()-start;assert totals['offline']+elapsed<120 and totals['build']<120 and totals['query']<60
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<128*1024**2
    report=dict(passed=True,goalStatus='active',physicalQualification=False,newQualifiedAreaM2=0,newHeatAreaM2=0,
      realRocketFlowUpdates=0,smallControlFlowUpdates=12,successfulControlUpdatesEach=5,
      observedAndDisabledFramesHistoryRestartBitEqual=True,nativeStages=60,sameStageWallRecords=10,
      readRejections=8,newRejectionCallbackDynamicallyExercised=False,
      firstBoundaryInconsistentPairAndSourceMappingFailuresPreserved=True,
      boundNativeSources=11,compilerDependencies=724,observationObjectsBuilt=2,fullEngineBuilt=False,
      originalWallFaces=19752,originalAreaM2=1.9646292565715553,assets=assets,
      protectedCount=len(protected),engineeringArtifactsUnchanged=len(eng['files']),candidateCount=len(candidates['files']),
      processSeconds=totals,finalizationSeconds=elapsed,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      outputBytes=size,allRecordedInputsStable=True,allRecordedGroupsEmpty=True,historicalSourceBindings=historical,
      closedCLIs=closed,fullRocketAdmitted=False,MATLABExecuted=False,sharedWebChanged=False,browserInteractionTested=False,
      uniqueAssetsDeleted=False,newRocketMeshGenerated=False,unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v235/evidence_manifest.json']+['docs/validation/v235/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,outputBytes=size)))


if __name__=='__main__':main()
