"""模块职责：绑定真实FAS控制、失败和原面出口，保护并发工作并归档小证据。"""
from pathlib import Path
import ast,os,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v246_validation';DOC=ROOT/'docs/validation/v246'
CLIS=['scripts/'+n+'.py' for n in ('runFASRocketPreflightV246','auditFASRocketPreflightV246','runFASRocketControlV246',
 'auditFASRocketControlV246','exportFASRocketV246','reviewFASRealTransfersV246')]
OWN=CLIS+['scripts/finalizeRealFASV246.py','native/su2PhaseAwareAir/fasRocketPreflightV246.cpp','native/su2PhaseAwareAir/fasRocketControlV246.cpp',
 'tests/testFASRocketV246.js','docs/REAL_FAS_CONTROL_V246.md','docs/validation/v246/PLAN.md','docs/validation/v246/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
SPLIT=('AGENTS.md','docs/PROJECT_GOALS.md')


def source(n):return BASE/'owned_source'/n if n in SPLIT else ROOT/n


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    prev=read(ROOT/'docs/validation/v245/evidence/final_report.json');assets=prev['assets'];unrelated=prev['unrelatedTrackedSHA256']
    for n,h in dict(assets,**unrelated).items():assert sha(ROOT/n)==h,n
    private=read(BASE/'concurrent_prefixes.json');split={}
    for n in SPLIT:
        t=(ROOT/n).read_text(encoding='utf8');prefix=private[n]['prefix'];assert t.count(prefix)==1 and t.replace(prefix,'',1)==source(n).read_text(encoding='utf8')
        assert 'public_benchmarks_20261116' not in source(n).read_text(encoding='utf8')
        split[n]=dict(localSupplementPreserved=True,publishedOwnVersionSHA256=sha(source(n)),combinedWorktreeSHA256=sha(ROOT/n))
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==14
    failures={'process_build':1,'process_control_audit':1};totals={};bindings=[];cache={};required=set()
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert r['exitCode']==failures.get(name,0) and r['passed']==(name not in failures)
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for n,h in r['inputSHA256'].items():
            if n not in cache:cache[n]=sha(n)
            if cache[n]==h:continue
            matches=list(BASE.glob('process_*/execution_sources/'+h+Path(n).suffix))
            matches += [p for p in (BASE/'failed_source01').iterdir() if sha(p)==h]
            assert matches,(name,n);snapshot=matches[0];assert sha(snapshot)==h;required.add(snapshot)
            bindings.append(dict(process=name,path=n,executedSHA256=h,currentSHA256=cache[n],snapshot=snapshot.relative_to(BASE).as_posix()))
    assert totals['offline']<120 and totals['audit']<240 and totals['query']<120 and totals['cfd']<180 and totals['build']<120
    for n in ('audit01','control_audit02','transfer_audit01','workbench01'):assert read(BASE/n/'report.json')['passed'],n
    control=read(BASE/'control_audit02/report.json');assert control['completeNativeFiveStepCompatibilityPassed'] and not control['physicalQualification']
    assert control['savedFineFlowArrayIsCorrection'] and control['finalCycleResidualNotSaved'] and not control['physicalMinus8All']
    geometry=read(BASE/'audit01/report.json');assert geometry['geometryPreflightPassed'] and geometry['coarsePoints']==44484
    dedup=read(BASE/'deduplication_corrected.json');assert dedup['duplicatedContentBytes']==535916990 and not dedup['uniqueDataDeleted']
    for r in dedup['files']:
        assert os.path.samefile(r['source'],r['target']) and sha(r['target'])==r['sha256']
    closed={}
    for n in OWN:
        if n.endswith('.py'):ast.parse(source(n).read_text(encoding='utf8'))
    for n in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        msg=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in msg;closed[n]=msg.strip()
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        dest=evidence/(rel if rel is not None else p.relative_to(BASE));dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=p in required or p.suffix in ('.json','.txt','.csv','.png') or top.startswith('process_') or top=='failed_source01'
        if p.name=='deduplicate_verified.ps1':keep=True
        if p.suffix in ('.bin','.exe','.d','.su2','.pyc') or top=='owned_source':keep=False
        if p.name.startswith('original_faces_') or p.name in ('concurrent_prefixes.json','update_docs.py','rocket_fas_control.json'):keep=False
        if p.stat().st_size>=2*1024**2:keep=False
        if keep:copy(p)
        else:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for n in OWN:copy(source(n),Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,originalFaces=19752,originalAreaM2=1.9646292565715553,
      realGeometryPreflightPassed=True,finePoints=335726,coarsePoints=44484,minimumCoarseWallNormalRatio=geometry['boundary'][0]['minimumNormalRatio'],
      actualFiveStepFlowSACompatibilityPassed=True,observationNeutralBitEqual=True,originalRocketFlowUpdatesPerRun=5,controlRuns=2,
      nativeJointDomainStateVisits=19010500,independentControlSavedStateChecks=1096146,independentPreflightSavedStateChecks=760420,
      trueSARestrictionAndWallProjectionPassed=True,realTransferPhysicalEnergyPassed=True,
      savedFineArrayIsCorrection=True,completeFinalCycleResidualMissing=True,physicalQualification=False,newQualifiedAreaM2=0,
      actualWholeSurfaceResearchStages=[0,1,5],exportedFields=14,originalFaceValuesVerified=237024,
      sourceAndNumericTestsPassed=True,numericalObjectsRecompiled=0,sharedWebGLChanged=False,MATLABExecuted=False,
      browserInteractionTested=False,scientificImageVisuallyInspected=True,
      protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      assets=assets,unrelatedTrackedSHA256=unrelated,splitDocumentationOwnership=split,concurrentLocalOnlyDataPublished=False,
      allReceiptsTerminalStableGroupsEmpty=True,executionFailuresPreserved=failures,exportParseFailureBeforeSupervisionPreserved=True,
      duplicateContentLinkedBytes=dedup['duplicatedContentBytes'],uniqueDataDeleted=False,originalNullDedupTotalPreserved=True,
      actualAllocatedDiskSavingsMeasured=False,processSeconds=totals,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      historicalSourceBindings=bindings,closedCLIs=closed,unmeteredDevelopmentTargetedReadsAndDedupMaintenanceDisclosed=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputLogicalBytesBeforeFinalReport']<1536*1024**2
    report['finalizationSeconds']=time.monotonic()-start;save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    digest={p.relative_to(DOC).as_posix():sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()};save(DOC/'evidence_manifest.json',digest)
    files=OWN+['docs/validation/v246/evidence_manifest.json']+[p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(dict(passed=True,files=len(files),evidenceFiles=len(digest),evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),processSeconds=totals))


if __name__=='__main__':main()
