"""模块职责：封存真实失败、独立审核和科研出口，核对保护资产并生成精确公开清单。"""
from pathlib import Path
import ast,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v247_validation';DOC=ROOT/'docs/validation/v247'
CLIS=['scripts/'+n+'.py' for n in ('runFASDevelopmentV247','auditFASDevelopmentV247','exportFASDevelopmentV247','reviewFASRecoveryV247')]
OWN=CLIS+['scripts/fasDevelopmentFieldsV247.py','scripts/finalizeFASDevelopmentV247.py',
 'native/su2PhaseAwareAir/fasDevelopmentV247.cpp','tests/testFASDevelopmentFramesV247.py','tests/testFASDevelopmentV247.js',
 'docs/FAS_BOUNDED_DEVELOPMENT_V247.md','docs/validation/v247/PLAN.md','docs/validation/v247/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
SPLIT=('AGENTS.md','docs/PROJECT_GOALS.md')


def source(n):return BASE/'owned_source'/n if n in SPLIT else ROOT/n


def main():
    begin=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    prev=read(ROOT/'docs/validation/v246/evidence/final_report.json');assets=prev['assets'];unrelated=prev['unrelatedTrackedSHA256']
    for n,h in dict(assets,**unrelated).items():assert sha(ROOT/n)==h,n
    private=read(BASE/'concurrent_prefixes.json');split={}
    for n in SPLIT:
        live=(ROOT/n).read_text(encoding='utf8');prefix=private[n]['prefix']
        assert live.count(prefix)==1 and live.replace(prefix,'',1)==source(n).read_text(encoding='utf8')
        assert 'public_benchmarks_20261116' not in source(n).read_text(encoding='utf8')
        split[n]=dict(localSupplementPreserved=True,publishedOwnVersionSHA256=sha(source(n)),combinedWorktreeSHA256=sha(ROOT/n))
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==11
    failures={'process_prepare':1,'process_rocket':86,'process_audit':1,'process_audit2':1,'process_method_review':1}
    totals={};cache={};bindings=[];required=set()
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert r['exitCode']==failures.get(name,0) and r['passed']==(name not in failures)
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for n,h in r['inputSHA256'].items():
            if n not in cache:cache[n]=sha(n)
            if cache[n]==h:continue
            matches=list(BASE.glob('process_*/execution_sources/'+h+Path(n).suffix))
            matches += [p for p in (BASE/'failed_source01').iterdir() if p.is_file() and sha(p)==h]
            assert matches,(name,n);snapshot=matches[0];required.add(snapshot)
            bindings.append(dict(process=name,path=n,executedSHA256=h,currentSHA256=cache[n],snapshot=snapshot.relative_to(BASE).as_posix()))
    assert totals['offline']<120 and totals['build']<120 and totals['cfd']<1800 and totals['audit']<240
    audit=read(BASE/'audit03/report.json');method=read(BASE/'method_review/report.json');ui=read(BASE/'workbench01/report.json')
    assert audit['passed'] and not audit['developmentSucceeded'] and audit['completedUpdates']==60 and audit['terminal'] is None
    assert audit['completeSavedUpdates']==[0,1] and audit['originalFaceUpdates']==[0,1,50]
    assert method['passed'] and method['actualNativeRecoverySucceededBeforeObserverAbort'] and ui['passed']
    assert not (BASE/'rocket01/native_completion.json').exists() and not (BASE/'rocket01/restart_flow.csv').exists()
    closed={}
    for n in OWN:
        if n.endswith('.py'):ast.parse(source(n).read_text(encoding='utf8'))
    for n in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        message=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in message;closed[n]=message.strip()
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        dest=evidence/(rel if rel is not None else p.relative_to(BASE));dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=p in required or p.suffix in ('.json','.txt','.csv','.png') or top.startswith('process_') or top=='failed_source01'
        if p.suffix in ('.bin','.exe','.d','.su2','.pyc','.c','.extra_newline') or top=='owned_source':keep=False
        if p.name in ('concurrent_prefixes.json','update_docs.py','write_native.py','write_export.py','rocket_fas_control.json') or p.name.startswith('original_faces_'):keep=False
        if p.stat().st_size>=2*1024**2:keep=False
        if keep:copy(p)
        else:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for n in OWN:copy(source(n),Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,developmentSucceeded=False,
      originalFaces=19752,originalAreaM2=1.9646292565715553,finePoints=335726,coarsePoints=44484,completedUpdates=60,
      completeSavedUpdates=[0,1],originalFaceUpdates=[0,1,50],failureUpdate=61,nativeRecoveredLegalOldFiveState=True,
      observerAbortedByContract=True,intrinsicUnrecoverableSolverCrashProved=False,actualPreProlongationStateMissing=True,
      terminalResidualCompiledNotExecuted=True,normalFinalCSVNotReached=True,benefitDecisionNotReached=True,
      initialZeroBitEqual=True,firstUpdateBitEqual=False,firstUpdateMaximumScaledDifference=audit['errorMetrics']['crossThinExecutableFirstUpdate'],
      independentSavedStates=audit['independentSavedStateChecks'],brentCrosschecks=audit['brentCrosschecks'],
      scientificImageVisuallyInspected=True,sourceAndNumericTestsPassed=True,physicalQualification=False,newQualifiedAreaM2=0,
      exportedFields=ui['fields'],originalFaceValuesVerified=ui['exactResearchValuesChecked'],
      numericalObjectsRecompiled=0,sharedWebGLChanged=False,MATLABExecuted=False,browserInteractionTested=False,
      protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      assets=assets,unrelatedTrackedSHA256=unrelated,splitDocumentationOwnership=split,concurrentLocalOnlyDataPublished=False,
      allReceiptsTerminalStableGroupsEmpty=True,executionFailuresPreserved=failures,documentParseFailureBeforeExecutionPreserved=True,
      processSeconds=totals,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      historicalSourceBindings=bindings,closedCLIs=closed,uniqueDataDeleted=False,newPackageInstalled=False,
      nativeRecoveryNextAdopted=True,PETScSourceOnlyAlternative=True,unmeteredDevelopmentSourceReadsNetworkAndPublicationDisclosed=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputLogicalBytesBeforeFinalReport']<1536*1024**2
    report['finalizationSeconds']=time.monotonic()-begin;save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    files=OWN+['docs/validation/v247/evidence_manifest.json']+[p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(dict(passed=True,files=len(files),evidenceFiles=len(manifest),evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),processSeconds=totals))


if __name__=='__main__':main()
