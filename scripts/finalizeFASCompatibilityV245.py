"""模块职责：核对保护资产、真实终态和准确执行源，分离并发本地补充并归档本轮交付。"""
from pathlib import Path
import ast,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v245_validation';DOC=ROOT/'docs/validation/v245'
CLIS=['scripts/'+n+'.py' for n in ('runFASCompatibilityV245','auditFASCompatibilityV245','renderFASCompatibilityV245')]
OWN=CLIS+['scripts/finalizeFASCompatibilityV245.py','native/su2PhaseAwareAir/FASHooksV245.hpp',
 'native/su2PhaseAwareAir/fasCompatibilityV245.cpp','docs/FAS_COMPATIBILITY_V245.md',
 'docs/validation/v245/PLAN.md','docs/validation/v245/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
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
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    prev=read(ROOT/'docs/validation/v244/evidence/final_report.json');assets=prev['assets'];unrelated=prev['unrelatedTrackedSHA256']
    for n,h in dict(assets,**unrelated).items():assert sha(ROOT/n)==h,n
    private=read(BASE/'concurrent_prefixes.json');split={}
    for n in SPLIT:
        t=(ROOT/n).read_text(encoding='utf8');p=private[n]['prefix'];assert t.count(p)==1 and t.replace(p,'',1)==source(n).read_text(encoding='utf8')
        assert 'public_benchmarks_20261116' not in source(n).read_text(encoding='utf8')
        split[n]=dict(localSupplementPreserved=True,publishedOwnVersionSHA256=sha(source(n)),combinedWorktreeSHA256=sha(ROOT/n))
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==7
    totals={};cache={};bindings=[];required=set()
    for name,r in records.items():
        assert r['passed'] and r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'] and r['exitCode']==0,name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for n,h in r['inputSHA256'].items():
            if n not in cache:cache[n]=sha(n)
            if cache[n]==h:continue
            snapshots=list(BASE.glob('process_*/execution_sources/'+h+Path(n).suffix));assert snapshots,(name,n)
            snapshot=snapshots[0];assert sha(snapshot)==h;required.add(snapshot)
            bindings.append(dict(process=name,path=n,executedSHA256=h,currentSHA256=cache[n],snapshot=snapshot.relative_to(BASE).as_posix()))
    assert totals['offline']<120 and totals['build']<180 and totals['query']<60 and totals['audit']<240
    audit=read(BASE/'audit01/report.json');assert audit['passed'] and not audit['integratedThermalCompatibilityPassed']
    assert audit['independentStateDomainChecks']==3548 and audit['independentComparedScalars']==3272
    assert audit['observationNeutralBitEqual'] and audit['wallCancellation'][0]['numericallyUnresolved']
    closed={}
    for name in OWN:
        if name.endswith('.py'):ast.parse(source(name).read_text(encoding='utf8'))
    for n in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        message=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in message;closed[n]=message.strip()
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        dest=evidence/(rel if rel is not None else p.relative_to(BASE));dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=p in required or p.suffix in ('.json','.txt','.csv','.png','.patch') or top.startswith('process_')
        if top=='stage01' and p.suffix in ('.cpp','.hpp'):keep=True
        if p.suffix in ('.bin','.exe','.obj','.a','.d','.su2','.pyc') or top=='owned_source':keep=False
        if p.name in ('concurrent_prefixes.json','update_docs.py'):keep=False
        if keep:assert p.stat().st_size<2*1024**2,str(p);copy(p)
        else:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for n in OWN:copy(source(n),Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,originalFaces=19752,originalAreaM2=1.9646292565715553,
      softwareAlgebraPassed=True,integratedThermalCompatibilityPassed=False,productionAdmitted=False,
      observedFrames=232,independentStates=3548,independentComparedScalars=3272,observedFlowSAUpdates=5,neutralControlFlowSAUpdates=5,
      observationNeutralBitEqual=True,coarseWallNormalRatio=audit['wallCancellation'][0]['ratio'],
      maxScaledComparisonError=max(audit['maxScaledErrors'].values()),scientificRejection='closed-shell wall agglomerated to a directionless coarse node',
      originalRocketFlowUpdates=0,newQualifiedAreaM2=0,physicalQualification=False,realInitializationPerformed=False,
      numericalObjectsChanged=0,readOnlyObservationObjectsCompiled=1,sharedWebGLChanged=False,MATLABExecuted=False,
      browserInteractionTested=False,scientificImageVisuallyInspected=True,firstPlotPreserved=True,
      protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      assets=assets,unrelatedTrackedSHA256=unrelated,splitDocumentationOwnership=split,concurrentLocalOnlyDataPublished=False,
      allReceiptsTerminalStableGroupsEmpty=True,executionFailures=0,scientificFailurePreserved=True,processSeconds=totals,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),historicalSourceBindings=bindings,
      closedCLIs=closed,uniqueAssetsDeleted=False,unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    report['outputBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputBytesBeforeFinalReport']<256*1024**2
    report['finalizationSeconds']=time.monotonic()-start;save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    digest={p.relative_to(DOC).as_posix():sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()};save(DOC/'evidence_manifest.json',digest)
    files=OWN+['docs/validation/v245/evidence_manifest.json']+[p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(dict(passed=True,files=len(files),evidenceFiles=len(digest),evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),processSeconds=totals))


if __name__=='__main__':main()
