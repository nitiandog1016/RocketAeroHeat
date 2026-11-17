"""模块职责：按SHA去重封存真实证据，保留大连续解索引，保护并发内容，仅本地交付。"""
import ast,json,subprocess,sys,shutil,time
from pathlib import Path
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v255_validation';DOC=ROOT/'docs/validation/v255'
CLIS=['scripts/'+s+'.py' for s in ['prepareRegularizedEntropyV255','runRegularizedEntropyV255','auditRegularizedEntropyV255','exportRegularizedEntropyV255']]
OWN=CLIS+['scripts/finalizeRegularizedEntropyV255.py','solveRegularizedEntropyIntegral.m','phaseAwareAirCaloricDrop.m',
    'tests/testRegularizedEntropyV255.m','tests/testRegularizedEntropyViewerV255.js','web/regularized-entropy-v255.js','web/regularized-entropy-v255.html',
    'docs/REGULARIZED_ENTROPY_INTEGRAL_V255.md','docs/validation/v255/PLAN.md','docs/validation/v255/NEXT_DECISION.md',
    'docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md','.gitattributes']


def main():
    start=time.monotonic();out=BASE/'final01';ev=DOC/'evidence';assert not out.exists() and not ev.exists()
    baseline='fe1434e30005c33a4fb5cf12c3fcf48e40cf5756';assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==baseline
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for p,h in protected.items():assert sha(ROOT/p)==h,p
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    previous=read(ROOT/'docs/validation/v254/evidence/final_report.json')
    for p,h in previous['unrelatedTrackedSHA256'].items():assert sha(ROOT/p)==h,p
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==4
    totals={};sources={};cache={};changed=[]
    for p in BASE.glob('source_*/*'):
        if p.is_file():sources[sha(p)]=p
    for name,r in records.items():
        assert r['passed'] and r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert read(BASE/name/'terminal_capacity.json')['passed']
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if p not in cache:cache[p]=sha(p)
            if cache[p]!=h:
                assert h in sources,(p,h)
                changed.append(dict(process=name,path=p,executedSHA256=h,currentSHA256=cache[p]))
    assert totals['offline']<=60 and totals['matlab']<=180 and totals['audit']<=120
    matlab=read(BASE/'matlab01/report.json');audit=read(BASE/'audit01/report.json');ui=read(BASE/'workbench01/report.json')
    assert all(r['passed'] for r in [matlab,audit,ui])
    assert ui['candidateFaces']==audit['admittedOriginalFaces']==1306 and len(audit['recoveredPreviouslyWithheldFaceIDs'])==31
    assert not ui['physicalQualification'] and not ui['profileReadyForCFD']
    closed={}
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        text=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in text;closed[name]=text.strip()
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    out.mkdir();ev.mkdir();stored={};assets={}
    def copy(p,rel):
        h=sha(p)
        if h in stored:return stored[h]
        target=ev/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
        stored[h]=rel.as_posix();return rel.as_posix()
    # 计算输出先入库；同SHA执行输入可指向该唯一副本，不再次复制。
    for folder in ['prepare01','matlab01','audit01','workbench01']:
        for p in (BASE/folder).glob('*'):
            keep=p.suffix=='.json' and p.name!='coupled.json'
            keep=keep or p.suffix=='.png' or p.name=='frontend_test.txt'
            if p.is_file() and keep:copy(p,p.relative_to(BASE))
    for p in BASE.glob('process_*/*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    sourceIndex=[]
    for h,p in sorted(sources.items()):
        archived=copy(p,Path('execution_sources')/(h+p.suffix));sourceIndex.append(dict(sha256=h,originalSnapshot=p.relative_to(BASE).as_posix(),archivedPath=archived))
    for r in changed:r['archivedPath']=stored[r['executedSHA256']]
    for p in BASE.rglob('*'):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        h=sha(p);assets[p.relative_to(BASE).as_posix()]=dict(bytes=p.stat().st_size,sha256=h,archivedSmallEvidencePath=stored.get(h))
    save(ev/'local_asset_index.json',assets);save(ev/'execution_source_index.json',sourceIndex)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,baseCommit=baseline,
        originalFaces=19752,originalAreaM2=1.9646292565715553,candidateFaces=1306,recoveredPriorRejectedFaces=31,
        actualProfilesFilled=0,profileReadyForCFD=False,physicalQualification=False,newQualifiedAreaM2=0,
        fullSampleNumericalGatesPassed=True,positiveBelowSeedNumericallySolved=False,matlabSummary=matlab,independentAudit=audit,workbenchSummary=ui,
        allReceiptsTerminalStableGroupsEmpty=True,executionFailures={},processSeconds=totals,maximumSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
        historicalSourceBindings=changed,allRequiredChangedSourcesSnapshottedBeforeExecution=True,
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],
        scientificImagesVisuallyInspected=True,browserInteractionTested=False,sharedWebGLChanged=False,
        localOnly=True,githubPushed=False,newCFD=False,newPackageInstalled=False,uniqueDataDeleted=False,
        nativeContinuousMATArchivedInGit=False,oldV252ParentTerminalKnown=False,closedCLIs=closed,ownedFiles=OWN,
        currentOwnedFileSHA256={p:sha(ROOT/p) for p in OWN},developmentReadsArchivingSeparate=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputLogicalBytesBeforeFinalReport']<128*1024**2
    report['finalizationSeconds']=time.monotonic()-start
    save(out/'report.json',report);save(ev/'final_report.json',report)
    save(ev/'manifest.json',{p.relative_to(ev).as_posix():sha(p) for p in ev.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ['passed','candidateFaces','recoveredPriorRejectedFaces','actualProfilesFilled','processSeconds','finalizationSeconds']}))


if __name__=='__main__':main()
