"""模块职责：封存固定KINSOL控制的真实失败/回执/准确源，完整场与运行库只作本地索引。"""
import ast,json,subprocess,sys,time,shutil
from pathlib import Path
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v257_validation';DOC=ROOT/'docs/validation/v257';RT=ROOT/'tools/runtime/sundials-7.4.0-8e17876'
CLIS=['scripts/'+n+'.py' for n in ['runKinsolCoupledV257','runKinsolAuditV257','auditKinsolCoupledV257','exportKinsolCoupledV257']]
OWN=CLIS+['scripts/finalizeKinsolCoupledV257.py','native/su2PhaseAwareAir/FullResidualProviderV257.hpp',
    'native/su2PhaseAwareAir/kinsolCoupledV257.cpp','web/kinsol-coupled-v257.js','tests/testKinsolViewerV257.js',
    'docs/KINSOL_COUPLED_CONTROL_V257.md','docs/validation/v257/PLAN.md','docs/validation/v257/NEXT_DECISION.md',
    'docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md','.gitattributes']


def main():
    start=time.monotonic();ev=DOC/'evidence';out=BASE/'final01';assert not ev.exists() and not out.exists()
    baseline='83d382203653effaee0d1f117a21b69bb6de6653';assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==baseline
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    previous=read(ROOT/'docs/validation/v256/evidence/final_report.json')
    for name,h in previous['unrelatedTrackedSHA256'].items():assert sha(ROOT/name)==h,name
    receipts={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(receipts)==7
    snapshots={sha(p):p for p in BASE.glob('source_*/*') if p.is_file()}
    for p in BASE.glob('process_*/execution_sources/*'):snapshots.setdefault(sha(p),p)
    totals={};failures={};changes=[];cached={};ownExecuted=set();sourceIndex=[]
    for name,r in receipts.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert read(BASE/name/'terminal_capacity.json')['passed']
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        if not r['passed']:failures[name]=r['exitCode']
        for path,h in r['inputSHA256'].items():
            p=Path(path)
            if path not in cached:cached[path]=sha(p)
            if cached[path]!=h:
                assert h in snapshots,(path,h);changes.append(dict(process=name,path=path,executedSHA256=h,currentSHA256=cached[path]))
            if p in [ROOT/n for n in OWN]:
                assert h in snapshots,(path,h);ownExecuted.add(h)
            if h in snapshots:sourceIndex.append(dict(process=name,path=path,sha256=h,snapshot=snapshots[h].relative_to(BASE).as_posix()))
    assert set(failures)=={'process_prepare'}
    assert totals['offline']<=120 and totals['build']<=300 and totals['control']<=900 and totals['audit']<=120
    native=read(BASE/'rocket01/native_completion.json');audit=read(BASE/'audit01/report.json');ui=read(BASE/'workbench01/report.json')
    assert audit['passed'] and ui['passed'] and native['normalReturn'] and native['initialStateAndResidualBitwiseV256']
    assert native['nativeFlag']==-13 and native['nativeLastLinearFlag']==801 and native['linearGateRefused']
    assert native['observedCurrentStateChanges']==native['timeUpdates']==0 and native['nativeLinearIterations']==20
    assert audit['rawResidualRatios']==[1.]*6 and native['recoverableRejectedStates']==0
    assert sha(BASE/'rocket01/initial_state.bin')==sha(BASE/'rocket01/terminal_state.bin')
    assert sha(BASE/'rocket01/initial_full_residual.bin')==sha(BASE/'rocket01/terminal_full_residual.bin')
    closed={}
    for name in CLIS:
        result=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        text=result.stdout.decode('utf8',errors='replace');assert result.returncode!=0 and 'closed' in text;closed[name]=text.strip()
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    out.mkdir();ev.mkdir();archived={}
    def copy(p,rel):
        h=sha(p)
        if h in archived:return archived[h]
        q=ev/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);archived[h]=rel.as_posix();return rel.as_posix()
    for folder in ['stage01','stage02','rocket01','build01','audit01','workbench01']:
        for p in (BASE/folder).glob('*'):
            keep=p.suffix in ['.json','.csv','.cfg','.txt','.log','.png'] and p.name not in ['coupled.json','members.json'] and not p.name.endswith('_wall.csv')
            if p.is_file() and keep and p.stat().st_size<2*1024**2:copy(p,p.relative_to(BASE))
    for p in BASE.glob('process_*/*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    for h in ownExecuted:
        p=snapshots[h];copy(p,Path('execution_sources')/p.name)
    for item in sourceIndex:item['archivedSmallSource']=archived.get(item['sha256'])
    for item in changes:item['archivedSmallSource']=archived.get(item['executedSHA256']);assert item['archivedSmallSource']
    for p in (BASE/'research').rglob('*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    for name in ['LICENSE','NOTICE']:copy(RT/'source'/name,Path('external_license')/name)
    save(ev/'external_source.json',dict(repository='https://github.com/LLNL/sundials',tag='v7.4.0',commit='8e17876d3b4d682b4098684b07a85b005a122f81',
        license='BSD-3-Clause',retrievedVia='official fixed-commit codeload archive; exact GitHub connector comparisons',
        projectLocalInstalled=True,systemInstalled=False,linked=True,executed=True,libraryBuild=read(BASE/'stage02/library_build.json'),
        sourceManifest='stage02/library_source.json',libraryAndFullSourceArchivedInGit=False,
        previousKINSourceEvidence='docs/validation/v256/evidence/research/sundials',runtimePath=str(RT)))
    assets={}
    for p in BASE.rglob('*'):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        h=sha(p);assets[p.relative_to(BASE).as_posix()]=dict(bytes=p.stat().st_size,sha256=h,archivedSmallEvidencePath=archived.get(h))
    save(ev/'local_asset_index.json',assets);save(ev/'execution_source_index.json',sourceIndex)
    report=dict(passed=True,deliveryIntegrityPassed=True,numericalCombinationPassed=False,goalStatus='active',wholeSurfaceGoalComplete=False,baseCommit=baseline,
        originalFaces=19752,originalAreaM2=1.9646292565715553,originalNodes=335726,unknowns=1974844,
        actualObservedStateChanges=0,timeUpdates=0,newHeatSolutions=0,actualProfilesFilled=0,physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,
        KINSOLExecuted=True,KINSOLProjectLocalInstalled=True,nativeSummary=native,independentAudit=audit,workbenchSummary=ui,
        allReceiptsTerminalStableGroupsEmpty=True,executionFailures=failures,processSeconds=totals,
        maximumSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in receipts.values()),
        forecastMemoryBytes=read(BASE/'stage02/report.json')['forecastBytes'],historicalSourceBindings=changes,
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
        unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],scientificImagesVisuallyInspected=True,browserInteractionTested=False,
        localOnly=True,githubPushed=False,sharedWebGLChanged=False,uniqueDataDeleted=False,physicsObjectsRecompiled=0,
        fullFieldsArchivedInGit=False,nativeLibraryAndToolchainArchivedInGit=False,closedCLIs=closed,ownedFiles=OWN,
        # 当前入口和索引在证据落盘后更新，属于可变交接入口，不放进本快照SHA。
        currentOwnedFileSHA256={n:sha(ROOT/n) for n in OWN if n not in ['docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md']},
        softReserveTimerExcludesNativeConstruction=True,softReserveRiskDidNotTrigger=True,developmentReadsArchivingSeparate=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputLogicalBytesBeforeFinalReport']<1024**3
    report['finalizationSeconds']=time.monotonic()-start;save(out/'report.json',report);save(ev/'final_report.json',report)
    save(ev/'manifest.json',{p.relative_to(ev).as_posix():sha(p) for p in ev.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ['deliveryIntegrityPassed','numericalCombinationPassed','actualObservedStateChanges','processSeconds','finalizationSeconds']}))


if __name__=='__main__':main()
