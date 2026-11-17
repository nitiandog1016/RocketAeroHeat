"""模块职责：封存一次方向诊断的真实结果/准确源和失败结论，保留大场本地索引。"""
import ast,json,subprocess,sys,time,shutil
from pathlib import Path
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v258_validation';DOC=ROOT/'docs/validation/v258'
CLIS=['scripts/'+n+'.py' for n in ['runDirectionConsistencyV258','runDirectionAuditV258','auditDirectionConsistencyV258','exportDirectionConsistencyV258']]
OWN=CLIS+['scripts/finalizeDirectionConsistencyV258.py','native/su2PhaseAwareAir/directionConsistencyV258.cpp',
    'web/direction-consistency-v258.js','tests/testDirectionViewerV258.js','docs/DIRECTION_CONSISTENCY_V258.md',
    'docs/validation/v258/PLAN.md','docs/validation/v258/NEXT_DECISION.md','docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md','.gitattributes']


def main():
    start=time.monotonic();ev=DOC/'evidence';out=BASE/'final01';assert not ev.exists() and not out.exists()
    baseline='29d7f95a2a01ba61cdf3db46799e88fab96f94ca';assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==baseline
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    previous=read(ROOT/'docs/validation/v257/evidence/final_report.json')
    for name,h in previous['unrelatedTrackedSHA256'].items():assert sha(ROOT/name)==h,name
    receipts={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(receipts)==5
    snapshots={sha(p):p for p in BASE.glob('source_*/*') if p.is_file()}
    for p in BASE.glob('process_*/execution_sources/*'):snapshots.setdefault(sha(p),p)
    totals={};changes=[];cached={};ownExecuted=set();sourceIndex=[]
    for name,r in receipts.items():
        assert r['passed'] and r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert read(BASE/name/'terminal_capacity.json')['passed'];totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for path,h in r['inputSHA256'].items():
            p=Path(path)
            if path not in cached:cached[path]=sha(p)
            if cached[path]!=h:
                assert h in snapshots,(path,h);changes.append(dict(process=name,path=path,executedSHA256=h,currentSHA256=cached[path]))
            if p in [ROOT/n for n in OWN]:assert h in snapshots,(path,h);ownExecuted.add(h)
            if h in snapshots:sourceIndex.append(dict(process=name,path=path,sha256=h,snapshot=snapshots[h].relative_to(BASE).as_posix()))
    assert totals['offline']<=60 and totals['build']<=120 and totals['control']<=180 and totals['audit']<=60
    native=read(BASE/'rocket01/native_completion.json');audit=read(BASE/'audit01/report.json');ui=read(BASE/'workbench01/report.json')
    assert native['normalReturn'] and native['fullPhysicalEvaluations']==4 and native['nonlinearUpdates']==0
    assert audit['passed'] and not audit['fixedDirectionHalfStepConsistent'] and audit['strongRowsExactlyZero']
    assert not audit['sourceSpecificInterfaceBugProved'] and ui['passed'] and ui['differenceIsNotExactJacobianError']
    closed={}
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        text=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in text;closed[name]=text.strip()
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    out.mkdir();ev.mkdir();archived={}
    def copy(p,rel):
        h=sha(p)
        if h in archived:return archived[h]
        q=ev/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);archived[h]=rel.as_posix();return rel.as_posix()
    for folder in ['stage01','rocket01','build01','audit01','workbench01']:
        for p in (BASE/folder).glob('*'):
            keep=p.suffix in ['.json','.csv','.cfg','.txt','.png'] and p.name!='direction.json'
            if p.is_file() and keep and p.stat().st_size<2*1024**2:copy(p,p.relative_to(BASE))
    for p in BASE.glob('process_*/*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    for h in ownExecuted:
        p=snapshots[h];copy(p,Path('execution_sources')/p.name)
    for item in sourceIndex:item['archivedSmallSource']=archived.get(item['sha256'])
    for item in changes:item['archivedSmallSource']=archived.get(item['executedSHA256']);assert item['archivedSmallSource']
    assets={}
    for p in BASE.rglob('*'):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        h=sha(p);assets[p.relative_to(BASE).as_posix()]=dict(bytes=p.stat().st_size,sha256=h,archivedSmallEvidencePath=archived.get(h))
    save(ev/'local_asset_index.json',assets);save(ev/'execution_source_index.json',sourceIndex)
    report=dict(passed=True,deliveryIntegrityPassed=True,fixedDirectionConsistencyPassed=False,goalStatus='active',wholeSurfaceGoalComplete=False,baseCommit=baseline,
        originalFaces=19752,originalAreaM2=1.9646292565715553,originalNodes=335726,unknowns=1974844,nonlinearUpdates=0,timeUpdates=0,
        newHeatSolutions=0,physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,KINSOLExecuted=False,
        nativeSummary=native,independentAudit=audit,workbenchSummary=ui,allReceiptsTerminalStableGroupsEmpty=True,
        executionFailures={},processSeconds=totals,maximumSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in receipts.values()),
        forecastMemoryBytes=read(BASE/'stage01/report.json')['forecastBytes'],historicalSourceBindings=changes,
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
        unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],scientificImagesVisuallyInspected=True,browserInteractionTested=False,
        localOnly=True,githubPushed=False,sharedWebGLChanged=False,uniqueDataDeleted=False,physicsObjectsRecompiled=0,matureLibrariesRecompiled=0,
        fullFieldsArchivedInGit=False,nativeLibraryAndToolchainArchivedInGit=False,closedCLIs=closed,ownedFiles=OWN,
        currentOwnedFileSHA256={n:sha(ROOT/n) for n in OWN if n not in ['docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md']},
        softwarePreparationAndSourceChecksNotPhysicalQualification=True,softBudgetIncludesNativeConstruction=True,developmentReadsArchivingSeparate=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputLogicalBytesBeforeFinalReport']<512*1024**2
    report['finalizationSeconds']=time.monotonic()-start;save(out/'report.json',report);save(ev/'final_report.json',report)
    save(ev/'manifest.json',{p.relative_to(ev).as_posix():sha(p) for p in ev.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ['deliveryIntegrityPassed','fixedDirectionConsistencyPassed','nonlinearUpdates','processSeconds','finalizationSeconds']}))


if __name__=='__main__':main()
