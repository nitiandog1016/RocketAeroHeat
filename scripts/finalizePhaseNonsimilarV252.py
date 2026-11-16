"""模块职责：完整绑定实际源/恢复缺口/修复回归，发布小证据并保留并发本地内容。"""
from pathlib import Path
import ast,json,subprocess,sys,shutil,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v252_validation';DOC=ROOT/'docs/validation/v252'
CLIS=['scripts/'+n+'.py' for n in ('runPhaseNonsimilarV252','runPhaseNonsimilarAuditV252','auditPhaseNonsimilarV252','exportPhaseNonsimilarV252')]
OWN=CLIS+['scripts/finalizePhaseNonsimilarV252.py','solvePhaseAwareNonsimilarStation.m',
 'marchPhaseAwareLaminarProfiles.m','phaseAwareNoseMarchInputs.m','mapPhaseAwareMarchToSurface.m',
 'tests/testPhaseNonsimilarV252.m','tests/testPhaseNonsimilarGuardsV252.m','tests/testPhaseNonsimilarViewerV252.js',
 'web/phase-profiles-v252.js','web/phase-profiles-v252.html','docs/PHASE_NONSIMILAR_PROFILES_V252.md',
 'docs/validation/v252/PLAN.md','docs/validation/v252/NEXT_DECISION.md','docs/validation/v252/update_goals.py',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
SPLIT=('AGENTS.md','docs/PROJECT_GOALS.md')


def source(n):return BASE/'owned_source'/n if n in SPLIT else ROOT/n


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    previous=read(ROOT/'docs/validation/v251/evidence/final_report.json')
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for group in (previous['assets'],previous['unrelatedTrackedSHA256'],protected):
        for n,h in group.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for item in engineering['files']:assert sha(ROOT/item['path'])==item['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for item in candidates['files']:
        assert sha(ROOT/item['file'])==item['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/item['file'])==item['candidateSHA256']
    split={};private=read(BASE/'concurrent_prefixes.json')
    for n in SPLIT:
        current=(ROOT/n).read_text(encoding='utf8');prefix=private[n]['prefix']
        assert current.count(prefix)==1 and current.replace(prefix,'',1)==source(n).read_text(encoding='utf8')
        assert 'public_benchmarks_20261116' not in source(n).read_text(encoding='utf8')
        split[n]=dict(localSupplementPreserved=True,publishedSHA256=sha(source(n)),worktreeSHA256=sha(ROOT/n))
    # 第一出口JS随后修复了切换分支的位置选择；按准确旧片段重建，必须与准入SHA一致。
    js=ROOT/'web/phase-profiles-v252.js';current=js.read_text(encoding='utf8')
    begin=current.index('  let lastBranch=branch.value;');end=current.index('\n  document.getElementById',begin)
    old=current[:begin]+"  branch.onchange=()=>{slider.value=Math.min(Number(slider.value),data[branch.value].stations.length-1);render();};slider.oninput=render;window.addEventListener('resize',render);"+current[end:]
    oldpath=BASE/'prior_source/phase-profiles-v252.js';oldpath.parent.mkdir(exist_ok=True);oldpath.write_bytes(old.encode('utf8'))
    assert sha(oldpath)==read(BASE/'process_export/admission.json')['inputSHA256'][str(js)]
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')}
    assert len(records)==6 and not (BASE/'process_matlab01/report.json').exists()
    totals={};bindings=[];cache={}
    snapshotpaths=list(BASE.glob('process_*/execution_sources/*'))+list(BASE.glob('matlab01/executed_source/**/*.m'))+list(BASE.glob('guards01/executed_source/**/*.m'))+[oldpath]
    snapshotmap={}
    for p in snapshotpaths:
        if p.is_file():snapshotmap.setdefault(sha(p),p)
    # 缺原父终态不绕过输入绑定，仍核对该任务准入及保存原文。
    for folder in sorted(BASE.glob('process_*')):
        if not folder.is_dir():continue
        admission=read(folder/'admission.json')
        for n,h in admission['inputSHA256'].items():
            if n not in cache:cache[n]=sha(n)
            if cache[n]!=h:
                assert h in snapshotmap,(folder.name,n,h)
                bindings.append(dict(process=folder.name,path=n,executedSHA256=h,currentSHA256=cache[n],snapshot=snapshotmap[h].relative_to(BASE).as_posix()))
        if folder.name in records:
            r=records[folder.name];assert r['passed'] and r['inputsStable'] and not r['liveGroupAfterCleanup']
            assert read(folder/'terminal_capacity.json')['passed']
            totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
    assert totals['offline']<=120 and totals['audit']<=180 and totals['matlab']<=60
    audit=read(BASE/'audit01/report.json');ui=read(BASE/'workbench02/report.json');guard=read(BASE/'guards01/report.json');recovery=read(BASE/'recovery01/report.json')
    assert audit['allNumericalGatesPassed'] and ui['passed'] and guard['passed'] and recovery['passed']
    assert guard['maximumScaledDifference']==0 and guard['mappingRejections']==5
    assert all(sha(BASE/'matlab01'/n)==h for n,h in recovery['outputSHA256'].items())
    closed={}
    for n in CLIS:
        proc=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15)
        msg=proc.stdout.decode('utf8',errors='replace');assert proc.returncode!=0 and 'closed' in msg;closed[n]=msg.strip()
    for n in OWN:
        if n.endswith('.py'):ast.parse(source(n).read_text(encoding='utf8'))
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel):
        dest=evidence/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=top.startswith('process_') or top in ('method01','recovery01','audit01','prior_source','finalization_attempt01')
        keep=keep or (top in ('matlab01','guards01') and ('executed_source' in rel.parts or p.name=='report.json'))
        keep=keep or (top.startswith('workbench') and (p.suffix=='.png' or p.name in ('report.json','frontend_test.txt')))
        if keep:assert p.stat().st_size<2*1024**2;copy(p,rel)
        elif top!='owned_source' and p.name!='concurrent_prefixes.json':local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for n in OWN:copy(source(n),Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,originalFaces=19752,
        originalAreaM2=1.9646292565715553,candidateFaces=1306,candidateAreaM2=.003518440977521075,
        newQualifiedAreaM2=0,physicalQualification=False,independentNumericalAudit=audit,
        originalMatlabSupervisorTerminalKnown=False,originalMatlabExitCode=None,originalMatlabPeakBytes=None,
        recoveryObservationOnly=True,originalMatlabBudgetDebitedSeconds=300,
        correctedSourceSingleStationRegression=guard,processSeconds=totals,
        conservativelyDebitedMatlabSeconds=300+totals['matlab'],
        existingReceiptsStableGroupsEmpty=True,missingOriginalReceiptNotFabricated=True,
        originalResultsUnchanged=True,sharedWebGLChanged=False,browserInteractionTested=False,
        scientificImagesVisuallyInspected=True,firstPlotLayoutFailurePreserved=True,
        fields=ui['fields'],originalFaceValuesVerified=ui['exactSurfaceValuesChecked'],
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
        assets=previous['assets'],unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],
        splitDocumentationOwnership=split,concurrentLocalOnlyDataPublished=False,
        historicalSourceBindings=bindings,firstProfileJSSourceRestoredWithExactAdmissionSHA=True,
        closedCLIs=closed,ownedFiles=OWN,uniqueDataDeleted=False,newPackageInstalled=False,newCFD=False,
        maxKnownSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
        unmeteredDevelopmentReadsAndPublicationDisclosed=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputLogicalBytesBeforeFinalReport']<256*1024**2
    report['finalizationSeconds']=time.monotonic()-start
    save(out/'report.json',report);save(evidence/'final_report.json',report)
    save(evidence/'manifest.json',{p.relative_to(evidence).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ('passed','fields','newQualifiedAreaM2','processSeconds','outputLogicalBytesBeforeFinalReport','finalizationSeconds')}))


if __name__=='__main__':main()
