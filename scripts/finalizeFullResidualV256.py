"""模块职责：封存完整域接口的成功与失败，准确源按范围去重，大场/内核源留本地。"""
import ast,json,subprocess,sys,time,shutil
from pathlib import Path
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v256_validation';DOC=ROOT/'docs/validation/v256'
CLIS=['scripts/'+n+'.py' for n in ['runFullResidualV256','runFullResidualAuditV256','auditFullResidualV256','exportFullResidualV256']]
OWN=CLIS+['scripts/finalizeFullResidualV256.py','native/su2PhaseAwareAir/fullDomainResidualV256.cpp',
    'web/full-residual-v256.js','tests/testFullResidualViewerV256.js','docs/FULL_DOMAIN_RESIDUAL_V256.md',
    'docs/validation/v256/PLAN.md','docs/validation/v256/NEXT_DECISION.md','docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md','.gitattributes']


def main():
    started=time.monotonic();out=BASE/'final01';ev=DOC/'evidence';assert not out.exists() and not ev.exists()
    baseline='34643c2d0fa9e950cc94587bd7f2e5ee2fb642fe';assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==baseline
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    previous=read(ROOT/'docs/validation/v255/evidence/final_report.json')
    for name,h in previous['unrelatedTrackedSHA256'].items():assert sha(ROOT/name)==h,name
    receipts={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(receipts)==7
    snapshots={sha(p):p for p in BASE.glob('source_*/*') if p.is_file()};totals={};failures={};changes=[];cached={}
    ownExecuted=set();sourceIndex=[]
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
    assert set(failures)=={'process_prepare','process_control'}
    assert totals['offline']<=120 and totals['build']<=120 and totals['control']<=180 and totals['audit']<=120
    native=read(BASE/'rocket02/native_completion.json');audit=read(BASE/'audit01/report.json');ui=read(BASE/'workbench01/report.json')
    assert audit['passed'] and ui['passed'] and native['normalReturn']
    assert native['repeatBitMismatches']==native['restoredBitMismatches']==native['actualTimeUpdates']==native['nonlinearIterations']==0
    assert audit['nativeReconstructionWarningCounts']==[3904] and not audit['allSixResidualsBelow1eMinus8']
    assert not audit['physicalQualification'] and not ui['newHeatSolutions']
    closed={}
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        s=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in s;closed[name]=s.strip()
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    out.mkdir();ev.mkdir();archived={}
    def copy(p,rel):
        h=sha(p)
        if h in archived:return archived[h]
        q=ev/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);archived[h]=rel.as_posix();return rel.as_posix()
    # 只打包本轮小证据、准确适配器/执行计划；整套工具链及成熟对象仍在E盘索引。
    for folder in ['stage01','stage02','rocket01','rocket02','build01','audit01','workbench01']:
        for p in (BASE/folder).glob('*'):
            keep=p.suffix in ['.json','.csv','.cfg','.txt','.png'] and p.name not in ['residual.json','dependencies.json','members.json']
            if p.is_file() and keep and p.stat().st_size<2*1024**2:copy(p,p.relative_to(BASE))
    for p in BASE.glob('process_*/*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    for h in ownExecuted:
        p=snapshots[h];copy(p,Path('execution_sources')/p.name)
    for item in sourceIndex:item['archivedSmallSource']=archived.get(item['sha256'])
    for item in changes:item['archivedSmallSource']=archived.get(item['executedSHA256']);assert item['archivedSmallSource']
    # 所有下载原文保留版权/许可证/NOTICE；这是来源证据，未作为新运行库。
    for p in (BASE/'research/sundials').rglob('*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    save(ev/'external_source.json',dict(repository='https://github.com/LLNL/sundials',tag='v7.4.0',commit='8e17876d3b4d682b4098684b07a85b005a122f81',
        license='BSD-3-Clause',retrievedVia='GitHub connector; tag sources compared with fixed commit sources',installed=False,linked=False,executed=False,
        files={p.relative_to(BASE/'research/sundials').as_posix():sha(p) for p in (BASE/'research/sundials').rglob('*') if p.is_file()}))
    assets={}
    for p in BASE.rglob('*'):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        h=sha(p);assets[p.relative_to(BASE).as_posix()]=dict(bytes=p.stat().st_size,sha256=h,archivedSmallEvidencePath=archived.get(h))
    save(ev/'local_asset_index.json',assets);save(ev/'execution_source_index.json',sourceIndex)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,baseCommit=baseline,
        originalFaces=19752,originalAreaM2=1.9646292565715553,originalNodes=335726,unknowns=1974844,
        currentProfileAssemblyRejected=True,reason='CR112180 not a detailed near-wall field; TP3271 simplified equation requires constant edge entropy',
        reusableNativeResidualImplemented=True,nonlinearIterations=0,timeUpdates=0,newHeatSolutions=0,actualProfilesFilled=0,
        physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,KINSOLExecuted=False,KINSOLInstalled=False,
        nativeSummary=native,independentAudit=audit,workbenchSummary=ui,allReceiptsTerminalStableGroupsEmpty=True,
        executionFailures=failures,processSeconds=totals,maximumSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in receipts.values()),
        historicalSourceBindings=changes,protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
        unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],scientificImagesVisuallyInspected=True,browserInteractionTested=False,
        localOnly=True,githubPushed=False,sharedWebGLChanged=False,uniqueDataDeleted=False,physicsObjectsRecompiled=0,
        fullFieldsArchivedInGit=False,nativeLibraryAndToolchainArchivedInGit=False,closedCLIs=closed,ownedFiles=OWN,
        currentOwnedFileSHA256={n:sha(ROOT/n) for n in OWN},developmentReadsArchivingSeparate=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputLogicalBytesBeforeFinalReport']<512*1024**2
    report['finalizationSeconds']=time.monotonic()-started;save(out/'report.json',report);save(ev/'final_report.json',report)
    save(ev/'manifest.json',{p.relative_to(ev).as_posix():sha(p) for p in ev.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ['passed','unknowns','nonlinearIterations','processSeconds','finalizationSeconds']}))


if __name__=='__main__':main()
