"""模块职责：归档真实区域证据，保留失败及并发文档的本地补充。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v244_validation';DOC=ROOT/'docs/validation/v244'
CLIS=['scripts/'+n+'.py' for n in ('prepareViscousRegionV244','auditViscousRegionV244',
    'reviewViscousCouplingV244','renderViscousRegionV244')]
OWN=CLIS+['scripts/finalizeViscousRegionV244.py','docs/VISCOUS_REGIONAL_BOUNDARY_V244.md',
    'docs/validation/v244/PLAN.md','docs/validation/v244/NEXT_DECISION.md',
    'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
SPLIT_DOCS=('AGENTS.md','docs/PROJECT_GOALS.md')


def source(name):
    return BASE/'owned_source'/name if name in SPLIT_DOCS else ROOT/name


def main():
    """模块一：保护原资产，确认并发补充仍在工作区，交付稿仅含本任务修改。"""
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    prev=read(ROOT/'docs/validation/v243/evidence/final_report.json')
    assets=prev['assets'];unrelated=prev['unrelatedTrackedSHA256']
    for n,h in dict(assets,**unrelated).items():assert sha(ROOT/n)==h,n
    private=read(BASE/'concurrent_prefixes.json');split={}
    for n in SPLIT_DOCS:
        text=(ROOT/n).read_text(encoding='utf8');prefix=private[n]['prefix']
        assert text.count(prefix)==1 and text.replace(prefix,'',1)==source(n).read_text(encoding='utf8'),n
        assert 'public_benchmarks_20261116' not in source(n).read_text(encoding='utf8')
        split[n]=dict(pendingLocalAdditionPreserved=True,publishedOwnVersionSHA256=sha(source(n)),
            worktreeCombinedSHA256=sha(ROOT/n),publishedOwnSource=source(n).relative_to(ROOT).as_posix())
    records={p.parent.name:read(p) for p in sorted(BASE.glob('process_*/report.json'))}
    assert len(records)==6;failures={'process_prepare':1,'process_audit':15}
    totals={};bindings=[];required=set();cache={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert r['exitCode']==failures.get(name,0) and r['passed']==(name not in failures)
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if p not in cache:cache[p]=sha(p)
            if cache[p]==h:continue
            matches=list(BASE.glob('process_*/execution_sources/'+h+Path(p).suffix))
            assert matches,(name,p);snapshot=matches[0];assert sha(snapshot)==h;required.add(snapshot)
            bindings.append(dict(process=name,path=p,executedSHA256=h,currentSHA256=cache[p],
                snapshot=snapshot.relative_to(BASE).as_posix()))
    assert totals['offline']<120 and totals['audit']<240
    for name in ('region02','audit02','coupling01','visual01'):
        assert read(BASE/name/'report.json')['passed'],name
    audit=read(BASE/'audit02/report.json');fixed=read(BASE/'audit02/corrected_region_report.json')
    assert audit['exactTopologicalWallContactFaces']==46
    assert audit['strictlyZeroStoredDistanceFaces']==14 and fixed['facesRequiringDownstreamInformation']==827
    assert audit['originalStateScalarsBitEqual']==2602700 and not fixed['physicalQualification']
    closed={}
    for name in OWN:
        if name.endswith('.py'):ast.parse(source(name).read_text(encoding='utf8'))
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,timeout=10)
        message=p.stdout.decode('utf8',errors='replace')
        assert p.returncode!=0 and 'closed' in message,(name,message);closed[name]=message.strip()

    """模块二：完整保留小回执、错误、接口CSV、源摘录和图，大数组只登记本地。"""
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        dest=evidence/(rel if rel is not None else p.relative_to(BASE))
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=p in required or p.suffix in ('.json','.txt','.png','.csv') or top.startswith('process_')
        if p.name=='update_delivery_docs.py':keep=True
        if p.suffix in ('.npz','.su2','.bin','.pyc') or top=='owned_source':keep=False
        if p.name=='concurrent_prefixes.json':keep=False
        if keep:
            assert p.stat().st_size<2*1024**2,('LARGE_EVIDENCE',str(p));copy(p)
        else:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for name in OWN:
        copy(source(name),Path('final_source')/(name if name!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,
        originalFaces=19752,originalAreaM2=1.9646292565715553,regionalNodes=76550,regionalCells=160304,
        regionWallParents=4434,internalInterfaces=1984,interfacesRequiringDownstreamInformation=827,
        originalStateScalarsBitEqual=2602700,independentMaxScaledError=audit['independentThermodynamicAndNormalMaxScaledError'],
        exactWallContactInterfaces=46,strictZeroDistanceInterfaces=14,firstIncorrectMetadataPreserved=True,
        correctedMetadataFunctionActuallyTested=True,sourceState='v236 update 300 research',
        physicalCutProductionRejected=True,physicalClosure='same native physical zone retaining downstream',
        newFlowUpdates=0,newBuilds=0,newNativeQueries=0,fasCandidateOnly=True,mpiInstalledOrBuilt=False,
        physicalQualification=False,edgeAccepted=False,newHeatAreaM2=0,newQualifiedAreaM2=0,
        originalNearWallAssetsChanged=False,sharedWebGLChanged=False,MATLABExecuted=False,
        browserInteractionTested=False,scientificImageVisuallyInspected=True,
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
        assets=assets,unrelatedTrackedSHA256=unrelated,splitDocumentationOwnership=split,
        concurrentLocalOnlyDataPublished=False,allReceiptsTerminalStableGroupsEmpty=True,
        executionFailuresPreserved=failures,firstAuditTimedOut=True,auditRecoveryNoNumericalChanges=True,
        processSeconds=totals,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
        historicalSourceBindings=bindings,closedCLIs=closed,uniqueAssetsDeleted=False,
        unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    report['outputBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputBytesBeforeFinalReport']<256*1024**2
    report['finalizationSeconds']=time.monotonic()-start;save(out/'report.json',report)
    copy(out/'report.json',Path('final_report.json'))
    digest={p.relative_to(DOC).as_posix():sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()}
    save(DOC/'evidence_manifest.json',digest)
    files=OWN+['docs/validation/v244/evidence_manifest.json']
    files += [p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(json.dumps(dict(passed=True,ownFiles=len(files),evidenceFiles=len(digest),
        evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),
        processSeconds=totals,finalizationSeconds=report['finalizationSeconds'])))


if __name__=='__main__':main()
