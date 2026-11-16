"""模块职责：核对保护资产，归档准确执行源及小证据，不重跑数值任务。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha

BASE=ROOT/'outputs/v243_validation'
DOC=ROOT/'docs/validation/v243'
CLIS=['scripts/'+name+'.py' for name in (
    'inspectForebodyTransportV243','refineForebodyV243',
    'reviewRefinementDiagonalV243','reconnectRefinementV243',
    'runForebodySpatialV243','auditForebodySpatialV243',
    'exportForebodySpatialV243','renderForebodySpatialV243')]
OWN=CLIS+[
    'scripts/finalizeForebodySpatialV243.py','tests/testForebodyWorkbenchV243.js',
    'docs/FOREBODY_SPATIAL_COMPARISON_V243.md',
    'docs/validation/v243/PLAN.md','docs/validation/v243/METHOD_DECISION.md',
    'docs/validation/v243/CFD_PLAN.md','docs/validation/v243/NEXT_DECISION.md',
    'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
OWN += [p.relative_to(ROOT).as_posix() for p in sorted((DOC/'research').rglob('*')) if p.is_file()]
EXTRA_UNRELATED={
    'docs/ENGINEERING_ALGORITHMS_REPORT_2026-11-12.md':
        'a108ba026a6767cd9179d13ed31196386b75c6f9f0d0ab192bc52b52ba91ca07',
    'scripts/buildEngineeringAlgorithmReport.py':
        '7e233eb57e850dd948d0a7b50a2c7a330030177c708e4ee3d56b9c919f4eaab0'}


def main():
    """模块一：先完整预检保护资产、终态及历史源绑定，再创建交付目录。"""
    start=time.monotonic()
    out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    previous=read(ROOT/'docs/validation/v242/evidence/final_report.json')
    assets=previous['assets'];unrelated=previous['unrelatedTrackedSHA256']
    for name,h in dict(assets,**unrelated).items():assert sha(ROOT/name)==h,name
    # 并发报告文件不归本任务；记录观察，不把另一任务的变化当成本轮破坏。
    concurrent_changes={name:dict(initialObservedSHA256=h,currentSHA256=sha(ROOT/name))
        for name,h in EXTRA_UNRELATED.items()}
    concurrent_tracked={}
    changed=subprocess.check_output(['git','diff','--name-only','-z'],cwd=str(ROOT))
    for name in changed.decode('utf8').split('\0'):
        if name and name not in OWN and name not in unrelated:
            concurrent_tracked[name]=sha(ROOT/name)
    records={p.parent.name:read(p) for p in sorted(BASE.glob('process_*/report.json'))}
    assert len(records)==14
    failures={'process_inspect':1,'process_flow_audit':1}
    totals={};bindings=[];required_snapshots=set();current_hashes={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert r['exitCode']==failures.get(name,0) and r['passed']==(name not in failures),name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for path,h in r['inputSHA256'].items():
            if path not in current_hashes:current_hashes[path]=sha(path)
            if current_hashes[path]==h:continue
            paths=list(BASE.glob('inputs_*/'+h+Path(path).suffix))
            paths+=list(BASE.glob('process_*/execution_sources/'+h+Path(path).suffix))
            assert paths,(name,path)
            snapshot=paths[0];assert sha(snapshot)==h;required_snapshots.add(snapshot)
            bindings.append(dict(process=name,path=path,executedSHA256=h,
                snapshot=snapshot.relative_to(BASE).as_posix(),currentSHA256=current_hashes[path]))
    for kind,limit in dict(offline=120,audit=240,mesh=60,build=60,cfd=600).items():
        assert totals[kind]<limit,(kind,totals[kind])
    for name in ('inspect02','refine01','diagonal01',
                 'reconnected_audit01','flow_prepare01','flow_audit02','workbench01','visual01'):
        assert read(BASE/name/'report.json')['passed'],name
    assert (BASE/'reconnected01/local_proof.json').is_file()
    rejected=read(BASE/'refine_audit01/report.json')
    assert not rejected['passed'] and rejected['geometryPassed'] and not rejected['qualityPassed']
    flow=read(BASE/'flow_audit02/report.json')
    mesh=read(BASE/'reconnected_audit01/report.json')
    ui=read(BASE/'workbench01/report.json')
    assert flow['completedUpdates']==800 and not flow['benefitStop']
    assert not flow['physicalMinus8Passed'] and not flow['netMassEnergyPassed']
    assert flow['pressureStabilityPassed'] and not flow['physicalQualification']
    assert not flow['spatialComparison']['gridIndependenceProved']
    assert mesh['cells']==162216 and mesh['nodes']==34956 and mesh['changedCentralCells']==16
    assert ui['exactResearchValues']==79812 and ui['fields']==21
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    closed={}
    for name in CLIS:
        r=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        message=r.stdout.decode('utf8',errors='replace')
        assert r.returncode!=0 and 'closed' in message,(name,message)
        closed[name]=message.strip()

    """模块二：保留失败和实际源字节；大体场及重复载荷只登记本地SHA。"""
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        dest=evidence/(rel if rel is not None else p.relative_to(BASE))
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=p in required_snapshots or p.suffix in ('.png','.patch')
        if top.startswith('process_'):keep=True
        if top.startswith('archive_preflight_'):keep=True
        if top.startswith('inputs_') and p.suffix in ('.py','.md') and p.stat().st_size<200000:keep=True
        if top=='flow_prepare01':keep=True
        if top in ('inspect01','inspect02') and p.suffix in ('.csv','.txt'):keep=True
        if top=='run01':
            keep=p.name in ('case.cfg','native_completion.json','benefit.csv',
                'pressure_trend.csv','outlet_history.csv','physical_history.csv','history.csv')
        if top=='flow_audit02':keep=p.name in ('report.json','coarse_fine_original_faces.csv')
        if p.suffix=='.json' and p.name!='forebody_research.json':keep=True
        if p.name in ('build.txt','frontend_test.txt','document_command_failure.txt','update_delivery_docs.py'):keep=True
        if p.suffix in ('.npz','.exe','.obj','.a','.su2','.dat','.pyc'):keep=False
        if p.name in ('case-data.js','forebody_research.json','index.html'):keep=False
        if keep:
            assert p.stat().st_size<2*1024**2,('LARGE_EVIDENCE',str(p));copy(p)
        else:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for name in OWN:
        copy(ROOT/name,Path('final_source')/(name if name!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)

    """模块三：科学拒绝与软件通过分开，生成明确白名单供Git交付。"""
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,
        CFDUpdates=800,stopReason=flow['stopReason'],thinBuilds=1,matureNumericalObjectsRecompiled=0,
        sourceAnalysisActualCells=210,sourceAnalysisStoppedWallFaces=6,
        rtResearchOnly=True,rtInstalledOrExecuted=False,gmshRefineCalls=1,
        rawRefinementQualityRejected=True,finiteCustomGeometryAdapterCells=16,
        meshNodes=mesh['nodes'],meshCells=mesh['cells'],minimumIndependentSICN=mesh['minimumIndependentSICN'],
        originalFaces=19752,originalAreaM2=1.9646292565715553,
        originalWallParents=4434,computationalChildWalls=17736,originalViscousMeshChanged=False,
        physicalFinalLogRMS=flow['physicalFinalLogRMS'],netMassEnergyPassed=False,
        finalNetMassEnergyScaled=flow['finalNetMassEnergyScaled'],pressureStabilityPassed=True,
        spatialComparison=flow['spatialComparison'],physicalQualification=False,edgeAccepted=False,
        newHeatAreaM2=0,newQualifiedAreaM2=0,UIResearchValues=ui['exactResearchValues'],UIFields=ui['fields'],
        savedUpdates=flow['savedUpdates'],nativeFinalCSVBitEqual=flow['nativeFinalCSVBitEqual'],
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),
        candidateCount=len(candidates['files']),assets=assets,unrelatedTrackedSHA256=unrelated,
        concurrentTrackedSHA256AtFinalization=concurrent_tracked,
        concurrentObservedChanges=concurrent_changes,unrelatedChangesExcludedFromDelivery=True,
        processSeconds=totals,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
        CFDPeakSampledRSSBytes=records['process_flow_run']['peakObservedGroupRSSBytes'],
        allTerminalReceiptsInputsStable=True,allTerminalReceiptsGroupsEmpty=True,
        executionFailuresPreserved=failures,scientificRejectionsPreserved=['refine_audit01','flow_audit02 physical gates'],
        administrativeCommandFailure='file invocation reported source encoding SyntaxError; explicit UTF-8 decode/compile/exec recovered; cause unproved',
        archivePreflightFailure='concurrent unrelated report script changed; initial source and refusal preserved; ownership-aware recovery',
        historicalSourceBindings=bindings,closedCLIs=closed,sharedWebGLChanged=False,
        MATLABExecuted=False,browserInteractionTested=False,scientificImagesVisuallyInspected=True,
        uniqueAssetsDeleted=False,unmeteredDevelopmentTargetedReadsAndNetworkDisclosed=True)
    report['outputBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputBytesBeforeFinalReport']<384*1024**2
    report['finalizationSeconds']=time.monotonic()-start
    save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    digest={p.relative_to(DOC).as_posix():sha(p) for p in sorted(evidence.rglob('*')) if p.is_file()}
    save(DOC/'evidence_manifest.json',digest)
    files=OWN+['docs/validation/v243/evidence_manifest.json']
    files += [p.relative_to(ROOT).as_posix() for p in sorted(evidence.rglob('*')) if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(json.dumps(dict(passed=True,files=len(files),evidenceFiles=len(digest),
        evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),
        outputBytes=report['outputBytesBeforeFinalReport'],processSeconds=totals,
        finalizationSeconds=report['finalizationSeconds'])))


if __name__=='__main__':main()
