"""模块职责：保护现有资产，归档准确执行源和小证据；不重新执行数值任务。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v242_validation';DOC=ROOT/'docs/validation/v242'
CLIS=['scripts/'+name+'.py' for name in ('prepareForebodyPLCV242','surfaceForebodyPLCV242','buildForebodyPLCV242',
 'meshForebodyPLCV242','auditForebodyPLCV242','runForebodyPLCV242','runForebodyEulerV242','auditForebodyEulerV242',
 'exportForebodyEulerV242','renderForebodyEulerV242','compareForebodyTubesV242')]
OWN=CLIS+['scripts/finalizeForebodyPLCV242.py','tests/testForebodyWorkbenchV242.js',
 'third_party/tetgen/single_shell_plc_v242.cpp','docs/FOREBODY_FIXED_WALL_PLC_V242.md',
 'docs/validation/v242/PLAN.md','docs/validation/v242/CFD_PLAN.md','docs/validation/v242/STREAMTUBE_PLAN.md',
 'docs/validation/v242/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    previous=read(ROOT/'docs/validation/v241/evidence/final_report.json');assets=previous['assets'];unrelated=previous['unrelatedTrackedSHA256']
    for name,h in dict(assets,**unrelated).items():assert sha(ROOT/name)==h,name
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==13
    totals={};bindings=[];required_snapshots=set()
    for name,r in records.items():
        assert r['passed'] and r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for path,h in r['inputSHA256'].items():
            if sha(path)==h:continue
            paths=list(BASE.glob('inputs_*/'+h+Path(path).suffix))+list(BASE.glob('process_*/execution_sources/'+h+Path(path).suffix))
            assert paths,(name,path);snapshot=paths[0];assert sha(snapshot)==h;required_snapshots.add(snapshot)
            bindings.append(dict(process=name,path=path,executedSHA256=h,snapshot=snapshot.relative_to(BASE).as_posix(),currentSHA256=sha(path)))
    for key,limit in dict(offline=120,audit=240,cfd=600,build=60,**{'build-flow':60,'mesh-surface':120,'mesh-volume':120}).items():assert totals[key]<limit
    for name in ('prepare01','surface01','tetgen01','mesh_audit01','flow_prepare01','flow_audit01','workbench01','visual01','tube_audit01'):
        assert read(BASE/name/'report.json')['passed'],name
    flow=read(BASE/'flow_audit01/report.json');mesh=read(BASE/'mesh_audit01/report.json');tube=read(BASE/'tube_audit01/report.json');ui=read(BASE/'workbench01/report.json')
    assert flow['completedUpdates']==400 and not flow['benefitStop'] and not flow['physicalMinus8Passed'] and flow['netMassEnergyPassed']
    assert not flow['physicalQualification'] and tube['samePhysicalNineSeeds'] and not tube['edgeAccepted']
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    closed={}
    for name in CLIS:
        result=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        message=result.stdout.decode('utf8',errors='replace');assert result.returncode!=0 and 'closed' in message,(name,message)
        closed[name]=message.strip()
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        dest=evidence/(rel if rel is not None else p.relative_to(BASE));dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    # 大体场/网格/执行文件和重复网页载荷只登记SHA；保存小回执、实际驱动和科研图。
    for p in list(BASE.rglob('*')):
        if not p.is_file() or out in p.parents:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=p in required_snapshots or p.suffix in ('.png','.patch') or p.name in ('document_command_failure.txt','update_delivery_docs.py')
        if top.startswith('process_'):keep=True
        if top.startswith('inputs_') and p.suffix in ('.py','.md') and p.stat().st_size<200000:keep=True
        if top=='flow_prepare01':keep=True
        if top=='run01':keep=p.name in ('case.cfg','native_completion.json','benefit.csv','pressure_trend.csv','outlet_history.csv','physical_history.csv','history.csv')
        if top=='flow_audit01':keep=p.name in ('report.json','original_faces_000400.csv')
        if top=='tube_audit01':keep=True
        if p.suffix=='.json' and p.name!='forebody_research.json':keep=True
        if p.name in ('build.txt','tetgen.txt','frontend_test.txt'):keep=True
        if p.suffix in ('.npz','.exe','.obj','.a','.su2','.dat') or p.name in ('case-data.js','forebody_research.json','index.html'):keep=False
        if keep:
            assert p.stat().st_size<2*1024**2,('LARGE_EVIDENCE',str(p));copy(p)
        else:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for name in OWN:copy(ROOT/name,Path('final_source')/(name if name!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,CFDUpdates=400,
      thinBuilds=2,matureNumericalObjectsRecompiled=0,meshNodes=5354,meshCells=20277,minimumIndependentSICN=mesh['independentMinSICN'],
      originalFaces=19752,originalAreaM2=1.9646292565715553,originalWallFacesPreserved=4434,originalViscousMeshChanged=False,
      physicalFinalLogRMS=flow['physicalFinalLogRMS'],netMassEnergyPassed=True,pressureStabilityPassed=True,
      physicalQualification=False,edgeAccepted=False,newHeatAreaM2=0,newQualifiedAreaM2=0,
      actualResearchTubes=4,auditedRayStations=tube['savedRayStationsAudited'],commonStations=tube['totalCommonStations'],
      tubeCSVRoundtripValues=tube['exactCSVRoundtripValues'],UIResearchValues=ui['exactResearchValues'],UIFields=ui['fields'],
      protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      assets=assets,unrelatedTrackedSHA256=unrelated,processSeconds=totals,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      allTerminalReceiptsInputsStable=True,allTerminalReceiptsGroupsEmpty=True,numericalFailures={},
      administrativeCommandFailure='stdin encoding parse error, no changes executed; recovered from saved UTF-8 script',
      historicalSourceBindings=bindings,closedCLIs=closed,sharedWebGLChanged=False,MATLABExecuted=False,browserInteractionTested=False,
      scientificImagesVisuallyInspected=True,uniqueAssetsDeleted=False,unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    report['outputBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputBytesBeforeFinalReport']<256*1024**2
    report['finalizationSeconds']=time.monotonic()-start
    save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    digest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()};save(DOC/'evidence_manifest.json',digest)
    files=OWN+['docs/validation/v242/evidence_manifest.json']+[p.relative_to(ROOT).as_posix() for p in evidence.rglob('*') if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(json.dumps(dict(passed=True,files=len(files),evidenceFiles=len(digest),evidenceBytes=sum(p.stat().st_size for p in evidence.rglob('*') if p.is_file()),
      processSeconds=totals,finalizationSeconds=report['finalizationSeconds'])))


if __name__=='__main__':main()
