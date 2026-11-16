"""模块职责：绑定准确执行源/失败/保护资产，把本版小证据交付GitHub；不启动数值任务。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha

BASE=ROOT/'outputs/v241_validation';DOC=ROOT/'docs/validation/v241'
CLIS=['scripts/'+name+'.py' for name in ['inspectForebodyBufferV241','runStreamtubeV241','produceForebodyStreamtubesV241',
       'auditForebodyStreamtubesV241','compareStreamtubeRecoveryV241','exportForebodyStreamtubesV241']]
OWN=CLIS+['scripts/tetraStreamtubes.py','scripts/finalizeStreamtubesV241.py','tests/testTetraStreamtubesV241.py',
 'tests/testStreamtubeWorkbenchV241.js','web/streamtube-research-v241.html','web/streamtube-research-v241.js',
 'docs/INTERNAL_VOLUME_STREAMTUBES_V241.md','docs/validation/v241/PLAN.md','docs/validation/v241/STREAMTUBE_PLAN.md',
 'docs/validation/v241/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    tick=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidate_manifest=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidate_manifest['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    assets=read(ROOT/'docs/validation/v240/evidence/final_report.json')['assets']
    for name,h in assets.items():assert sha(ROOT/name)==h,name
    unrelated={'tests/verifyWebSectionExport.py':'6f67f9fecc3e806ed74ab879ff8a5ba1b745a239158babba55d1856d4d54225b',
       'web/section.js':'74fd8f6cd6f017ecb60a578ebedfe5aabfda2cde8ce5477fdd5b4b9bf24f9169',
       'web/viewer.js':'90d73959c1d39bf7fd0919d1a034e4994c24edd90a922270064a3c5892088ca5'}
    for name,h in unrelated.items():assert sha(ROOT/name)==h,name
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')}
    assert len(records)==10;failed={name:r for name,r in records.items() if not r['passed']}
    assert set(failed)=={'process_test','process_audit'}
    totals={};bindings=[]
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            candidates=list(BASE.glob('inputs_*/'+h+Path(p).suffix))+list(BASE.glob('process_*/execution_sources/'+h+Path(p).suffix))
            assert candidates,(name,p,h)
            snapshot=candidates[0];assert sha(snapshot)==h
            bindings.append(dict(process=name,path=p,executedSHA256=h,snapshot=snapshot.relative_to(BASE).as_posix(),currentSHA256=sha(p)))
    assert set(totals)=={'offline','audit'} and totals['offline']<120 and totals['audit']<240
    for name in ['buffer01','test03','audit02','recovery01','workbench01']:assert read(BASE/name/'report.json')['passed'],name
    manifest=read(BASE/'streamtubes02/manifest.json');audit=read(BASE/'audit02/report.json');ui=read(BASE/'workbench01/report.json')
    assert manifest['passed'] and not manifest['edgeAccepted'] and not manifest['physicalQualification']
    assert audit['savedRayStationsAudited']==6748 and ui['exactCSVRoundtripNumbers']==9746
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    closed={}
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=8)
        message=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in message,(name,message)
        closed[name]=message.strip()
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        target=evidence/(rel if rel is not None else p.relative_to(BASE));target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    for p in list(BASE.rglob('*')):
        if not p.is_file() or out in p.parents:continue
        large=p.suffix in ('.npz','.exe','.obj','.a','.su2') or p.name in ('streamtubes.json','streamtubes-data.js','index.html')
        if large:local[p.relative_to(BASE).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p));continue
        assert p.stat().st_size<4*1024**2,('unclassified-large-file',str(p));copy(p)
    for name in OWN:copy(ROOT/name,Path('final_source')/(name if name!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,CFDUpdates=0,nativeBuilds=0,nativeQueries=0,
      actualResearchTubes=4,commonStations=443,auditedRayStations=6748,exactCSVValues=9746,
      originalFaces=19752,originalAreaM2=1.9646292565715553,originalViscousMeshChanged=False,
      newHeatAreaM2=0,newQualifiedAreaM2=0,physicalQualification=False,edgeAccepted=False,
      fixedBufferRejected=True,impossibleFixedTriangleCount=13,
      failedTestDerivativeAssertionPreserved=True,originalTube2BoundaryCrossingPreserved=True,
      recoveredBoundaryStops={'wall':43,'outlet':2},allSavedNineRaySegmentsAudited=True,
      protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidate_manifest['files']),
      assets=assets,unrelatedTrackedSHA256=unrelated,processSeconds=totals,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      allAvailableTerminalReceiptsInputsStable=True,allAvailableTerminalReceiptsGroupsEmpty=True,
      failedProcesses={n:dict(exitCode=r['exitCode'],stopReason=r['stopReason']) for n,r in failed.items()},
      historicalSourceBindings=bindings,closedCLIs=closed,sharedWebGLChanged=False,MATLABExecuted=False,
      browserInteractionTested=False,scientificImageVisuallyInspected=True,uniqueAssetsDeleted=False,
      unmeteredDevelopmentAndTargetedReadsDisclosed=True,pauseResumeNoActiveNumericalProcessObserved=True)
    report['outputBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputBytesBeforeFinalReport']<256*1024**2
    report['finalizationSeconds']=time.monotonic()-tick
    save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    digest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()};save(DOC/'evidence_manifest.json',digest)
    files=OWN+['docs/validation/v241/evidence_manifest.json']+[p.relative_to(ROOT).as_posix() for p in evidence.rglob('*') if p.is_file()]
    assert len(files)==len(set(files));save(out/'delivery_files.json',files)
    print(json.dumps(dict(passed=True,files=len(files),evidenceFiles=len(digest),processSeconds=totals,outputBytes=report['outputBytesBeforeFinalReport'],finalizationSeconds=report['finalizationSeconds'])))


if __name__=='__main__':main()
