"""模块职责：保留失败和准确执行源、核对保护资产，将小证据交付Git而大场留E盘。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v240_validation';DOC=ROOT/'docs/validation/v240'
CLIS=['scripts/prepareForebodyEulerV240.py','scripts/runForebodyEulerV240.py','scripts/recoverForebodyMapV240.py','scripts/auditForebodyEulerV240.py']
OWN=CLIS+['scripts/clipEulerForebody.py','scripts/exportForebodyEulerV240.py','scripts/renderForebodyEulerV240.py',
 'scripts/localizeForebodyResidualV240.py','scripts/finalizeForebodyEulerV240.py','native/su2PhaseAwareAir/RiemannOutletObserverV240.hpp',
 'native/su2PhaseAwareAir/forebodyEulerV240.hpp','tests/testForebodyClipV240.py','tests/testForebodyWorkbenchV240.js',
 'docs/FOREBODY_PLANAR_EULER_V240.md','docs/validation/v240/PLAN.md','docs/validation/v240/REGIONAL_PLAN.md',
 'docs/validation/v240/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for p,h in protected.items():assert sha(ROOT/p)==h,p
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    assets=read(ROOT/'docs/validation/v222/evidence/final_report.json')['assets']
    for p,h in assets.items():assert sha(ROOT/p)==h,p
    untouched={'tests/verifyWebSectionExport.py':'6f67f9fecc3e806ed74ab879ff8a5ba1b745a239158babba55d1856d4d54225b',
      'web/section.js':'74fd8f6cd6f017ecb60a578ebedfe5aabfda2cde8ce5477fdd5b4b9bf24f9169',
      'web/viewer.js':'90d73959c1d39bf7fd0919d1a034e4994c24edd90a922270064a3c5892088ca5'}
    for p,h in untouched.items():assert sha(ROOT/p)==h,p
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')}
    failures={n:r for n,r in records.items() if not r['passed']}
    assert set(failures)=={'process_prepare','process_audit','process_audit2'}
    totals={};historical=[]
    manual={str(ROOT/'native/su2PhaseAwareAir/forebodyEulerV240.hpp'):BASE/'fix01/forebodyEulerV240.executed.hpp',
      str(DOC/'REGIONAL_PLAN.md'):BASE/'fix01/REGIONAL_PLAN.executed.md'}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            snapshot=manual.get(p,BASE/name/'execution_sources'/(h+Path(p).suffix))
            assert snapshot.exists() and sha(snapshot)==h,(name,p,h)
            historical.append(dict(process=name,path=p,executedSHA256=h,snapshot=snapshot.relative_to(BASE).as_posix(),currentSHA256=sha(p)))
    assert totals['offline']<=120 and totals['audit']<=240 and totals['build']<=120 and totals['query']<=15 and totals['cfd']<=600
    for folder in ('test01','prepare02','audit02','audit_accept01','workbench01','visual02','diagnostic01'):
        assert read(BASE/folder/'report.json')['passed'],folder
    a=read(BASE/'audit02/report.json');ui=read(BASE/'workbench01/report.json')
    assert a['completedUpdates']==200 and a['nativeFinalCSVBitEqual'] and a['independentStopPredicatePassed']
    assert a['correctedReadbackFlowUpdates']==0 and read(BASE/'process_run/report.json')['passed']
    assert not a['physicalMinus8Passed'] and not a['physicalQualification'] and a['newQualifiedAreaM2']==0
    for p in OWN:
        if p.endswith('.py'):ast.parse((ROOT/p).read_text(encoding='utf8'))
    closed={}
    for p in CLIS:
        result=subprocess.run([sys.executable,str(ROOT/p),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=8)
        message=result.stdout.decode('utf8',errors='replace');assert result.returncode!=0 and 'closed' in message,(p,message)
        closed[p]=message.strip()
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel=None):
        target=evidence/(rel if rel is not None else p.relative_to(BASE));target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    for p in list(BASE.rglob('*')):
        if not p.is_file() or out in p.parents:continue
        rel=p.relative_to(BASE)
        large=(p.suffix in ('.exe','.obj','.a','.npz','.su2') or p.name.startswith(('state_','restart_flow','original_faces_')) or
          p.name in ('forebody_research.json','case-data.js','index.html','saved_final.csv'))
        if large:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p));continue
        assert p.stat().st_size<4*1024**2,('UNCLASSIFIED_LARGE_FILE',str(p))
        copy(p)
    # 精确最终源与旧执行源分别归档，不能把已修复头标成原200更新使用的头。
    for p in OWN:
        rel=Path('final_source')/(p if p!='.gitattributes' else '.gitattributes.snapshot');copy(ROOT/p,rel)
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,completedUpdates=200,
      nativeOutletCalls=a['actualBoundaryCalls'],nativeOutletBitEqualComponents=a['actualBoundaryFiveComponentsBitEqual'],
      minimumNativeOutletNormalMach=a['minimumActualOutletNormalMach'],nativeFinalCSVBitEqual=a['nativeFinalCSVBitEqual'],
      originalMappingFailurePreserved=True,correctedNativeMappingReadbackPassed=True,correctedMappingReadbackFlowUpdates=0,
      originalPressureChange=a['pressureTrendCorrection'][-1]['originalIncorrectNativeChange'],correctedPressureChange=a['finalPressureAreaChange'],
      physicalFinalLogRMS=a['physicalFinalLogRMS'],finalNetMassEnergyScaled=a['finalNetMassEnergyScaled'],
      researchRegionFaces=4434,researchRegionAreaM2=a['regionAreaM2'],originalFaces=19752,originalAreaM2=a['originalAreaM2'],
      workbenchStages=ui['updates'],workbenchFields=ui['fields'],exactResearchValues=ui['exactResearchValues'],
      newHeatAreaM2=0,newQualifiedAreaM2=0,physicalQualification=False,
      protectedCount=len(protected),engineeringArtifactsUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
      assets=assets,unrelatedTrackedSHA256=untouched,processSeconds=totals,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
      CFDPeakSampledRSSBytes=records['process_run']['peakObservedGroupRSSBytes'],
      allAvailableTerminalReceiptsInputsStable=True,allAvailableTerminalReceiptsGroupsEmpty=True,
      failedProcesses={n:dict(exitCode=r['exitCode'],stopReason=r['stopReason']) for n,r in failures.items()},
      historicalSourceBindings=historical,closedCLIs=closed,MATLABExecuted=False,sharedWebChanged=False,browserInteractionTested=False,
      scientificImageVisuallyInspected=True,originalViscousMeshChanged=False,uniqueAssetsDeleted=False,
      unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    report['outputBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputBytesBeforeFinalReport']<256*1024**2
    report['finalizationSeconds']=time.monotonic()-start
    save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    paths=OWN+['docs/validation/v240/evidence_manifest.json']+sorted(p.relative_to(ROOT).as_posix() for p in evidence.rglob('*') if p.is_file())
    assert len(paths)==len(set(paths));save(out/'delivery_files.json',paths)
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=report['finalizationSeconds'],processSeconds=totals,
      outputBytes=report['outputBytesBeforeFinalReport'],failedProcesses=report['failedProcesses'])))


if __name__=='__main__':main()
