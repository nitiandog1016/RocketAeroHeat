"""模块职责：核对保护资产和准确执行源，归档失败及必要证据；大场只建本地索引。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v238_validation';DOC=ROOT/'docs/validation/v238'
CLIS=['scripts/runHLLCCompatibilityV238.py','scripts/runHLLCNativeV238.py','scripts/runHLLCRocketV238.py',
 'scripts/auditHLLCInterfaceV238.py','scripts/auditHLLCNativeV238.py','scripts/auditHLLCRocketV238.py','scripts/auditHLLCAcceptanceV238.py']
OWN=['native/su2PhaseAwareAir/hllcInterfaceV238.cpp','native/su2PhaseAwareAir/PhaseStepObservationV238.hpp',
 'native/su2PhaseAwareAir/uniformHLLCBranchV238.cpp']+CLIS+[
 'scripts/exportHLLCRocketV238.py','scripts/renderHLLCRocketV238.py','scripts/reviewHLLCOutcomeV238.py',
 'scripts/finalizeHLLCV238.py','tests/testRestartBranchV238.js','docs/GENERAL_HLLC_COMPATIBILITY_V238.md',
 'docs/validation/v238/PLAN.md','docs/validation/v238/ROCKET_PLAN.md','docs/validation/v238/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence'
    assert not out.exists() and not evidence.exists()
    # 模块一：既有保护清单与独有资产不受本轮源码/交付影响。
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    eng=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in eng['files']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    assets=read(ROOT/'docs/validation/v222/evidence/final_report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};totals={};historical=[]
    assert len(records)==23
    assert {n for n,r in records.items() if not r['passed']}=={'process_audit01','process_outcome'}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            snapshot=BASE/name/'execution_sources'/(h+Path(p).suffix)
            assert snapshot.exists() and sha(snapshot)==h,(name,p,h)
            historical.append(dict(process=name,path=p,executedSHA256=h,snapshot=snapshot.relative_to(BASE).as_posix(),currentSHA256=sha(p)))
    for folder in ('audit02','native_audit01','real_audit01','acceptance_audit01','workbench01','visual01','outcome02'):
        assert read(BASE/folder/'report.json')['passed'],folder
    audit=read(BASE/'real_audit01/report.json');accepted=read(BASE/'acceptance_audit01/report.json');ui=read(BASE/'workbench01/report.json')
    assert audit['successfulUpdates']==200 and accepted['actualAcceptedEvents']==16 and accepted['nativeFinalSixStateCoordinatesBitEqual']
    # 模块二：关闭旧准入并实际检验拒绝；不改变已经归档的执行源或回执。
    for n in OWN:
        if n.endswith('.py'):ast.parse((ROOT/n).read_text(encoding='utf8'))
    closed={}
    for n in CLIS:
        r=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
        value=r.stdout.decode('utf8',errors='replace');assert r.returncode!=0 and 'closed' in value,(n,value)
        closed[n]=value.strip()
    out.mkdir();evidence.mkdir()
    def copy(p,relative=None):
        target=evidence/(relative or p.relative_to(BASE));target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    # 模块三：必要小源、控制、观察、回执全部保留；大数据不重复放入Git。
    local={}
    for p in list(BASE.rglob('*')):
        if not p.is_file() or out in p.parents:continue
        rel=p.relative_to(BASE)
        large=(p.suffix in ('.exe','.obj','.a') or (p.suffix=='.bin' and rel.parts[0]=='rocket01') or
               p.name.startswith('original_faces_') or (rel.parts[0]=='rocket01' and p.name.startswith(('wall_000','restart_flow'))) or
               p.name in ('rocket_restart_comparison.json','case-data.js','index.html'))
        if large:
            local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p));continue
        assert p.stat().st_size<4*1024**2,('UNCLASSIFIED_LARGE_FILE',str(p),p.stat().st_size)
        copy(p)
    save(evidence/'local_asset_index.json',dict(base='outputs/v238_validation',files=local,largeAssetsBackedUpToGitHub=False,
        requiredRuntime='v238 HLLC density-pressure derivative fix/config/failure observer; v237 guard; v235 flow-SA observer; v220/v167 exact objects; new budget required'))
    for n in OWN:copy(ROOT/n,Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    elapsed=time.monotonic()-start
    assert totals['offline']+elapsed<120 and totals['build']<180 and totals['query']<60 and totals['audit']<240 and totals['cfd']<2405
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<1024**3
    report=dict(passed=True,goalStatus='active',physicalQualification=False,newQualifiedAreaM2=0,newHeatAreaM2=0,
      realRocketUpdates=200,originalFaces=19752,originalAreaM2=1.9646292565715553,meshCells=705490,points=335726,
      numericalObjectsCompiled=3,originalHLLCControlFailed=True,correctedControl=read(BASE/'audit02/report.json'),
      acceptance=accepted,realBranch=audit['nativeCompletion'],physicalMinus8Passed=audit['physicalMinus8All'],
      workbenchStages=ui['updates'],workbenchFields=ui['fields'],protectedCount=len(protected),engineeringArtifactsUnchanged=len(eng['files']),
      candidateCount=len(candidates['files']),assets=assets,processSeconds=totals,finalizationSeconds=elapsed,
      maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),outputBytes=size,
      allRecordedInputsStable=True,allRecordedGroupsEmpty=True,historicalSourceBindings=historical,closedCLIs=closed,
      failedProcesses={n:{k:r[k] for k in ('exitCode','stopReason','elapsedSeconds')} for n,r in records.items() if not r['passed']},
      MATLABExecuted=False,sharedWebChanged=False,browserInteractionTested=False,scientificImageVisuallyInspected=True,
      uniqueAssetsDeleted=False,newRocketMeshGenerated=False,unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v238/evidence_manifest.json']+['docs/validation/v238/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,outputBytes=size,processSeconds=totals)))


if __name__=='__main__':main()
