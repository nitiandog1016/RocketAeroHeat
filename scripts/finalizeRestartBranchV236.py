"""模块职责：归档真实运行、离线失败和原面出口，核对保护资产；不运行新的数值任务。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v236_validation';DOC=ROOT/'docs/validation/v236'
OWN=['native/su2PhaseAwareAir/restartBranchV236.cpp','scripts/runRestartBranchV236.py','scripts/auditRestartBranchV236.py',
 'scripts/exportRestartBranchV236.py','scripts/reviewRestartOutcomeV236.py','scripts/renderRestartBranchV236.py',
 'scripts/finalizeRestartBranchV236.py','tests/testRestartBranchV236.js','docs/RESTART_BRANCH_V236.md',
 'docs/validation/v236/PLAN.md','docs/validation/v236/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


def main():
    start=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
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
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==9
    failed={'process_audit':(1,None),'process_audit2':(15,'output')};historical=[];totals={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        code,reason=failed.get(name,(0,None));assert r['exitCode']==code and r['stopReason']==reason and r['passed']==(name not in failed),name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            snapshot=BASE/name/'execution_sources'/(h+Path(p).suffix);assert sha(snapshot)==h,(name,p)
            historical.append(dict(process=name,path=p,executedSHA256=h,currentSHA256=sha(p)))
    bound=read(BASE/'build01/inputs.json')
    for n,h in bound.items():assert sha(n)==h,n
    audit=read(BASE/'audit03/report.json');extra=read(BASE/'outcome01/report.json');ui=read(BASE/'workbench01/report.json');plot=read(BASE/'visual02/report.json')
    assert all(x['passed'] for x in (audit,extra,ui,plot))
    assert audit['successfulUpdates']==300 and audit['newRejectionEvents']==466 and extra['uniqueRejectedNodes']==19
    assert audit['benefitDecisions']==[dict(checkpoint=200,stop=False,reason='continue'),dict(checkpoint=300,stop=True,reason='significant-deterioration')]
    assert extra['nativeFinalRestartMissing'] and not extra['plannedNativeRestartAcceptance']
    assert ui['fields']==34 and not ui['sharedGUIChanged'] and ui['dataRoundTripExact']
    closed={}
    for n in OWN:
        if n.endswith('.py'):ast.parse((ROOT/n).read_text(encoding='utf8'))
    for n in ('scripts/runRestartBranchV236.py','scripts/auditRestartBranchV236.py'):
        r=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
        value=r.stdout.decode('utf8',errors='replace');assert r.returncode!=0 and 'closed' in value;closed[n]=value.strip()
    out.mkdir();evidence.mkdir()
    def copy(p,relative=None):
        target=evidence/(relative or p.relative_to(BASE));target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    for folder in BASE.glob('process_*'):
        for p in folder.rglob('*'):
            if p.is_file():copy(p)
    for name in ('contract.json','prepare_report.json'):copy(BASE/name)
    for p in (BASE/'build01').iterdir():
        if p.suffix in ('.json','.d'):copy(p)
    for name in ('case.cfg','initial_check.json','native_completion.json','schema.json','history.csv','stages.csv','wall_stages.csv',
                 'steps.csv','admissibility.csv','rejections.csv','surface_trend.csv','benefit_decisions.csv','original_faces.txt'):
        copy(BASE/'rocket01'/name)
    for name in ('audit01','audit02','audit03','outcome01','workbench01','visual02'):
        for p in (BASE/name).iterdir():
            # 完整面元CSV、工作台数据及大场仅保留本地；紧凑报告、真实回退证据和图公开。
            if p.is_file() and not p.name.startswith('original_faces_') and p.name not in ('rocket_restart_comparison.json','case-data.js','index.html'):
                copy(p)
    local={p.relative_to(BASE).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*')
        if p.is_file() and (p.suffix in ('.bin','.exe') or p.name.startswith('wall_000') or p.name.startswith('original_faces_')
          or p.name in ('rocket_restart_comparison.json','case-data.js','index.html'))}
    save(evidence/'local_asset_index.json',dict(base='outputs/v236_validation',files=local,largeAssetsBackedUpToGitHub=False,
        fullBinaryFrames=[0,1,100,200,300],nativeFinalCSVGenerated=False,requiredRuntime='v235 observed thin archives; v220 identity15 / v167 exact objects; new execution budget required'))
    for n in OWN:copy(ROOT/n,Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    elapsed=time.monotonic()-start;assert totals['offline']+elapsed<120 and totals['build']<60 and totals['audit']<240 and totals['cfd']<2400
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<1024**3
    report=dict(passed=True,goalStatus='active',physicalQualification=False,newQualifiedAreaM2=0,newHeatAreaM2=0,
        realRocketUpdates=300,originalFaces=19752,originalAreaM2=1.9646292565715553,meshCells=705490,points=335726,
        initialSixStateCoordinatesBitEqual=True,nativeCFLReset=.1,exactHistoricalContinuation=False,
        fullFrames=[0,1,100,200,300],workbenchStages=ui['updates'],workbenchFields=34,
        nativeStageRows=3600,wallStageRows=600,rejectionEvents=466,uniqueRejectedNodes=19,
        stopReason='significant-deterioration',physicalMinus8Passed=False,nativeFinalRestartCSVWritten=False,nativeRestartPlanAcceptance=False,
        completeBinaryFinalFramePreserved=True,firstAuditFailuresAndFirstPlotPreserved=True,momentumSignSupplementRequired=True,
        protectedCount=len(protected),engineeringArtifactsUnchanged=len(eng['files']),candidateCount=len(candidates['files']),assets=assets,
        processSeconds=totals,finalizationSeconds=elapsed,maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),
        outputBytes=size,allRecordedInputsStable=True,allRecordedGroupsEmpty=True,historicalSourceBindings=historical,
        closedCLIs=closed,MATLABExecuted=False,sharedWebChanged=False,browserInteractionTested=False,
        revisedScientificImageVisuallyInspected=True,uniqueAssetsDeleted=False,newRocketMeshGenerated=False,
        unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v236/evidence_manifest.json']+['docs/validation/v236/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,outputBytes=size)))


if __name__=='__main__':main()
