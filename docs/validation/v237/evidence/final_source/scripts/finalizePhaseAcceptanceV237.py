"""模块职责：核对保护资产、归档准确执行源和失败、生成小证据与本地大资产索引。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v237_validation';DOC=ROOT/'docs/validation/v237'
OWN=['native/su2PhaseAwareAir/PhaseStepAcceptanceV237.hpp','native/su2PhaseAwareAir/phaseStepControlV237.cpp',
 'scripts/runPhaseAcceptanceV237.py','scripts/auditPhaseAcceptanceV237.py','scripts/runPhaseBranchV237.py',
 'scripts/auditPhaseBranchV237.py','scripts/auditPhaseProductionV237.py','scripts/exportPhaseBranchV237.py',
 'scripts/renderPhaseBranchV237.py','scripts/reviewPhaseOutcomeV237.py','scripts/finalizePhaseAcceptanceV237.py','tests/testRestartBranchV237.js',
 'docs/PHASE_STEP_ACCEPTANCE_V237.md','docs/validation/v237/PLAN.md','docs/validation/v237/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']


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
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};historical=[];totals={}
    expectedFailures={'process_prepare','process_prepare2','process_dependencies','process_build','process_rocket','process_production_audit','process_export'}
    assert {n for n,r in records.items() if not r['passed']}==expectedFailures
    # 原失败不改写；按准确退出码登记，不把数值任务exit0等同科学资格。
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        for p,h in r['inputSHA256'].items():
            if sha(p)==h:continue
            snapshot=BASE/name/'execution_sources'/(h+Path(p).suffix)
            if not snapshot.exists():
                matches=[q for folder in ('source_history','stage01','stage02','build01') for q in (BASE/folder).glob('*')
                         if q.is_file() and q.suffix not in ('.exe','.obj','.a') and sha(q)==h]
                assert matches,(name,p,h);snapshot=matches[0]
            assert sha(snapshot)==h,(name,p)
            historical.append(dict(process=name,path=p,executedSHA256=h,snapshot=snapshot.relative_to(BASE).as_posix(),currentSHA256=sha(p)))
    control=read(BASE/'audit01/report.json');audit=read(BASE/'real_audit01/report.json');actual=read(BASE/'production_audit02/report.json')
    ui=read(BASE/'workbench01/report.json');plot=read(BASE/'visual01/report.json')
    assert all(x['passed'] for x in (control,audit,actual,ui,plot,read(BASE/'outcome01/report.json')))
    closed={}
    for n in OWN:
        if n.endswith('.py'):ast.parse((ROOT/n).read_text(encoding='utf8'))
    for n in ('scripts/runPhaseAcceptanceV237.py','scripts/auditPhaseAcceptanceV237.py','scripts/runPhaseBranchV237.py',
              'scripts/auditPhaseBranchV237.py','scripts/auditPhaseProductionV237.py'):
        r=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
        value=r.stdout.decode('utf8',errors='replace');assert r.returncode!=0 and 'closed' in value;closed[n]=value.strip()
    out.mkdir();evidence.mkdir()
    def copy(p,relative=None):
        target=evidence/(relative or p.relative_to(BASE));target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    for folder in BASE.glob('process_*'):
        for p in folder.rglob('*'):
            if p.is_file():copy(p)
    for name in ('contract.json','candidate_inputs.json','candidate_admission.json'):copy(BASE/name)
    for folder in ('stage01','stage02','build01','source_history'):
        for p in (BASE/folder).iterdir():
            if p.suffix in ('.json','.d','.patch') or folder=='source_history' or p.name=='restartBranchV237.cpp':copy(p)
    for folder in ('control01','audit01','real_audit01','production_audit01','production_audit02','workbench01','visual01','outcome01'):
        for p in (BASE/folder).iterdir():
            if p.is_file() and not p.name.startswith('original_faces_') and p.name not in ('rocket_restart_comparison.json','case-data.js','index.html'):copy(p)
    for p in (BASE/'rocket01').iterdir():
        if p.is_file() and p.suffix in ('.json','.txt','.csv','.cfg') and not p.name.startswith(('wall_000','restart_flow')):copy(p)
    local={p.relative_to(BASE).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BASE.rglob('*')
        if p.is_file() and (p.suffix in ('.bin','.exe','.obj','.a') or p.name.startswith(('wall_000','original_faces_','restart_flow'))
          or p.name in ('rocket_restart_comparison.json','case-data.js','index.html','CEulerSolver.cpp'))}
    save(evidence/'local_asset_index.json',dict(base='outputs/v237_validation',files=local,largeAssetsBackedUpToGitHub=False,
        requiredRuntime='v237 guarded CEuler + v235 observation archives + v220/v167 exact objects; new budget required'))
    for n in OWN:copy(ROOT/n,Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    elapsed=time.monotonic()-start
    assert totals['offline']+elapsed<120 and totals['build']<180 and totals['query']<60 and totals['audit']<240 and totals['cfd']<2405
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<1024**3
    report=dict(passed=True,goalStatus='active',physicalQualification=False,newQualifiedAreaM2=0,newHeatAreaM2=0,
        realRocketUpdates=audit['successfulUpdates'],originalFaces=19752,originalAreaM2=1.9646292565715553,meshCells=705490,points=335726,
        guardControlEvents=466,production=actual,realBranch=audit['nativeCompletion'],physicalMinus8Passed=audit['physicalMinus8All'],
        workbenchStages=ui['updates'],workbenchFields=ui['fields'],protectedCount=len(protected),engineeringArtifactsUnchanged=len(eng['files']),
        candidateCount=len(candidates['files']),assets=assets,processSeconds=totals,finalizationSeconds=elapsed,
        maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),outputBytes=size,
        allRecordedInputsStable=True,allRecordedGroupsEmpty=True,historicalSourceBindings=historical,closedCLIs=closed,
        failedProcesses={name:{k:r[k] for k in ('exitCode','stopReason','elapsedSeconds')} for name,r in records.items() if not r['passed']},
        MATLABExecuted=False,sharedWebChanged=False,browserInteractionTested=False,scientificImageVisuallyInspected=True,
        uniqueAssetsDeleted=False,newRocketMeshGenerated=False,unmeteredDevelopmentAndTargetedReadsDisclosed=True)
    save(out/'report.json',report);copy(out/'report.json',Path('final_report.json'))
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()};save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v237/evidence_manifest.json']+['docs/validation/v237/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),finalizationSeconds=elapsed,outputBytes=size,processSeconds=totals)))


if __name__=='__main__':main()
