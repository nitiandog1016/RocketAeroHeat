"""模块职责：归档新工程层的完整MATLAB文件及小证据，保护共享/原生资产；不制作工程包。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from runPhaseAwareLayerV224 import SOURCES
BASE=ROOT/'outputs/v224_validation';DOC=ROOT/'docs/validation/v224'
OWN=['phaseAwareAirState.m','phaseAwareAirReferenceState.m','phaseAwareMeshIdentity.m',
    'calcPhaseAwareLaminarHeating.m','calcSurfaceHeating.m','buildPhaseAwareLaminarViewerData.m',
    'buildAerothermalViewerData.m','runAeroHeatSingleCase.m','tests/testPhaseAwareLaminarHeating.m',
    'tests/testPhaseAwareLayerViewerV224.js','scripts/runPhaseAwareLayerV224.py',
    'scripts/auditPhaseAwareLayerV224.py','scripts/finalizePhaseAwareLayerV224.py',
    'docs/HYBRID_PHASE_LAMINAR_CHAIN_V224.md','docs/validation/v224/PLAN.md',
    'docs/validation/v224/NEXT_DECISION.md','AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    evidence=DOC/'evidence';evidence.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json');assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    # 只用已封存清单核对并发任务，不复制或发布该任务的未提交代码。
    other=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in other['files']:assert sha(ROOT/row['path'])==row['sha256'],row['path']
    candidate=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidate['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    previous=read(ROOT/'outputs/v223_validation/final01/report.json');assets=previous['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    for name in OWN:
        p=ROOT/name
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
    for name in SOURCES:assert sha(ROOT/name)==sha(BASE/'matlab01/executed_source'/name),name
    assert read(BASE/'audit01/report.json')['passed'] and read(BASE/'matlab01/report.json')['passed']
    closed=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseAwareLayerV224.py'),'matlab01'],
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert closed.returncode!=0 and b'closed' in closed.stdout
    records=[]
    for p in BASE.glob('process_*/report.json'):
        r=read(p);assert r['passed'] and r['inputsStable'] and not r['liveGroupAfterCleanup'];records.append(r)
        shutil.copytree(p.parent,evidence/p.parent.name)
    for folder in ('matlab01','audit01'):shutil.copytree(BASE/folder,evidence/folder)
    # MATLAB原始日志混合ASCII及本机中文编码；保留原字节，UTF8版仅为带SHA的派生阅读件。
    raw=BASE/'process_matlab01/stdout.txt';decoded=raw.read_bytes().decode('gb18030',errors='replace')
    (evidence/'process_matlab01/stdout_utf8.txt').write_text(decoded,encoding='utf8')
    save(evidence/'process_matlab01/log_decoding.json',dict(sourceSHA256=sha(raw),
        sourceEncodingAssumption='gb18030',replacementCharacters=decoded.count('\ufffd'),
        derivedFile='stdout_utf8.txt',originalBytesRetained=True,numericalRerun=False))
    for name in OWN:
        rel=Path(name) if name!='.gitattributes' else Path('.gitattributes.snapshot')
        dest=evidence/'final_source'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,dest)
    totals={k:sum(r['elapsedSeconds'] for r in records if r['budgetKind']==k) for k in ('offline','matlab')}
    assert totals['offline']<120 and totals['matlab']<180
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<128*1024**2
    elapsed=time.monotonic()-start;assert totals['offline']+elapsed<120
    result=dict(passed=True,goalStatus='active',physicalQualification=False,newQualifiedAreaM2=0,
        protectedCount=29,sharedEngineeringArtifactsUnchanged=len(other['files']),sharedCandidatesUnchanged=True,
        oldObjectsStable=len(old),dependenciesStable=len(deps),assets=assets,
        actualMatlabSourceFiles=len(SOURCES),completeReplacementMatlabFiles=[n for n in OWN if n.endswith('.m')],
        processBudgetSeconds=totals,finalizationSeconds=elapsed,outputBytesBeforeFinalReport=size,
        validationEntryClosed=True,closedEntryOutput=closed.stdout.decode('utf8',errors='replace'),
        CFDRuns=0,nativeQueries=0,engineBuilds=0,sharedWebSourceCopied=False,browserInteractionTested=False,
        testedControlFaces=81,realRocketNewCandidateAreaM2=0,internalPhaseEdgeProviderConnected=False,
        localCFDProfileReady=False,wholeSurfaceGoalComplete=False)
    save(out/'report.json',result);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v224/evidence_manifest.json']+
        ['docs/validation/v224/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,files=len(manifest),protected=29,
        sharedEngineeringArtifactsUnchanged=len(other['files']),seconds=elapsed)))

if __name__=='__main__':main()
