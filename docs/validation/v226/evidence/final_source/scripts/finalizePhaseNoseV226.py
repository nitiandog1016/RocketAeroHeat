"""模块职责：交付原面外缘/工作台和实际测试证据，保留失败并保护共享工程。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from runPhaseNoseV226 import SOURCES
BASE=ROOT/'outputs/v226_validation';DOC=ROOT/'docs/validation/v226'
OWN=['phaseAwarePayloadSHA.m','solvePhaseAwareNoseEdge.m','buildPhaseAwareNoseViewerData.m',
 'runPhaseAwareNoseCase.m','runAeroHeatSingleCase.m','buildAerothermalViewerData.m',
 'tests/testPhaseAwareNoseV226.m','tests/testPhaseNoseViewerV226.js','scripts/runPhaseNoseV226.py',
 'scripts/auditPhaseNoseV226.py','scripts/finalizePhaseNoseV226.py',
 'docs/HYBRID_PHASE_NOSE_EDGE_V226.md','docs/validation/v226/PLAN.md','docs/validation/v226/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']
def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);evidence=DOC/'evidence';evidence.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json');assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    other=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in other['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    assets=read(ROOT/'outputs/v225_validation/final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    for n in SOURCES:assert sha(ROOT/n)==sha(BASE/'matlab02/executed_source'/n),n
    for stage in ('matlab02','audit02','render01'):assert read(BASE/stage/'report.json')['passed']
    records=[]
    for p in sorted(BASE.glob('process_*/report.json')):
        r=read(p);assert r['inputsStable'] and not r['liveGroupAfterCleanup'];records.append(r)
    assert not read(BASE/'process_matlab01/report.json')['passed']
    assert not read(BASE/'process_audit01/report.json')['passed']
    assert read(BASE/'process_matlab02/report.json')['passed']
    assert read(BASE/'process_audit02/report.json')['passed']
    skipped=[]
    for folder in ('matlab01','matlab02','audit01','audit02','render01'):
        src=BASE/folder;target=evidence/folder;target.mkdir()
        for p in src.iterdir():
            if folder=='matlab01' and p.name in ('rocket_workbench.json','rocket_faces.csv'):
                skipped.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),retainedLocally=True));continue
            if p.is_dir():shutil.copytree(p,target/p.name)
            else:shutil.copyfile(p,target/p.name)
    # 失败首轮的大型重复字段保留E盘，公开其哈希与确切对照，失败元数据和执行源仍公开。
    assert sha(BASE/'matlab01/rocket_faces.csv')==sha(BASE/'matlab02/rocket_faces.csv')
    before=read(BASE/'matlab01/rocket_workbench.json');after=read(BASE/'matlab02/rocket_workbench.json')
    assert before['fields']==after['fields'] and before['vertices']==after['vertices'] and before['faces']==after['faces']
    save(evidence/'failed_duplicate_assets_retained_locally.json',dict(files=skipped,numericalFieldsIdentical=True,
      sourceIdentityDiffers=True,reason='failed lower name collision created invalid code-hash metadata; successful second run is authoritative'))
    for src in BASE.glob('process_*'):
        if src.is_dir():shutil.copytree(src,evidence/src.name)
    for name in ('process_matlab01','process_matlab02'):
        raw=BASE/name/'stdout.txt';decoded=raw.read_bytes().decode('gb18030',errors='replace')
        (evidence/name/'stdout_utf8.txt').write_text(decoded,encoding='utf8')
        save(evidence/name/'log_decoding.json',dict(sourceSHA256=sha(raw),sourceEncodingAssumption='gb18030',
          replacementCharacters=decoded.count('\ufffd'),originalBytesRetained=True,numericalRerun=False))
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        rel=Path(n) if n!='.gitattributes' else Path('.gitattributes.snapshot')
        q=evidence/'final_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    closed=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseNoseV226.py'),'matlab01'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert closed.returncode!=0 and b'closed' in closed.stdout
    total={kind:sum(r['elapsedSeconds'] for r in records if r['budgetKind']==kind) for kind in ('matlab','audit')}
    assert total['matlab']<240 and total['audit']<120
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<128*1024**2
    elapsed=time.monotonic()-start;assert elapsed<120
    candidate=read(BASE/'matlab02/report.json')['realRocket']
    result=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,physicalQualification=False,
      newQualifiedAreaM2=0,newHeatAreaM2=0,candidate=candidate,mainSingleCaseAndViewerUsed=True,
      protectedCount=29,sharedEngineeringArtifactsUnchanged=len(other['files']),sharedCandidatesUnchanged=True,
      oldObjectsStable=len(old),dependenciesStable=len(deps),assets=assets,executedMatlabSourceFiles=len(SOURCES),
      matlabFullRecoveryPassed=True,firstFailuresRetained=True,processBudgetSeconds=total,
      finalizationSeconds=elapsed,outputBytesBeforeFinalReport=size,validationEntryClosed=True,
      completeReplacementMatlabFiles=[n for n in OWN if n.endswith('.m')],
      CFDRuns=0,nativeQueries=0,engineBuilds=0,sharedWebSourceCopied=False,browserInteractionTested=False,
      finalScientificImageVisuallyInspected=True,originalGeometryUnchanged=True,largeDuplicateAssetsDeleted=False)
    save(out/'report.json',result);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v226/evidence_manifest.json']+
      ['docs/validation/v226/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,evidenceFiles=len(manifest),protected=29,sharedEngineeringArtifactsUnchanged=len(other['files']),
      seconds=elapsed,candidate=candidate,processBudgetSeconds=total)))
if __name__=='__main__':main()
