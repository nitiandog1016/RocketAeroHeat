"""模块职责：封存波系已测部分与原失败，核对共享/原资产并生成准确Git交付清单。"""
from pathlib import Path
import ast,shutil,subprocess,sys,time,json
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v225_validation';DOC=ROOT/'docs/validation/v225'
OWN=['phaseAwareAirCaloric.m','solvePhaseAwareWave.m','marchPhaseAwareWavePath.m',
 'inspectPhaseAwareRocketInlet.m','tests/testPhaseAwareWaveV225.m','tests/completePhaseWaveGeometryV225.m',
 'scripts/runPhaseWaveV225.py','scripts/auditPhaseWaveV225.py','scripts/finalizePhaseWaveV225.py',
 'docs/HYBRID_PHASE_WAVES_V225.md','docs/validation/v225/PLAN.md','docs/validation/v225/NEXT_DECISION.md',
 'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);evidence=DOC/'evidence';evidence.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json');assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    other=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in other['files']:assert sha(ROOT/row['path'])==row['sha256'],row['path']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    assets=read(ROOT/'outputs/v224_validation/final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    numerical=['phaseAwareAirCaloric.m','solvePhaseAwareWave.m','marchPhaseAwareWavePath.m',
      'phaseAwareAirState.m','phaseAwareAirReferenceState.m','phaseAwareMeshIdentity.m','data/thermo/phase_aware_dry_air_v1.json']
    for n in numerical:assert sha(ROOT/n)==sha(BASE/'matlab02/executed_source'/n),n
    geometry=['inspectPhaseAwareRocketInlet.m','phaseAwareMeshIdentity.m','tests/completePhaseWaveGeometryV225.m']
    for n in geometry:assert sha(ROOT/n)==sha(BASE/'geometry01/executed_source'/n),n
    assert read(BASE/'audit01/report.json')['passed'] and read(BASE/'geometry01/report.json')['passed']
    observation=read(BASE/'audit01/matlab02_terminal_observation.json')
    assert not observation['originalSupervisorReportPresent'] and not observation['MATLABFullTestPassed']
    first=read(BASE/'process_matlab01/report.json');assert not first['passed'] and first['exitCode']==1
    for name in ('process_inspect01','process_matlab01','process_audit01','process_geometry01'):
        r=read(BASE/name/'report.json');assert r['inputsStable'] and not r['liveGroupAfterCleanup']
    closed=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseWaveV225.py'),'matlab01'],
      stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert closed.returncode!=0 and b'closed' in closed.stdout
    for p in BASE.iterdir():
        if p.is_dir() and p.name!='final01':shutil.copytree(p,evidence/p.name)
        elif p.is_file():shutil.copyfile(p,evidence/p.name)
    for name in ('process_matlab01','process_matlab02','process_geometry01'):
        raw=BASE/name/'stdout.txt';decoded=raw.read_bytes().decode('gb18030',errors='replace')
        (evidence/name/'stdout_utf8.txt').write_text(decoded,encoding='utf8')
        save(evidence/name/'log_decoding.json',dict(sourceSHA256=sha(raw),sourceEncodingAssumption='gb18030',
          replacementCharacters=decoded.count('\ufffd'),originalBytesRetained=True,numericalRerun=False))
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
        rel=Path(n) if n!='.gitattributes' else Path('.gitattributes.snapshot')
        dest=evidence/'final_source'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    geometryRun=read(BASE/'process_geometry01/report.json')
    charged=first['elapsedSeconds']+observation['budgetConservativelyChargedSeconds']+geometryRun['elapsedSeconds']
    assert charged<240
    offline=read(BASE/'process_inspect01/report.json')['elapsedSeconds']
    audit=read(BASE/'process_audit01/report.json')['elapsedSeconds'];elapsed=time.monotonic()-start
    assert offline<120 and audit<120 and elapsed<120
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert size<128*1024**2
    result=dict(passed=True,scope='delivery preservation and component evidence; not full MATLAB suite or physical acceptance',
      goalStatus='active',wholeSurfaceGoalComplete=False,physicalQualification=False,newQualifiedAreaM2=0,
      newCandidateHeatAreaM2=0,originalFaceStatesProduced=False,fullMatlabSuitePassed=False,
      waveSourceBoundToSavedResults=numerical,geometrySourceBoundToPassedRun=geometry,
      correctedFullSuiteFileNotRerun='tests/testPhaseAwareWaveV225.m',
      matlab02OriginalSupervisorCompletionMissing=True,matlabConservativelyChargedSeconds=charged,
      offlineSupervisedSeconds=offline,auditSupervisedSeconds=audit,finalizationSeconds=elapsed,
      protectedCount=29,sharedEngineeringArtifactsUnchanged=len(other['files']),sharedCandidatesUnchanged=True,
      oldObjectsStable=len(old),dependenciesStable=len(deps),assets=assets,
      outputBytesBeforeFinalReport=size,validationEntryClosed=True,
      CFDRuns=0,nativeQueries=0,engineBuilds=0,sharedWebSourceCopied=False,browserInteractionTested=False,
      completeReplacementMatlabFiles=[n for n in OWN if n.endswith('.m')])
    save(out/'report.json',result);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v225/evidence_manifest.json']+
      ['docs/validation/v225/'+n for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,files=len(manifest),protected=29,sharedEngineeringArtifactsUnchanged=len(other['files']),
      seconds=elapsed,matlabConservativelyChargedSeconds=charged)))

if __name__=='__main__':main()
