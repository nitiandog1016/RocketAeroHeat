"""模块职责：封存唯一块推进候选的准确源码、明确拒绝与下一行动，不授新物理资格。"""
import json,sys,subprocess,shutil,time,hashlib
from pathlib import Path
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v259_validation';DOC=ROOT/'docs/validation/v259'
OWN=['scripts/reviewEilmerMarchV259.py','scripts/finalizeMarchReviewV259.py','docs/EILMER_MARCH_REVIEW_V259.md',
    'docs/validation/v259/PLAN.md','docs/validation/v259/NEXT_DECISION.md','docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md','.gitattributes']


def main():
    start=time.monotonic();ev=DOC/'evidence';out=BASE/'final01';assert not ev.exists() and not out.exists()
    baseline='fdc26459232763e01346d7524614db9522d9b652';assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==baseline
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256'];assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    previous=read(ROOT/'docs/validation/v258/evidence/final_report.json')
    for name,h in previous['unrelatedTrackedSHA256'].items():assert sha(ROOT/name)==h,name
    process=read(BASE/'process_collect/report.json');source=read(BASE/'source01/report.json')
    assert process['passed'] and process['status']=='terminal' and process['inputsStable'] and not process['liveGroupAfterCleanup']
    assert read(BASE/'process_collect/terminal_capacity.json')['passed'] and source['passed'] and len(source['files'])==20
    assert len(source['connectorComparisons'])==7 and process['elapsedSeconds']<110
    for row in source['files']:assert sha(BASE/'source01'/row['path'])==row['sha256']
    licenseBytes=(BASE/'license/gpl.txt').read_bytes();assert hashlib.sha1(b'blob '+str(len(licenseBytes)).encode()+b'\0'+licenseBytes).hexdigest()==read(BASE/'license/reference.json')['gitBlobSHA1']
    cli=ROOT/'scripts/reviewEilmerMarchV259.py';closed=subprocess.run([sys.executable,str(cli)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
    assert closed.returncode!=0 and b'closed' in closed.stdout
    facts=[
        dict(id='lmr-unimplemented',path='src/lmr/commands/runsim.d',first=373,last=375,finding='fixed lmr block-marching run branch prints NOT IMPLEMENTED and exits 1'),
        dict(id='legacy-forward-only',path='src/eilmer/simcore.d',first=855,last=926,finding='two active slices, frozen older upstream slice, temporary downstream extrapolation; no reverse feedback sweep'),
        dict(id='legacy-extrapolate',path='src/eilmer/bc/boundary_condition.d',first=258,last=272,finding='original convective and diffusive actions removed; zero-order GhostCellExtrapolateCopy installed'),
        dict(id='legacy-all-block-residency',path='src/eilmer/simcore.d',first=182,last=308,finding='single process initializes workspace, grid/flow arrays and reads solution for every local block before marching'),
        dict(id='unstructured-propagation-missing',path='src/eilmer/ufluidblock.d',first=743,last=746,finding='optional west-to-east initial-state propagation not implemented for unstructured grid'),
        dict(id='actual-example-scope',path='examples/eilmer/3D/flat-plate-laminar/lam_flat_plate_march.lua',first=1,last=90,finding='3D ideal-gas laminar flat plate at approximately Mach4 and given Tw300; not alpha20 rocket/SA/experiment')]
    decision=dict(candidate='fixed-lmr-block-marching-provider',status='rejected',legacyDirectSubstituteRejected=True,
        legacyHasNoFeedbackAtAll=False,legacyAllUnstructuredMarchingCategoricallyRejected=False,legacyTwoActiveSlicesMeansOnlyTwoResidentSlices=False,
        reasons=['fixed lmr execution unimplemented','legacy temporary extrapolation/frozen upstream cannot prove required incoming/downstream closure','no original-wall/near-wall mapping and native phase-air/turbulence compatibility proof'],
        evidence=facts,noNewEngineBuildOrRun=True,noNewWallHeatSolution=True,physicalQualification=False,newQualifiedAreaM2=0,
        memoryPeakForNewCandidateMeasured=False,memoryOOMClaimed=False,fullDomainResponsibilityPreserved=True,
        nextAction='stable piecewise caloric E/H inversion in existing same-model interface, with independent precision evidence before any flow integration',
        inverseIsProvedUniqueCauseOfV258Failure=False)
    save(BASE/'decision.json',decision)
    snapshots={sha(p):p for p in BASE.glob('source_collect/*') if p.is_file()};changes=[]
    for name,h in process['inputSHA256'].items():
        if sha(name)!=h:
            assert h in snapshots,(name,h);changes.append(dict(path=name,executedSHA256=h,currentSHA256=sha(name),snapshot=snapshots[h].relative_to(BASE).as_posix()))
    ev.mkdir();out.mkdir();archived={}
    def copy(p,rel):
        h=sha(p)
        if h in archived:return
        target=ev/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target);archived[h]=rel.as_posix()
    for p in (BASE/'source01').rglob('*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    for p in (BASE/'license').glob('*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    for p in (BASE/'process_collect').glob('*'):
        if p.is_file():copy(p,p.relative_to(BASE))
    for p in (BASE/'source_collect').glob('*'):
        if p.is_file():copy(p,Path('execution_sources')/p.name)
    for name in ['official_tree.json','connector_refs.json','connector_extra_refs.json','decision.json']:copy(BASE/name,Path(name))
    for change in changes:change['archivedSmallSource']=archived[change['executedSHA256']]
    assets={}
    for p in BASE.rglob('*'):
        if not p.is_file() or out in p.parents:continue
        h=sha(p);assets[p.relative_to(BASE).as_posix()]=dict(bytes=p.stat().st_size,sha256=h,archivedSmallEvidencePath=archived.get(h))
    save(ev/'local_asset_index.json',assets)
    report=dict(passed=True,deliveryIntegrityPassed=True,methodAdopted=False,decision=decision,goalStatus='active',wholeSurfaceGoalComplete=False,baseCommit=baseline,
        originalFaces=19752,originalAreaM2=1.9646292565715553,newQualifiedAreaM2=0,physicalQualification=False,profileReadyForCFD=False,
        collectedPinnedSourceFiles=20,officialConnectorComparisons=7,upstreamCommit=source['commit'],license='GPL-3.0-or-later',
        sourceCollectionSeconds=process['elapsedSeconds'],WindowsSampledGroupRSSBytes=process['peakObservedGroupRSSBytes'],
        LinuxRSSMeasured=False,LinuxGitPerCommandAddressSpaceLimitBytes=536870912,LinuxGitCPU=1,
        CFDExecutions=0,meshExecutions=0,MATLABExecutions=0,engineBuilds=0,propertyQueries=0,UIChanged=False,
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
        unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],historicalSourceBindings=changes,
        allSupervisedReceiptsTerminalStableGroupsEmpty=True,supervisedExecutionFailures={},
        discoveryFailures=['stale E-drive wsl.exe path missing, no process launched','JS syntax error before tool execution','PowerShell brace expansion parse error before commands'],
        WSLWarningsPreservedInReport=True,sourceDiscoveryShellSecondsApprox=[3.734,1.750],networkAndReadingTimeNotIncludedInCollection=True,
        localOnly=True,githubPushed=False,uniqueDataDeleted=False,fullFieldsArchivedInGit=False,nativeLibraryAndToolchainArchivedInGit=False,
        closedCLIs={'scripts/reviewEilmerMarchV259.py':closed.stdout.decode('utf8',errors='replace').strip()},ownedFiles=OWN,
        currentOwnedFileSHA256={name:sha(ROOT/name) for name in OWN if name not in ['docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md']})
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file());assert report['outputLogicalBytesBeforeFinalReport']<32*1024**2
    report['finalizationSeconds']=time.monotonic()-start;save(out/'report.json',report);save(ev/'final_report.json',report)
    save(ev/'manifest.json',{p.relative_to(ev).as_posix():sha(p) for p in ev.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ['deliveryIntegrityPassed','methodAdopted','collectedPinnedSourceFiles','sourceCollectionSeconds','finalizationSeconds']}))


if __name__=='__main__':main()
