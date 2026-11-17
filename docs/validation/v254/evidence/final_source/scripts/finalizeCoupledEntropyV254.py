"""模块职责：封存真实终态/失败/数值范围与原面出口，保护并发资产，仅本地交付。"""
import ast,json,subprocess,sys,shutil,time
from pathlib import Path
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v254_validation';DOC=ROOT/'docs/validation/v254'
CLIS=['scripts/'+n+'.py' for n in ['prepareCoupledEntropyV254','runCoupledEntropyV254','auditCoupledEntropyV254','auditCoupledSurfaceV254','exportCoupledEntropyV254']]
OWN=CLIS+['scripts/finalizeCoupledEntropyV254.py','solvePhaseAwareEntropyIntegral.m','tests/testCoupledEntropyV254.m',
    'tests/testCoupledEntropyViewerV254.js','web/coupled-entropy-v254.js','web/coupled-entropy-v254.html',
    'docs/COUPLED_ENTROPY_INTEGRAL_V254.md','docs/validation/v254/PLAN.md','docs/validation/v254/NEXT_DECISION.md',
    'docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md','.gitattributes']


def main():
    started=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()=='b38eab95f9790fe4e46bbfc29a62ca74fd0725df'
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,h in protected.items():assert sha(ROOT/name)==h,name
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for r in engineering['files']:assert sha(ROOT/r['path'])==r['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for r in candidates['files']:
        assert sha(ROOT/r['file'])==r['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/r['file'])==r['candidateSHA256']
    previous=read(ROOT/'docs/validation/v253/evidence/final_report.json')
    for name,h in previous['unrelatedTrackedSHA256'].items():assert sha(ROOT/name)==h,name
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==9
    snapshots={sha(p):p for p in BASE.glob('process_*/execution_sources/*') if p.is_file()}
    snapshots.update({sha(p):p for p in BASE.glob('matlab*/executed_source/**/*') if p.is_file()})
    recovered=read(BASE/'source_recovery01/report.json');assert recovered['passed']
    snapshots.update({sha(p):p for p in (BASE/'source_recovery01').glob('*') if p.is_file()})
    changes=[];cache={};totals={};failures={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert read(BASE/name/'terminal_capacity.json')['passed']
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        if not r['passed']:failures[name]=r['exitCode']
        for p,h in r['inputSHA256'].items():
            if p not in cache:cache[p]=sha(p)
            if cache[p]!=h:
                assert h in snapshots,(p,h)
                changes.append(dict(process=name,path=p,executedSHA256=h,currentSHA256=cache[p],snapshot=snapshots[h].relative_to(BASE).as_posix()))
    assert set(failures)=={'process_matlab01','process_export'}
    assert totals['offline']<=120 and totals['matlab']<=600 and totals['audit']<=180
    matlab=read(BASE/'matlab02/report.json');audit=read(BASE/'audit01/report.json');screen=read(BASE/'surface_audit01/report.json');ui=read(BASE/'workbench04/report.json')
    assert all(r['passed'] for r in [matlab,audit,screen,ui])
    assert not screen['fullCapNumericalGatePassed'] and screen['retainedFaces']==ui['candidateFaces']==1275
    assert ui['startupWithheldFaces']==21 and ui['outerBoundWithheldFaces']==10 and not ui['physicalQualification']
    # 两次最终图数值相同，前一次已逐图看过；校验后沿用图审，不宣称浏览器交互。
    for name in ['coupled_integral.png','rocket_coupled_integral.png']:
        assert sha(BASE/'workbench03'/name)==sha(BASE/'workbench04'/name)
    closed={}
    for name in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/name),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
        message=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in message;closed[name]=message.strip()
    for name in OWN:
        if name.endswith('.py'):ast.parse((ROOT/name).read_text(encoding='utf8'))
    out.mkdir();evidence.mkdir();assets={}
    def copy(p,rel):
        target=evidence/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,target)
    researchKeep={'program_upright-1.png','program_upright-4.png','cr111921-51.png','cr111921-53.png','cr112180-32.png','cr112180-33.png','cr111921_source.json','cr111922_sources.json'}
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=top.startswith('process_') or top in ['prepare01','matlab01','matlab02','audit01','surface_audit01','source_recovery01']
        keep=keep or (top.startswith('workbench') and (p.suffix=='.png' or p.name in ['report.json','frontend_test.txt']))
        keep=keep or (top=='research' and p.name in researchKeep)
        if keep:
            assert p.stat().st_size<2*1024**2;copy(p,rel)
        else:assets[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for name in OWN:copy(ROOT/name,Path('final_source')/(name if name!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',assets)
    report=dict(passed=True,deliveryStatus='tested engineering module with full-cap numerical rejection retained',
        goalStatus='active',wholeSurfaceGoalComplete=False,baseCommit='b38eab95f9790fe4e46bbfc29a62ca74fd0725df',
        originalFaces=19752,originalAreaM2=1.9646292565715553,actualRingPositions=46,actualProfilesFilled=0,
        physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,fullCapNumericalGatePassed=False,
        matlabSummary=matlab,independentAudit=audit,surfaceScreen=screen,workbenchSummary=ui,
        allReceiptsTerminalStableGroupsEmpty=True,executionFailuresPreserved=failures,processSeconds=totals,
        maximumSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),historicalSourceBindings=changes,
        frontendHistoricalSourceRecovery=recovered,
        scientificImagesVisuallyInspected=True,browserInteractionTested=False,protectedCount=len(protected),
        engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],
        closedCLIs=closed,ownedFiles=OWN,localOnly=True,githubPushed=False,newCFD=False,newPackageInstalled=False,
        uniqueDataDeleted=False,sharedWebGLChanged=False,originalV252ParentTerminalKnown=False,developmentReadsArchivingSeparate=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputLogicalBytesBeforeFinalReport']<128*1024**2
    report['finalizationSeconds']=time.monotonic()-started
    save(out/'report.json',report);save(evidence/'final_report.json',report)
    save(evidence/'manifest.json',{p.relative_to(evidence).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ['passed','fullCapNumericalGatePassed','actualRingPositions','actualProfilesFilled','processSeconds','finalizationSeconds']}))


if __name__=='__main__':main()
