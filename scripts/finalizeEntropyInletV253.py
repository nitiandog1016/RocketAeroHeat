"""模块职责：封存真实入口/失败/独立审核和科研出口；本地交付，不推送。"""
from pathlib import Path
import ast,subprocess,sys,shutil,time,json
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v253_validation';DOC=ROOT/'docs/validation/v253'
CLIS=['scripts/'+n+'.py' for n in ('runEntropyInletV253','produceEntropyInletV253','auditEntropyInletV253','exportEntropyInletV253')]
OWN=CLIS+['scripts/phaseAwareEntropyInlet.py','scripts/finalizeEntropyInletV253.py',
    'web/entropy-inlet-v253.html','web/entropy-inlet-v253.js','tests/testEntropyInletViewerV253.js',
    'docs/ENTROPY_INLET_V253.md','docs/validation/v253/PLAN.md','docs/validation/v253/NEXT_DECISION.md',
    'docs/CURRENT_STATE.md','docs/EXECUTION_INDEX.md','.gitattributes']


def main():
    started=time.monotonic();out=BASE/'final01';evidence=DOC/'evidence';assert not out.exists() and not evidence.exists()
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    engineering=read(ROOT/'harness_tasks/engineering_bl_20261108/delivery_manifest.json')
    for row in engineering['files']:assert sha(ROOT/row['path'])==row['sha256']
    candidates=read(ROOT/'harness_tasks/engineering_bl_20261108/candidate_manifest.json')
    for row in candidates['files']:
        assert sha(ROOT/row['file'])==row['baselineSHA256']
        assert sha(ROOT/'harness_tasks/engineering_bl_20261108/candidates'/row['file'])==row['candidateSHA256']
    previous=read(ROOT/'docs/validation/v252/evidence/final_report.json')
    for n,h in previous['unrelatedTrackedSHA256'].items():assert sha(ROOT/n)==h,n
    records={p.parent.name:read(p) for p in BASE.glob('process_*/report.json')};assert len(records)==6
    snapshots={sha(p):p for p in BASE.glob('process_*/execution_sources/*') if p.is_file()};bindings=[];cache={};totals={};failures={}
    for name,r in records.items():
        assert r['status']=='terminal' and r['inputsStable'] and not r['liveGroupAfterCleanup'],name
        assert read(BASE/name/'terminal_capacity.json')['passed']
        totals[r['budgetKind']]=totals.get(r['budgetKind'],0)+r['elapsedSeconds']
        if not r['passed']:failures[name]=r['exitCode']
        for n,h in r['inputSHA256'].items():
            if n not in cache:cache[n]=sha(n)
            if cache[n]!=h:
                assert h in snapshots,(n,h)
                bindings.append(dict(process=name,path=n,executedSHA256=h,currentSHA256=cache[n],snapshot=snapshots[h].relative_to(BASE).as_posix()))
    assert set(failures)=={'process_produce','process_produce2','process_audit'}
    assert totals['offline']<=120 and totals['audit']<=180
    production=read(BASE/'inlet02/report.json');audit=read(BASE/'audit02/report.json');ui=read(BASE/'workbench01/report.json')
    assert production['passed'] and audit['passed'] and ui['passed']
    assert not audit['onePassJointInletAccepted'] and production['commonRingEdges']==46
    assert all(v['availableInlets']==41 and v['unavailableInlets']==5 and v['correctionsBelowOnePercent']==0 for v in production['branches'].values())
    assert ui['fullNormalSamples']==audit['physicalConservativeSamples']==73923
    closed={}
    for n in CLIS:
        p=subprocess.run([sys.executable,str(ROOT/n),'not-admitted'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15)
        text=p.stdout.decode('utf8',errors='replace');assert p.returncode!=0 and 'closed' in text;closed[n]=text.strip()
    for n in OWN:
        if n.endswith('.py'):ast.parse((ROOT/n).read_text(encoding='utf8'))
    # 方法证据定向复用一手报告原文本，不下载/复制整套历史。
    textpath=ROOT/'outputs/v52_validation/research/cr112180_text.txt';text=textpath.read_text(encoding='utf8')
    method=BASE/'method01';method.mkdir(exist_ok=False)
    excerpt=text[text.index('PDF PAGE 29'):text.index('PDF PAGE 33')]
    (method/'cr112180_mass_excerpt.txt').write_text(excerpt,encoding='utf8')
    save(method/'source.json',dict(url='https://ntrs.nasa.gov/api/citations/19730005231/downloads/19730005231.pdf',
        localTextPath=str(textpath),textSHA256=sha(textpath),equations=[65,66,67,68,69],stagnationSource='Appendix B16, B=1',
        decision='adopt explicit metric/label mass matching as engineering inlet candidate; reject unchanged-profile one-pass coupling',
        exactThreeDimensionalShockShape=False,smallCrossflowAssumption=True,externalCodeRedistributed=False))
    out.mkdir();evidence.mkdir();local={}
    def copy(p,rel):
        dest=evidence/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in sorted(BASE.rglob('*')):
        if not p.is_file() or out in p.parents or '__pycache__' in p.parts:continue
        rel=p.relative_to(BASE);top=rel.parts[0]
        keep=top.startswith('process_') or top in ('method01','inlet01','audit01','audit02')
        keep=keep or (top=='inlet02' and p.suffix=='.json')
        keep=keep or (top=='workbench01' and (p.suffix=='.png' or p.name in ('report.json','frontend_test.txt')))
        if keep:
            assert p.stat().st_size<2*1024**2;copy(p,rel)
        else:local[rel.as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
    for n in OWN:copy(ROOT/n,Path('final_source')/(n if n!='.gitattributes' else '.gitattributes.snapshot'))
    save(evidence/'local_asset_index.json',local)
    report=dict(passed=True,goalStatus='active',wholeSurfaceGoalComplete=False,baseCommit='65523caaedb00a0cceb466cab1a3e333d23599ce',
        originalFaces=19752,originalAreaM2=1.9646292565715553,commonRingEdges=46,availableProfiles=123,missingProfiles=15,
        physicalConservativeSamples=73923,onePassJointInletAccepted=False,profileReadyForCFD=False,
        newQualifiedAreaM2=0,newHeatAreaM2=0,physicalQualification=False,independentAudit=audit,
        productionSummary=production,workbenchSummary=ui,scientificImagesVisuallyInspected=True,browserInteractionTested=False,
        protectedCount=len(protected),engineeringFilesUnchanged=len(engineering['files']),candidateCount=len(candidates['files']),
        unrelatedTrackedSHA256=previous['unrelatedTrackedSHA256'],historicalSourceBindings=bindings,
        allReceiptsTerminalStableGroupsEmpty=True,executionFailuresPreserved=failures,processSeconds=totals,
        maxSampledGroupRSSBytes=max(r['peakObservedGroupRSSBytes'] for r in records.values()),closedCLIs=closed,ownedFiles=OWN,
        localOnly=True,githubPushed=False,newCFD=False,MATLABExecuted=False,sharedWebGLChanged=False,
        originalSourceProfilesUnchanged=True,uniqueDataDeleted=False,newPackageInstalled=False,
        developmentReadsAndArchivingSeparate=True)
    report['outputLogicalBytesBeforeFinalReport']=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert report['outputLogicalBytesBeforeFinalReport']<128*1024**2
    report['finalizationSeconds']=time.monotonic()-started
    save(out/'report.json',report);save(evidence/'final_report.json',report)
    save(evidence/'manifest.json',{p.relative_to(evidence).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()})
    print(json.dumps({k:report[k] for k in ('passed','availableProfiles','missingProfiles','onePassJointInletAccepted','processSeconds','finalizationSeconds')}))


if __name__=='__main__':main()
