"""模块职责：核验历史/并发资产、关闭执行入口、归档小证据；大场和原生程序只留E盘。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,sha,read,save
BASE=ROOT/'outputs/v223_validation';DOC=ROOT/'docs/validation/v223'
OWN=['scripts/inspectSU2WallBalanceV223.py','scripts/runSU2FrozenBalanceV223.py',
     'scripts/auditSU2FrozenBalanceV223.py','scripts/exportSU2PhaseComparisonV223.py',
     'scripts/finalizeSU2FrozenBalanceV223.py','native/su2PhaseAwareAir/frozenBalance.cpp',
     'tests/testSU2PhaseComparisonV223.js','docs/HYBRID_SU2_FROZEN_BALANCE_V223.md',
     'docs/validation/v223/PLAN.md','docs/validation/v223/NEXT_DECISION.md',
     'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False)
    evidence=DOC/'evidence';evidence.mkdir(exist_ok=False)
    # 模块一：只读核对独有大资产、旧对象及准确新包依赖；绝不删除或覆盖。
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json');assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    assets=read(ROOT/'outputs/v222_validation/final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    for rows in read(ROOT/'outputs/v220_validation/build02/members.json').values():
        for row in rows:assert sha(row['path'])==row['sha256']
    for n,h in read(BASE/'inspection01/report.json')['sourceFrameSHA256'].items():assert sha(n)==h,n
    # 当前共享前端只记录摘要，不复制其未提交源码到公开证据。
    shared={n:sha(ROOT/n) for n in ('web/viewer.js','web/section.js','tests/verifyWebSectionExport.py')}
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
    closed=subprocess.run([sys.executable,str(ROOT/'scripts/runSU2FrozenBalanceV223.py'),'run2'],
                          stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert closed.returncode!=0 and b'closed' in closed.stdout
    assert read(BASE/'audit02/report.json')['passed'] and read(BASE/'workbench01/report.json')['passed']
    # 模块二：保留所有成功/失败回执及准确执行源，体场field.csv/可执行程序不公开。
    records=[]
    for p in BASE.glob('process_*/report.json'):
        r=read(p);assert r['inputsStable'] and not r['liveGroupAfterCleanup'];records.append(r)
        shutil.copytree(p.parent,evidence/p.parent.name)
    for folder in ('inspection01','localize01','audit02','workbench01'):
        shutil.copytree(BASE/folder,evidence/folder)
    for folder in ('frozen01','frozen02'):
        dest=evidence/folder;dest.mkdir()
        for n in ('case.cfg','stages.csv','native.json','wall.csv'):
            p=BASE/folder/n
            if p.is_file():shutil.copyfile(p,dest/n)
    dest=evidence/'build01';dest.mkdir()
    for n in ('inputs.json','command.json','prepare_report.json'):
        shutil.copyfile(BASE/'build01'/n,dest/n)
    for n in OWN:
        rel=Path(n) if n!='.gitattributes' else Path('.gitattributes.snapshot')
        dest=evidence/'final_source'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/n,dest)
    for n in ('COPYING','AUTHORS.md'):
        shutil.copyfile(ROOT/'docs/validation/v220/evidence'/n,evidence/n)
    # 模块三：资源分类不把父子重复计为串行成本；读取/开发与监督运行分别说明。
    totals={kind:sum(r['elapsedSeconds'] for r in records if r['budgetKind']==kind)
            for kind in ('offline','build','query','parent')}
    assert totals['offline']<120 and totals['build']<60 and totals['query']<120
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert size<256*1024**2
    assert all(sha(ROOT/n)==h for n,h in shared.items())
    assert all(sha(ROOT/n)==h for n,h in protected.items())
    seconds=time.monotonic()-start
    assert totals['offline']+seconds<120
    report=dict(passed=True,goalStatus='active',physicalQualification=False,newQualifiedAreaM2=0,
        protectedCount=29,oldObjectsStable=len(old),dependenciesStable=len(deps),assets=assets,
        sharedSHA256=shared,sharedSourceCopied=False,candidateEntryClosed=True,
        closedEntryOutput=closed.stdout.decode('utf8',errors='replace'),flowUpdates=0,spaceEvaluations=2,
        processBudgetSecondsBeforeFinal=totals,finalizationSeconds=seconds,
        offlineIncludingFinalSeconds=totals['offline']+seconds,outputBytesBeforeFinalReport=size,
        outputUnder256MiB=True,initialInspectionOnlineOutputScopeMatched=False,
        nativePeakRSSBytes=read(BASE/'process_query2/report.json')['peakObservedGroupRSSBytes'],
        nextDecision='H1/H2 same-model engineering boundary-layer entry; no direct restart from 199',
        independentFullPhysicalConservationPassed=False,browserInteractionTested=False,
        failedAttemptsRetained=['process_query','process_parent','process_audit'])
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    save(out/'delivery_files.json',OWN+['docs/validation/v223/evidence_manifest.json']+
        [str(Path('docs/validation/v223')/n).replace('\\','/') for n in sorted(manifest)])
    print(json.dumps(dict(passed=True,files=len(manifest),protectedCount=29,
        offlineIncludingFinalSeconds=totals['offline']+seconds,physicalQualification=False)))

if __name__=='__main__':main()
