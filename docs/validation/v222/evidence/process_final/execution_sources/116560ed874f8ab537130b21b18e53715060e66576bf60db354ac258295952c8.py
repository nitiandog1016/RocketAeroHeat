"""模块职责：保护历史资产，归档小证据并关闭候选；大场和执行文件留E盘。"""
from pathlib import Path
import ast,json,shutil,subprocess,sys,time
from runSU2RocketPhaseV222 import ROOT,BASE,CASE,DOC,ENGINE,MESH,sha,read,save

OWN=['scripts/runSU2RocketPhaseV222.py','scripts/su2PhaseAirFields.py','scripts/auditSU2RocketPhaseV222.py',
     'scripts/renderSU2RocketPhaseV222.py','scripts/finalizeSU2RocketPhaseV222.py','scripts/exportSU2PhaseResearchV222.py',
     'tests/testSU2PhaseResearchV222.js',
     'docs/HYBRID_SU2_ROCKET_PHASE_V222.md','docs/validation/v222/PLAN.md',
     'docs/validation/v222/AUDIT_CONTRACT.md','docs/validation/v222/NEXT_DECISION.md',
     'AGENTS.md','docs/PROJECT_GOALS.md','README.md','.gitattributes']

def main():
    start=time.monotonic();out=BASE/'final01';out.mkdir(exist_ok=False);evidence=DOC/'evidence';evidence.mkdir(exist_ok=False)
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    assert len(protected)==29
    for n,h in protected.items():assert sha(ROOT/n)==h,n
    assets=read(ROOT/'outputs/v221_validation/final01/report.json')['assets']
    for n,h in assets.items():assert sha(ROOT/n)==h,n
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(ROOT/'tools/runtime/su2-thermal-build/build-v117'/n)==h,n
    deps=read(ROOT/'outputs/v220_validation/build02/dependencies.json')
    for n,h in deps.items():assert sha(n)==h,n
    for rows in read(ROOT/'outputs/v220_validation/build02/members.json').values():
        for row in rows:assert sha(row['path'])==row['sha256']
    for n in OWN:
        p=ROOT/n
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf8'))
    result=subprocess.run([sys.executable,str(ROOT/'scripts/runSU2RocketPhaseV222.py'),'run'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=5)
    assert result.returncode!=0 and b'closed' in result.stdout
    records=[]
    for p in BASE.glob('process_*/report.json'):
        r=read(p);assert r['inputsStable'] and not r['liveGroupAfterCleanup'];records.append(r)
        shutil.copytree(p.parent,evidence/p.parent.name)
    audit=read(BASE/'audit01/report.json')
    # 物理失败可正常交付；不将关闭候选或接口通过冒充物理资格。
    for folder in ('source01','source02','geometry01','software01','audit01','visual01','visual02','workbench01'):
        if not (BASE/folder).exists():continue
        dest=evidence/folder
        if folder.startswith('source'):shutil.copytree(BASE/folder,dest)
        else:
            dest.mkdir()
            for p in (BASE/folder).iterdir():
                if p.suffix in ('.json','.csv','.png','.js'):shutil.copyfile(p,dest/p.name)
    for p in BASE.glob('entry_*.py'):shutil.copyfile(p,evidence/p.name)
    target=evidence/'rocket01';target.mkdir()
    for n in ('case.cfg','case_manifest.json','history.csv'):shutil.copyfile(CASE/n,target/n)
    for n in OWN:
        p=ROOT/n;rel=Path(n) if n!='.gitattributes' else Path('.gitattributes.snapshot')
        q=evidence/'final_source'/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
    for n in ('COPYING','AUTHORS.md'):shutil.copyfile(ROOT/'docs/validation/v220/evidence'/n,evidence/n)
    totals={k:sum(r['elapsedSeconds'] for r in records if r['budgetKind']==k) for k in ('offline','audit','cfd','parent')}
    assert totals['offline']<=120 and totals['audit']<=240
    actual_size=sum(p.stat().st_size for p in CASE.iterdir() if p.is_file())
    offline_size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file() and CASE not in p.parents)
    report=dict(passed=True,protectedCount=29,oldObjectsStable=len(old),dependenciesStable=len(deps),assets=assets,
        processBudgetSecondsBeforeFinal=totals,caseBytes=actual_size,offlineBytes=offline_size,caseUnder2048MiB=actual_size<=2048*1024**2,
        offlineUnder256MiB=offline_size<=256*1024**2,originalPreparationOnlineOutputScopeMatched=False,
        interfacePassed=audit['passed'],planAcceptance=audit['planAcceptance'],physicalQualification=False,newQualifiedAreaM2=0,
        goalStatus='active',candidateEntryClosed=True,closedEntryOutput=result.stdout.decode('utf8',errors='replace'),
        engineSHA256=sha(ENGINE),meshSHA256=sha(MESH),seconds=time.monotonic()-start)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',evidence/'final_report.json')
    manifest={p.relative_to(DOC).as_posix():sha(p) for p in evidence.rglob('*') if p.is_file()}
    save(DOC/'evidence_manifest.json',manifest)
    print(json.dumps(dict(passed=True,files=len(manifest),protectedCount=29,physicalQualification=False)))

if __name__=='__main__':main()
