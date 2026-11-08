"""归档v211小证据与准确源摘录；原大网格、运行包及共享历史保持原位。"""
from pathlib import Path
import json,shutil,time,ast
from eilmerNKSAResourceBudget import ROOT,BASE,sha,save

DEST=ROOT/'docs/validation/v211/evidence'

def main():
    start=time.monotonic();assert json.loads((BASE/'final01/report.json').read_text())['passed']
    DEST.mkdir(exist_ok=False)
    # 模块一：保存布局执行时的原计划；后加扩域核算段不能追称事前已存在。
    plan=ROOT/'docs/validation/v211/PLAN.md';data=plan.read_bytes()
    lines=data.splitlines(keepends=True)
    before=b''.join(line for line in lines if not line.startswith('离线账同时读取'.encode('utf-8')))
    # 被移除段后多出的空行不属于原始计划。
    before=before.replace(b'\n\n\n',b'\n\n')
    snapshot=BASE/'layout01/plan_at_layout.md';snapshot.write_bytes(before)
    proc=json.loads((BASE/'layout_process01/report.json').read_text())
    plan_key='/mnt/e/RocketAeroHeat/docs/validation/v211/PLAN.md'
    assert sha(snapshot)==proc['inputSHA256'][plan_key]
    assert sha(BASE/'source01/executed_entry.py')==proc['inputSHA256']['/mnt/e/RocketAeroHeat/scripts/inspectEilmerNKSAResources.py']
    # 模块二：只摘录实际资源相关段落，绑定完整源SHA，不上传整份上游树。
    spans={
        'lmr/fluidblock.d':[(237,338),(1297,1346),(1348,1443),(1826,1863)],
        'lmr/ufluidblock.d':[(283,415),(418,462),(588,638)],
        'lmr/fluidfvcell.d':[(179,253),(1505,1596)],
        'lmr/flowgradients.d':[(31,107)],'lmr/lsqinterp.d':[(23,34),(252,316)],
        'lmr/fvinterface.d':[(112,141)],'lmr/fvvertex.d':[(26,57)],
        'lmr/coredata.d':[(20,52),(98,121)],'lmr/conservedquantities.d':[(154,205)],
        'lmr/turbulence.d':[(600,628)],
        'lmr/newtonkrylovsolver.d':[(670,687),(951,968),(1549,1556),(1592,1603),(2555,2562)],
        'nm/smla.d':[(60,84),(319,342)],'nm/bbla.d':[(24,37)],
        'geom/grid/usgrid.d':[(154,193),(218,270),(424,432)],'gas/gas_state.d':[(52,65)],
        'lmr/onedinterp.d':[(78,89)]}
    excerpts=[]
    for name,intervals in spans.items():
        p=BASE/'source01/source'/name;lines=p.read_text(encoding='utf-8').splitlines()
        excerpts.append(dict(path='src/'+name,sha256=sha(p),ranges=[dict(start=a,end=b,
            text='\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1))) for a,b in intervals]))
    p=BASE/'source02/jacobian.d';lines=p.read_text(encoding='utf-8').splitlines()
    excerpts.append(dict(path='src/lmr/jacobian.d',sha256=sha(p),ranges=[dict(start=24,end=91,
        text='\n'.join('%d: %s'%(i,lines[i-1]) for i in range(24,92)))]))
    save(BASE/'source_excerpts.json',excerpts)
    save(BASE/'manual_tool_notes.json',dict(scope='Tool-only failures, not native flow receipts',
        failures=['rg used nonexistent guessed script names in two read batches; continued with rg --files',
            'Optional Get-Content for nonexistent v210 GITHUB_DELIVERY.md; actual existing name DELIVERY.md',
            'PowerShell piped goal editor rejected by Python3.7 Non-UTF-8 SyntaxError before mutation',
            'UTF-8 file goal editor also rejected by Python3.7 tokenizer; same file executed successfully by Python3.12; exact tokenizer cause not investigated'],
        scientificFailure='resource01 report relative __file__ after four completed payload ledgers; retained false; final01 read-only recovery',
        CFDRuns=0,engineBuilds=0))
    rows=[];local=[]
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        if 'source01/source/' in p.relative_to(BASE).as_posix():
            local.append(dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size));continue
        assert p.stat().st_size<2*1024**2
        q=DEST/p.relative_to(BASE);q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(p)==sha(q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=p.relative_to(ROOT).as_posix(),sha256=sha(p),bytes=p.stat().st_size))
    names=['AGENTS.md','.gitattributes','docs/PROJECT_GOALS.md','docs/HYBRID_REGIONAL_NK_MEMORY_V211.md',
        'scripts/inspectEilmerNKSAResources.py','scripts/eilmerNKSAResourceBudget.py',
        'scripts/verifyEilmerNKSAResources.py','scripts/archiveEilmerNKSAResourcesV211.py']
    names += [p.relative_to(ROOT).as_posix() for p in (ROOT/'docs/validation/v211').glob('*.md')]
    for name in names:
        p=ROOT/name;q=DEST/'final_source'/('gitattributes.source' if name=='.gitattributes' else name)
        q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);assert sha(p)==sha(q)
        rows.append(dict(path=q.relative_to(ROOT).as_posix(),source=name,sha256=sha(p),bytes=p.stat().st_size))
        if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'))
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text())
    assert len(protected)==29 and all(sha(ROOT/p)==h for p,h in protected.items())
    size=sum(r['bytes'] for r in rows);assert size<6*1024**2
    result=dict(schema='RAH.V211Evidence.v1',files=rows,localSources=local,count=len(rows),bytes=size,
        sourceFiles=names,protectedFilesUnchanged=29,uniqueDeleted=0,CFDRuns=0,fullEngineBuilds=0,
        physicalQualification=False,goalComplete=False,seconds=time.monotonic()-start)
    save(ROOT/'docs/validation/v211/evidence_manifest.json',result)
    print(json.dumps({k:result[k] for k in ['count','bytes','protectedFilesUnchanged','seconds']}))

if __name__=='__main__':main()
