"""只读核对准确分配/索引源码，生成未构建的存储适配草稿；不改变运行时。"""
from pathlib import Path
import difflib,hashlib,json,re,subprocess,time
ROOT=Path(__file__).resolve().parents[2];BASE=Path(__file__).resolve().parent
start=time.monotonic();out=BASE/'memory01';out.mkdir(exist_ok=False)
src=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source');lmr=src/'src/lmr'
assert subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip()=='f53f4609a0331d48efee69a4e4f3c3598378cc03'
assert not subprocess.check_output(['git','-C',str(src),'diff','--name-only'],text=True).strip()
spans={'fluidblock.d':[(237,302),(2016,2058),(2070,2140)],'fluidfvcell.d':[(190,220)],
       'ufluidblock.d':[(283,340)],'globalconfig.d':[(2418,2430)]}
files=[]
for name,ranges in spans.items():
    p=lmr/name;raw=p.read_bytes();lines=raw.decode().splitlines()
    excerpt='\n\n'.join('\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1)) for a,b in ranges)+'\n'
    (out/(name+'.txt')).write_text(excerpt)
    files.append(dict(path=str(p),sha256=hashlib.sha256(raw).hexdigest(),spans=ranges))
# 模块一：登记全部直接命名数组访问。别名/调用链仍须构建前有限覆盖，不冒称完整静态证明。
hits=[]
for p in sorted(lmr.rglob('*.d')):
    for i,line in enumerate(p.read_text().splitlines(),1):
        if re.search(r'\.(?:U[0-4]|dUdt[0-4])\b',line):hits.append(dict(path=str(p.relative_to(src)),line=i,text=line))
(out/'direct_accesses.json').write_text(json.dumps(hits,indent=2)+'\n')
# 模块二：保留同样四层和所有算子，仅草拟缩短显式数组未被单元视图使用的尾部。
p=lmr/'fluidblock.d';old=p.read_text()
allocation='''        // 显式每个时间层已有独立数组；只压缩尾部容量，不改变层数/索引/数值更新。
        // 非显式路径保留原分配，待独立运行包和边界索引检查后使用。
        size_t stage_storage = (ncells + nghost)*neq;
        if (myConfig.solverMode != SolverMode.transient) stage_storage *= nftl;
'''
new=old.replace('        celldata.U0.length = (ncells + nghost)*neq*nftl;',allocation+'        celldata.U0.length = stage_storage;',1)
new,count=re.subn(r'(celldata\.(?:U[1-4]|dUdt[0-4])\.length = )\(ncells \+ nghost\)\*neq\*nftl;',r'\1stage_storage;',new)
assert count==9 and new!=old
patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/src/lmr/fluidblock.d',tofile='b/src/lmr/fluidblock.d'))
target=ROOT/'patches/eilmer/explicit-dense-storage-v1.patch';assert not target.exists();target.parent.mkdir(parents=True,exist_ok=True);target.write_text(patch)
N=705490+19752+6572;neq=7;levels=4;buffers=2*levels
current=N*neq*levels*8*buffers;proposed=N*neq*8*buffers
r=dict(sourceCommit='f53f4609a0331d48efee69a4e4f3c3598378cc03',files=files,directAccessCount=len(hits),
       totalCellsIncludingGhosts=N,equations=neq,timeLevels=levels,stageBuffers=buffers,
       currentStageArrayPayloadBytes=current,proposedStageArrayPayloadBytes=proposed,potentialPayloadSavingBytes=current-proposed,
       allocationFormulaObserved='Each U0..U3 and dUdt0..dUdt3 allocates N*neq*nftl; cell views use id*neq..(id+1)*neq.',
       patch=str(target.relative_to(ROOT)),patchSHA256=hashlib.sha256(target.read_bytes()).hexdigest(),
       status='unbuilt-allocation-only-draft',nativeRuns=0,engineBuilds=0,
       caveat='Payload arithmetic only, not whole-process footprint or proof of the failed allocation. Full alias/bounds/native field equality validation pending. No assertion that saving is sufficient for 8GiB.',
       elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
