"""准确源码的有限资源路线回执；不运行求解、不重编、不推断实数全域一定能装下。"""
from pathlib import Path
import hashlib,json,time
start=time.monotonic();source=Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source')
out=Path('outputs/v180_validation/memory_route');assert not out.exists();out.mkdir()
selection={
 'src/lmr/main.d':[(108,126),(169,179),(280,288)],
 'src/lmr/commands/runsim.d':[(145,190)],
 'src/lmr/fluidblock.d':[(265,348)],
 'src/lmr/newtonkrylovsolver.d':[(174,189),(241,263),(1370,1408),(3180,3240),(3272,3318)],
 'src/lmr/makefile':[(54,60),(286,295),(415,422)]}
records=[]
for rel,ranges in selection.items():
    p=source/rel;raw=p.read_bytes();lines=raw.decode().splitlines()
    excerpt='\n\n'.join('\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1)) for a,b in ranges)+'\n'
    dest=out/(p.name+'.excerpt.txt');dest.write_text(excerpt)
    records.append(dict(path=rel,sourceSHA256=hashlib.sha256(raw).hexdigest(),ranges=ranges,excerpt=dest.name,
                        excerptSHA256=hashlib.sha256(dest.read_bytes()).hexdigest()))
cfg=json.loads(Path('outputs/v180_validation/real01/lmrsim/config').read_text())
run=json.loads(Path('outputs/v180_validation/real01_run/process00/report.json').read_text())
report=dict(passed=True,kind='source-route-review-only',sourceCommit='f53f4609a0331d48efee69a4e4f3c3598378cc03',sources=records,
    cells=cfg['block_0']['ncells'],faces=cfg['block_0']['nfaces'],
    observedGroupRSSBytes=run['peakObservedGroupRSSBytes'],addressSpaceLimitMiB=8192,
    hostPhysicalMemoryBytes=16844877824,wslMemoryGiB=8,wslSwapGiB=2,
    selectedNextRoute='Pinned official real-valued Newton-Krylov / direct lmr-run; finite-difference Jacobian with explicitly justified perturbation scales',
    realModeExists=True,realModeTestedHere=False,guaranteedWholeDomainMemoryFit=False,
    criticalCondition='Defaults for Frechet and preconditioner perturbations are 1e-30 for complex step; never carry them into real finite differences without checking actual operators.',
    wrapperExitObservation='delegateAndExecute returns C system status directly; actual solver child exit code was not separately captured in v180, so raw child code remains unknown.',
    nativeRuns=0,builds=0,elapsedSeconds=time.monotonic()-start)
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sources'}))
