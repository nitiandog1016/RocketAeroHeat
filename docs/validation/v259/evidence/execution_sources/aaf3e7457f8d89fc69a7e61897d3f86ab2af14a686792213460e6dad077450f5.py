"""模块职责：只读固定git对象，绑定实际块推进源码/算例和官方对照；不启动求解器。"""
import sys,subprocess,hashlib,time,shutil
from pathlib import Path
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v259_validation';DOC=ROOT/'docs/validation/v259';CLOSED=False
COMMIT='f53f4609a0331d48efee69a4e4f3c3598378cc03';REPO='/opt/rocket-aeroheat/eilmer-v5.0.0/source'
FILES=['COPYING','src/lmr/commands/runsim.d','src/lmr/globalconfig.d','doc/lmr/config-options-user.adoc',
    'src/eilmer/main.d','src/eilmer/simcore.d','src/eilmer/block.d','src/eilmer/fluidblock.d',
    'src/eilmer/sfluidblock.d','src/eilmer/ufluidblock.d','src/eilmer/bc/boundary_condition.d',
    'src/eilmer/bc/ghost_cell_effect/extrapolate_copy.d','src/eilmer/fbarray.lua','src/eilmer/fluidblockarray.d',
    'examples/eilmer/3D/flat-plate-laminar/lam_flat_plate_march.lua','examples/eilmer/3D/flat-plate-laminar/ideal-air.inp',
    'examples/eilmer/3D/flat-plate-laminar/readme.txt','examples/eilmer/3D/flat-plate-laminar/post.sh',
    'examples/eilmer/3D/flat-plate-laminar/boundary-layer-profile-clbl.data','examples/eilmer/2D/flat-plate-laminar/usg/lam_flat_plate_march.lua']


def collect():
    out=BASE/'source01';out.mkdir(exist_ok=False);tree={r['path']:r for r in read(BASE/'official_tree.json')['tree']};rows=[]
    assert len(FILES)<=40
    for index,name in enumerate(FILES):
        assert name in tree and tree[name]['type']=='blob'
        argv=['C:/Windows/System32/wsl.exe','-d','RocketAeroHeat-Ubuntu-24.04','-u','rah','--exec',
            '/usr/bin/timeout','15','/usr/bin/taskset','-c','0','/usr/bin/prlimit','--as=536870912','--cpu=10','--',
            '/usr/bin/git','-C',REPO,'show',COMMIT+':'+name]
        started=time.monotonic();r=subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=25)
        (out/(f'stderr_{index:02d}.raw')).write_bytes(r.stderr)
        save(out/(f'command_{index:02d}.json'),dict(argv=argv,exitCode=r.returncode,elapsedSeconds=time.monotonic()-started,
            LinuxGitCPU=1,LinuxGitAddressSpaceLimitBytes=536870912,WindowsSupervisorRSSDoesNotMeasureLinuxRSS=True,
            stderrSHA256=hashlib.sha256(r.stderr).hexdigest()))
        assert r.returncode==0,(name,r.stderr);assert len(r.stdout)==tree[name]['size'] and len(r.stdout)<2*1024**2
        blob=hashlib.sha1(b'blob '+str(len(r.stdout)).encode()+b'\0'+r.stdout).hexdigest();assert blob==tree[name]['sha'],name
        p=out/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(r.stdout)
        rows.append(dict(path=name,gitBlobSHA1=blob,sha256=sha(p),bytes=len(r.stdout)))
    refs=read(BASE/'connector_refs.json')['refs']+read(BASE/'connector_extra_refs.json');comparisons=[]
    for item in refs:
        p=out/item['path'];lines=p.read_text(encoding='utf8').splitlines();excerpt='\n'.join(lines[item['first']-1:item['last']])
        assert excerpt.strip('\n')==item['content'].strip('\n'),item['path'];assert tree[item['path']]['sha']==item['sha']
        comparisons.append(dict(path=item['path'],first=item['first'],last=item['last'],connectorBlobMatches=True,selectedLinesMatch=True))
    save(out/'report.json',dict(passed=True,commit=COMMIT,files=rows,connectorComparisons=comparisons,
        actualEngineOrCFDExecuted=False,sourceModified=False,systemModified=False,methodAdopted=False,physicalQualification=False))


def main():
    if CLOSED:raise SystemExit('v259 fixed marching source review closed; follow NEXT_DECISION')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',BASE/'official_tree.json',BASE/'connector_refs.json',BASE/'connector_extra_refs.json',
        ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py']
    snap=BASE/'source_collect';snap.mkdir(exist_ok=False)
    for p in locks:
        if p==Path(__file__).resolve() or p==DOC/'PLAN.md':shutil.copyfile(p,snap/(sha(p)+p.suffix))
    return supervise(guard,BASE,'collect',[sys.executable,str(Path(__file__).resolve()),'task_collect'],110,32,locks,ROOT,'offline')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v259 fixed marching source review closed; follow NEXT_DECISION')
    if len(sys.argv)==2 and sys.argv[1]=='task_collect':collect()
    elif len(sys.argv)==1:raise SystemExit(0 if main() else 1)
    else:raise SystemExit('stage not admitted')
