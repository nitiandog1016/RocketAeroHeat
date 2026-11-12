"""模块职责：隔离构建和一次整体Euler控制的资源监督；不提供CFD推进入口。"""
from pathlib import Path
import json,os,shutil,subprocess,sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from prepareSU2PhaseAirV220 import OLD_BUILD,OLD_SOURCE

BASE=ROOT/'outputs/v229_validation';DOC=ROOT/'docs/validation/v229'
SOURCE=BASE/'source';BUILD=BASE/'build01';PREV=ROOT/'outputs/v220_validation'
HELPER=ROOT/'native/su2PhaseAwareAir/eulerCompatibility.cpp'
CLOSED=True

def size():
    # 外部只读junction不是本轮输出，禁止沿它遍历旧第三方树。
    total=0
    for folder in BASE.iterdir():
        if folder.is_file():total+=folder.stat().st_size
        elif folder==SOURCE:
            for part in ('Common','SU2_CFD'):total+=sum(p.stat().st_size for p in (folder/part).rglob('*') if p.is_file())
        else:total+=sum(p.stat().st_size for p in folder.rglob('*') if p.is_file())
    return total

def prep():
    """模块一：复用准确头文件/旧对象，将四文件补丁置于独立源目录并提取实际依赖。"""
    out=BASE/'prepare01';out.mkdir(exist_ok=False);SOURCE.mkdir();BUILD.mkdir()
    deps=read(PREV/'build02/dependencies.json');tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    for n,h in deps.items():assert sha(n)==h,n
    staged=read(ROOT/'outputs/v228_validation/stage01/manifest.json');copied={}
    for part in ('Common/include','SU2_CFD/include'):
        for p in (PREV/'source'/part).rglob('*'):
            if not p.is_file():continue
            rel=p.relative_to(PREV/'source').as_posix();h=sha(p)
            expected=deps.get(str(p.resolve()),tree.get(str(OLD_SOURCE/rel)))
            assert expected==h,rel
            dest=SOURCE/rel;dest.parent.mkdir(parents=True,exist_ok=True)
            if rel in staged:
                source=Path(staged[rel]['isolatedPath']);assert sha(source)==staged[rel]['patchedSHA256']
            else:source=p
            shutil.copyfile(source,dest);copied[str(dest)]=sha(dest)
    for rel,row in staged.items():
        p=SOURCE/rel
        if p.exists():continue
        assert sha(row['isolatedPath'])==row['patchedSHA256']
        p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(row['isolatedPath'],p);copied[str(p)]=sha(p)
    subprocess.check_call(['cmd','/c','mklink','/J',str(SOURCE/'externals'),str(OLD_SOURCE/'externals')],stdout=subprocess.DEVNULL)
    commands={}
    for rel,cmd in read(PREV/'prepare03/commands.json').items():
        if rel not in staged:continue
        cmd=[x.replace(str(PREV/'source'),str(SOURCE)).replace(str(PREV/'build'),str(BUILD)) for x in cmd]
        obj=cmd[cmd.index('-o')+1];commands[rel]=cmd[:-2]+['-MD','-MF',obj+'.d']+cmd[-2:]
    helpercmd=read(ROOT/'outputs/v221_validation/build/command.json')
    helpercmd=[x.replace(str(PREV/'source'),str(SOURCE)).replace(str(PREV/'build02'),str(BUILD))
               .replace(str(ROOT/'outputs/v221_validation/build'),str(BUILD))
               .replace(str(ROOT/'native/su2PhaseAwareAir/compatibility.cpp'),str(HELPER)) for x in helpercmd]
    commands['helper']=helpercmd
    bindings={};all_dep_commands=[]
    for index,(name,cmd) in enumerate(commands.items()):
        filtered=[];i=0
        while i<len(cmd):
            if cmd[i] in ('-MF','-o'):i+=2;continue
            if cmd[i] in ('-MD','-c','-Wl,--start-group','-Wl,--end-group') or cmd[i].endswith('.a'):i+=1;continue
            filtered.append(cmd[i]);i+=1
        dep=out/('unit%d.d'%index);actual=filtered+['-M','-MT','dep','-MF',str(dep)]
        subprocess.check_call(actual,cwd=str(OLD_BUILD),timeout=20);all_dep_commands.append(actual)
        for n in dep.read_text(encoding='utf8').replace('\\\n',' ')[4:].split():
            p=Path(n);p=(p if p.is_absolute() else OLD_BUILD/p).resolve();h=sha(p)
            if str(p) in copied:assert copied[str(p)]==h
            elif p==HELPER:pass
            elif str(p) in deps:assert deps[str(p)]==h
            else:
                try:rel=p.relative_to(SOURCE).as_posix()
                except ValueError:rel=None
                if rel is not None:assert sha(OLD_SOURCE/rel)==h,rel
            bindings[str(p)]=h
    members=read(PREV/'build02/members.json')
    for rows in members.values():
        for row in rows:assert sha(row['path'])==row['sha256']
    old=read(ROOT/'outputs/v167_validation/build01/objects_after.json')
    for n,h in old.items():assert sha(OLD_BUILD/n)==h,n
    shutil.copyfile(HELPER,out/'executed_helper.cpp')
    save(out/'source.json',copied);save(out/'dependencies.json',bindings);save(out/'commands.json',commands)
    save(out/'dependency_commands.json',all_dep_commands);save(out/'old_members.json',members)
    control=BASE/'control01';control.mkdir()
    shutil.copyfile(ROOT/'outputs/v166_validation/small01/mesh.su2',control/'mesh.su2')
    cfg={}
    for line in (ROOT/'outputs/v228_validation/prepare03/case.proposed.cfg').read_text(encoding='utf8').splitlines():
        if '=' in line:k,v=line.split('=',1);cfg[k.strip()]=v.strip()
    cfg.update(MESH_FILENAME='mesh.su2',KIND_TURB_MODEL='NONE',ITER='0',
      MARKER_RIEMANN='(farfield, STATIC_SUPERSONIC_INFLOW_PT, 1197, 350, 2.001, 0.401, 0.299)',
      OUTPUT_FILES='(RESTART_ASCII)',WRT_RESTART_OVERWRITE='YES',OUTPUT_WRT_FREQ='1')
    (control/'case.cfg').write_text('\n'.join(k+'= '+v for k,v in cfg.items())+'\n',encoding='utf8')
    save(out/'report.json',dict(passed=True,sourceFiles=len(copied),dependencies=len(bindings),oldObjects=len(old),
      translationUnits=3,helper=1,CFDRuns=0,actualEulerCompatibility=False))

def build():
    """模块二：只替换三对象并链接独立薄归档；旧归档/对象内容均只读。"""
    commands=read(BASE/'prepare01/commands.json');receipts=[]
    import time
    for name,cmd in commands.items():
        if name=='helper':continue
        t=time.monotonic();subprocess.check_call(cmd,cwd=str(OLD_BUILD))
        receipts.append(dict(name=name,seconds=time.monotonic()-t));save(BUILD/'compile_receipts.json',receipts)
    ar=str(Path(commands['helper'][0]).parent/'ar.exe');members=read(BASE/'prepare01/old_members.json');new={}
    names={'CConfig.obj','CEulerSolver.obj','roe.obj'}
    for group,rows in members.items():
        result=[]
        for row in rows:
            p=Path(row['path']);assert sha(p)==row['sha256']
            if p.name in names:p=BUILD/p.name
            result.append(dict(path=str(p),sha256=sha(p)))
        subprocess.check_call([ar,'crsT',str(BUILD/group)]+[r['path'] for r in result],cwd=str(OLD_BUILD))
        new[group]=result
    save(BUILD/'members.json',new)
    link=[x.replace(str(PREV/'build02'),str(BUILD)) for x in read(PREV/'build02/link_command.json')]
    save(BUILD/'link_command.json',link);subprocess.check_call(link,cwd=str(OLD_BUILD))
    subprocess.check_call(commands['helper'],cwd=str(OLD_BUILD))
    save(BUILD/'report.json',dict(passed=True,engineSHA256=sha(BUILD/'SU2_CFD.exe'),helperSHA256=sha(BUILD/'compatibility.exe'),
      identity='rah-su2-phase-air-euler-eplus300k-v1',compiled=True,nativeCompatibility=False,CFDRuns=0))

def locks():
    result=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',HELPER,
      ROOT/'outputs/v228_validation/stage01/manifest.json',ROOT/'outputs/v228_validation/prepare03/case.proposed.cfg',
      PREV/'build02/members.json',PREV/'build02/dependencies.json',ROOT/'outputs/v167_validation/build01/objects_after.json']
    for manifest in (PREV/'build02/dependencies.json',BASE/'prepare01/dependencies.json'):
        if manifest.exists():result += [Path(n) for n in read(manifest)]
    for rows in read(PREV/'build02/members.json').values():result += [Path(row['path']) for row in rows]
    for row in read(ROOT/'outputs/v228_validation/stage01/manifest.json').values():result += [Path(row['isolatedPath'])]
    return list(dict.fromkeys(result))

def main(stage):
    global BUILD
    if stage in ('audit2','restart','restart_audit'):BUILD=BASE/'build02'
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;guard.size_bytes=lambda _:size()
    inputs=locks();kind='offline';seconds=120;cwd=BASE
    if stage=='prepare':cmd=[sys.executable,str(Path(__file__).resolve()),'task_prepare']
    elif stage=='build':
        assert read(BASE/'prepare01/report.json')['passed'];kind='build';seconds=360
        cmd=[sys.executable,str(Path(__file__).resolve()),'task_build']
    elif stage in ('control','restart'):
        assert read(BUILD/'report.json')['passed'];kind='query';seconds=30;cwd=BASE/(stage+'01')
        if stage=='restart':
            assert read(BASE/'audit02/report.json')['passed']
            from su2PhaseEulerContract import bind,read_restart
            origin=BASE/'control02';metadata=read(origin/'restart_flow.rah.json')
            expected=bind(origin/'restart_flow.csv',origin/'mesh.su2',origin/'case.cfg',BUILD/'SU2_CFD.exe',BUILD/'compatibility.exe','manufactured-control')
            read_restart(origin/'restart_flow.csv',metadata,expected,origin/'mesh.su2')
            cwd.mkdir(exist_ok=False)
            for n in ('mesh.su2','restart_flow.csv','restart_flow.rah.json'):shutil.copyfile(origin/n,cwd/n)
            cfg=(origin/'case.cfg').read_text(encoding='utf8');assert cfg.count('RESTART_SOL= NO')==1
            (cwd/'case.cfg').write_text(cfg.replace('RESTART_SOL= NO','RESTART_SOL= YES'),encoding='utf8')
        cmd=[str(BUILD/'compatibility.exe'),'case.cfg','.',stage]
        inputs += [BUILD/'compatibility.exe',BUILD/'SU2_CFD.exe',cwd/'case.cfg',cwd/'mesh.su2']
        if stage=='restart':inputs += [cwd/'restart_flow.csv',cwd/'restart_flow.rah.json',ROOT/'scripts/su2PhaseEulerContract.py']
    elif stage in ('build_audit','build_audit2'):
        script=ROOT/'scripts/verifyPhaseEulerBuildV229.py';cmd=[sys.executable,str(script)];inputs += [script]
    elif stage in ('audit','audit2','restart_audit'):
        script=ROOT/'scripts/verifyPhaseEulerCompatibilityV229.py';cmd=[sys.executable,str(script),stage]
        inputs += [script,ROOT/'scripts/su2PhaseEulerContract.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py']
    else:raise ValueError('NO_CFD_ADMISSION')
    return guard.supervise(stage,cmd,seconds,512,inputs,cwd,kind)

if __name__=='__main__':
    if CLOSED:raise SystemExit('v229 execution closed')
    if sys.argv[1]=='task_prepare':prep()
    elif sys.argv[1]=='task_build':build()
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
