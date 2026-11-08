"""模块职责：唯一黏度对象适配恢复；保留初次对象和引擎，仅写独立build02。"""
from pathlib import Path
import json,shutil,subprocess,sys
from prepareSU2PhaseAirV220 import ROOT,BASE,SOURCE,BUILD,OLD_BUILD,BOOST,sha,read,save
from buildSU2PhaseAirV220 import compile_commands
from runSU2PhaseAirV220 import supervise
NEW=BASE/'build02'

def task(stage):
    if stage=='prepare':
        NEW.mkdir(exist_ok=False)
        old=SOURCE/'SU2_CFD/include/fluid/CSU2PhaseAwareAir.hpp'
        assert sha(old)==sha(BASE/'helper02_source/CSU2PhaseAwareAir.hpp')
        assert old.stat().st_nlink==1
        shutil.copyfile(ROOT/'native/su2PhaseAwareAir/CSU2PhaseAwareAir.hpp',old)
        deps=read(BASE/'dependencies01/inputs.json')
        deps[str(old.resolve())]=sha(old)
        assert all(sha(Path(n))==h for n,h in deps.items())
        save(NEW/'dependencies.json',deps)
        return
    if stage=='link':
        members=read(BASE/'link01/members.json');updated={}
        for name,rows in members.items():
            paths=[];out=[]
            for row in rows:
                p=Path(row['path']);assert sha(p)==row['sha256']
                if p.name in ('CEulerSolver.obj','CNSVariable.obj'):p=NEW/p.name
                paths.append(str(p));out.append(dict(path=str(p),sha256=sha(p)))
            subprocess.check_call([str(Path(compile_commands()[0][0]).parent/'ar.exe'),'crsT',str(NEW/name)]+paths)
            updated[name]=out
        command=read(BASE/'link01/command.json')
        command=[str(NEW/Path(x).name) if x in (str(BUILD/'SU2_CFD.exe'),str(BUILD/'libSU2core.a'),str(BUILD/'libSU2Common.a')) else x for x in command]
        save(NEW/'members.json',updated);save(NEW/'link_command.json',command)
        subprocess.check_call(command,cwd=str(OLD_BUILD))
        save(NEW/'link_report.json',dict(passed=True,executableSHA256=sha(NEW/'SU2_CFD.exe'),nativeVerification=False,CFDRuns=0))
        return
    raise ValueError(stage)

def main(stage):
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2PhaseAirV220.py',ROOT/'docs/validation/v220/PLAN.md',
           ROOT/'native/su2PhaseAwareAir/CSU2PhaseAwareAir.hpp',ROOT/'native/su2PhaseAwareAir/probe.cpp']
    if stage=='prepare':
        return supervise('recover_prepare',[sys.executable,__file__,'task','prepare'],'offline',30,512,locks,BASE)
    deps=read(NEW/'dependencies.json');assert all(sha(Path(n))==h for n,h in deps.items())
    locks += [Path(n) for n in deps]+[NEW/'dependencies.json']
    if stage in ('compile1','compile3'):
        command=compile_commands()[int(stage[-1])]
        for flag in ('-o','-MF'):
            i=command.index(flag)+1;command[i]=str(NEW/Path(command[i]).name)
        return supervise('recover_'+stage,command,'build',120,256,locks,BASE)
    if stage=='link':
        for name in ('compile1','compile3'):assert read(BASE/('process_recover_'+name)/'report.json')['passed']
        locks += [NEW/'CEulerSolver.obj',NEW/'CNSVariable.obj']
        return supervise('recover_link',[sys.executable,__file__,'task','link'],'build',30,256,locks,BASE)
    if stage=='helper':
        command=read(BASE/'process_helper2/report.json')['argv']
        command=[str(NEW/Path(x).name) if x in (str(BUILD/'probe02.d'),str(BUILD/'phase_probe02.exe'),str(BUILD/'libSU2core.a'),str(BUILD/'libSU2Common.a')) else x for x in command]
        # 小控制保留行号，便于失败时只读定位；不改变优化或物理算式。
        command.insert(1,'-g')
        locks += [NEW/'libSU2core.a',NEW/'libSU2Common.a']
        return supervise('recover_helper',command,'build',120,256,locks,BASE)
    if stage=='query':
        assert read(BASE/'process_recover_helper/report.json')['passed']
        out=BASE/'control03';out.mkdir(exist_ok=False);locks += [NEW/'phase_probe02.exe']
        return supervise('recover_query',[str(NEW/'phase_probe02.exe'),str(out)],'query',30,32,locks,out)
    raise ValueError('Stage not admitted')

if __name__=='__main__':
    if len(sys.argv)>2 and sys.argv[1]=='task':task(sys.argv[2])
    else:raise SystemExit(0 if main(sys.argv[1]) else 1)
