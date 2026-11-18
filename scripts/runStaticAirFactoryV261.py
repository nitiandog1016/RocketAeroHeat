"""从已校验的 SU2 对象集构建隔离后端，监督静止空气配置对照。"""
import difflib
import shutil
import subprocess
import sys
from pathlib import Path

import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT, read, save, sha
from terminalCapacityV249 import supervise

BASE = ROOT/'outputs/v261_validation'
DOC = ROOT/'docs/validation/v261'
BUILD = BASE/'build'
OVERLAY = BASE/'source'
V220 = ROOT/'outputs/v220_validation'
V235 = ROOT/'outputs/v235_validation'
V257 = ROOT/'outputs/v257_validation'
V260 = ROOT/'outputs/v260_validation'
WD = ROOT/'tools/runtime/su2-thermal-build/build-v117'
NATIVE = ROOT/'native/su2PhaseAwareAir'
FIXTURE = ROOT/'tests/data/static_air_v261/default.cfg'
GAS_HEADERS = ('CaloricOptionV261.hpp','CaloricFactoryV261.hpp','StableCaloricInverseV260.hpp','CSU2PhaseAwareAirStableV260.hpp')
SOURCES = {
    'CConfig': (V220/'source/Common/src/CConfig.cpp','241eb647b07a136eb87f5528960db0737ae933ef5e45579be97218d7d692a48e'),
    'CEulerSolver': (V220/'source/SU2_CFD/src/solvers/CEulerSolver.cpp','1d079fc23f96421b752333a46513a41d3d48f4785a6f34863fed57495581fff0'),
    'CNSVariable': (V235/'stage01/CNSVariable.cpp','e7b284dd4666685e9d91d28ff1a6a53cd2216815ea52c36bfaba045856400de7'),
}
VALID = ('default','stable')
INVALID = {'unknown':'RAH_CALORIC_INVERSE expects STABLE_PIECEWISE',
           'incompatible':'RAH_CALORIC_INVERSE requires FLUID_MODEL=RAH_PHASE_AIR'}
MANIFESTS = [V220/'build02/dependencies.json', V235/'stage01/dependencies.json',
             V260/'build02/dependencies.json', V257/'build01/members.json']
COMPILE_INPUTS = [V220/'process_compile0/admission.json',V220/'process_recover_compile1/admission.json',
                  V235/'stage01/commands.json',V260/'build02/command.json',V220/'build02/link_command.json']


def change_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('Upstream source context changed: '+old[:90])
    return text.replace(old,new)


def patched(unit, text):
    """只改变配置登记、两个工厂调用点及原生异常识别。"""
    if unit=='CConfig':
        text=change_once(text,'#undef ENABLE_MAPS','#undef ENABLE_MAPS\n#include "CaloricOptionV261.hpp"')
        anchor='  addEnumOption("FLUID_MODEL", Kind_FluidModel, FluidModel_Map, STANDARD_AIR);'
        text=change_once(text,anchor,anchor+'\n  all_options.emplace("RAH_CALORIC_INVERSE", true);\n'
                         '  option_map.emplace("RAH_CALORIC_INVERSE", new CStableCaloricOption());')
        return change_once(text,'  if (Kind_FluidModel == RAH_PHASE_AIR) {',
            '  if (OptionIsSet("RAH_CALORIC_INVERSE") && Kind_FluidModel != RAH_PHASE_AIR)\n'
            '    SU2_MPI::Error("RAH_CALORIC_INVERSE requires FLUID_MODEL=RAH_PHASE_AIR", CURRENT_FUNCTION);\n'
            '  if (Kind_FluidModel == RAH_PHASE_AIR) {')
    if unit=='CEulerSolver':
        anchor='#include "../../include/fluid/CSU2PhaseAwareAir.hpp"'
        text=change_once(text,anchor,anchor+'\n#include "CaloricFactoryV261.hpp"')
        text=change_once(text,'      auxFluidModel = new CSU2PhaseAwareAir();',
             '      auxFluidModel = rah_caloric_factory::create(*config);\n'
             '      if (rank == MASTER_NODE && MGLevel == MESH_0)\n'
             '        cout << "RAH caloric initialization: " << rah_caloric_factory::method(auxFluidModel) << endl;')
        text=change_once(text,'        FluidModel[thread] = new CSU2PhaseAwareAir();',
                              '        FluidModel[thread] = rah_caloric_factory::create(*config);')
        return change_once(text,'ModelTable << "RAH_PHASE_AIR 60..2000K phase domain e+300000";',
                          'ModelTable << "RAH_PHASE_AIR 60..2000K e+300000 " << rah_caloric_factory::method(GetFluidModel());')
    if unit!='CNSVariable':
        raise ValueError(unit)
    anchor='#include <array>'
    text=change_once(text,anchor,'#include "CaloricFactoryV261.hpp"\n'+anchor)
    text=change_once(text,'dynamic_cast<CSU2PhaseAwareAir*>(FluidModel) == nullptr',
                          '!rah_caloric_factory::isPhaseAir(FluidModel)')
    return change_once(text,'dynamic_cast<CSU2PhaseAwareAir*>(FluidModel) != nullptr',
                            'rah_caloric_factory::isPhaseAir(FluidModel)')


def prepare():
    OVERLAY.mkdir(exist_ok=False)
    BUILD.mkdir(exist_ok=False)
    changes=[]
    for unit,(path,expected) in SOURCES.items():
        assert sha(path)==expected,path
        old=path.read_text(encoding='utf8')
        new=patched(unit,old)
        (OVERLAY/path.name).write_text(new,encoding='utf8')
        changes.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),
                                           fromfile=unit+'.cpp (original)',tofile=unit+'.cpp (v261)'))
    (OVERLAY/'factory.patch').write_text(''.join(changes),encoding='utf8')
    for name in GAS_HEADERS+('staticAirFactoryV261.cpp',):
        shutil.copyfile(NATIVE/name,OVERLAY/name)
    default=FIXTURE.read_text(encoding='utf8')
    stable=default+'\nRAH_CALORIC_INVERSE= STABLE_PIECEWISE\n'
    cases={'default':default,'stable':stable,'unknown':stable.replace('STABLE_PIECEWISE','UNRECOGNIZED'),
           'incompatible':stable.replace('FLUID_MODEL= RAH_PHASE_AIR','FLUID_MODEL= STANDARD_AIR')}
    for name,config in cases.items():
        folder=BASE/name
        folder.mkdir(exist_ok=False)
        (folder/'case.cfg').write_text(config,encoding='utf8')
    commands={'CConfig':read(COMPILE_INPUTS[0])['argv'],'CEulerSolver':read(COMPILE_INPUTS[1])['argv'],
              'CNSVariable':read(COMPILE_INPUTS[2])['SU2_CFD/src/variables/CNSVariable.cpp']}
    for unit,cmd in commands.items():
        source=SOURCES[unit][0]
        assert str(source) in cmd,unit
        cmd[cmd.index(str(source))]=str(OVERLAY/source.name)
        cmd[cmd.index('-o')+1]=str(BUILD/(unit+'.obj'))
        cmd[cmd.index('-MF')+1]=str(BUILD/(unit+'.d'))
        cmd[1:1]=['-iquote',str(source.parent),'-I'+str(OVERLAY),'-I'+str(V220/'source/SU2_CFD/include/fluid')]
    normal=read(COMPILE_INPUTS[4])
    normal[normal.index('-o')+1]=str(BUILD/'SU2_CFD.exe')
    inspector=read(COMPILE_INPUTS[3])
    inspector[inspector.index(str(NATIVE/'caloricReplayV260.cpp'))]=str(OVERLAY/'staticAirFactoryV261.cpp')
    inspector[inspector.index('-MF')+1]=str(BUILD/'inspect.d')
    inspector[inspector.index('-o')+1]=str(BUILD/'inspect.exe')
    for cmd in (normal,inspector):
        for i,part in enumerate(cmd):
            if part.endswith('.a'):
                cmd[i]=str(BUILD/Path(part).name)
    save(BASE/'commands.json',dict(compile=commands,production=normal,inspector=inspector))
    save(BASE/'prepared.json',dict(source={k:h for k,(_,h) in SOURCES.items()},
         overlay={p.name:sha(p) for p in OVERLAY.iterdir()},configs={n:sha(BASE/n/'case.cfg') for n in cases}))


def dependencies(path):
    contents=path.read_text(encoding='utf8').replace('\\\n',' ').split(': ',1)[1]
    return [(WD/Path(token)).resolve() for token in contents.split()]


def build():
    commands=read(BASE/'commands.json')
    expected={}
    for path in MANIFESTS[:3]:
        for name,h in read(path).items():
            if name in expected and expected[name]!=h:
                raise ValueError('Conflicting dependencies: '+name)
            expected[name]=h
    expected.update({str(OVERLAY/name):h for name,h in read(BASE/'prepared.json')['overlay'].items()})
    def check_dependencies(filename):
        hashes={str(p):sha(p) for p in dependencies(BUILD/filename)}
        for name,h in hashes.items():
            assert expected.get(name)==h,('Unbound dependency',name)
        save(BUILD/(filename+'.json'),hashes)
    for name,cmd in commands['compile'].items():
        subprocess.check_call(cmd,cwd=WD)
        check_dependencies(name+'.d')
        save(BUILD/(name+'-built.json'),dict(objectSHA256=sha(BUILD/(name+'.obj'))))
    ar=Path(commands['compile']['CConfig'][0]).with_name('ar.exe')
    members=read(MANIFESTS[3])
    replacements=[]
    for lib,rows in members.items():
        paths=[]
        for row in rows:
            original=Path(row['path'])
            if original.stem in SOURCES:
                replacements.append(original.stem)
                paths.append(BUILD/original.name)
            else:
                paths.append(original)
        arguments=BUILD/(lib+'.rsp')
        arguments.write_text(''.join('"'+p.as_posix()+'"\n' for p in paths),encoding='utf8')
        subprocess.check_call([str(ar),'crsT',str(BUILD/lib),'@'+str(arguments)],cwd=WD)
    assert sorted(replacements)==sorted(SOURCES),replacements
    for entry in ('production','inspector'):
        subprocess.check_call(commands[entry],cwd=WD)
    check_dependencies('inspect.d')
    save(BUILD/'report.json',dict(rebuiltUnits=replacements,
         reusedObjects=sum(len(v) for v in members.values())-len(replacements),
         executables={p.name:sha(p) for p in BUILD.glob('*.exe')}))


def locks_for(stage):
    locks=[Path(__file__).resolve(),DOC/'PLAN.md',FIXTURE,Path(guard.__file__).resolve(),
           ROOT/'scripts/terminalCapacityV249.py',*(NATIVE/n for n in GAS_HEADERS+('staticAirFactoryV261.cpp',)),
           *(row[0] for row in SOURCES.values()),*MANIFESTS,*COMPILE_INPUTS]
    if stage=='prepare':
        return locks
    preparation=read(BASE/'prepared.json')
    for n,h in preparation['overlay'].items():
        p=OVERLAY/n
        assert sha(p)==h,p
        locks.append(p)
    locks.extend([BASE/'commands.json',BASE/'prepared.json'])
    for name,h in preparation['configs'].items():
        p=BASE/name/'case.cfg'
        assert sha(p)==h,p
        locks.append(p)
    if stage=='build':
        previous={}
        for path in MANIFESTS[:3]:
            previous.update(read(path))
        for name,h in previous.items():
            assert sha(name)==h,name
            locks.append(Path(name))
        for rows in read(MANIFESTS[3]).values():
            for row in rows:
                assert sha(row['path'])==row['sha256'],row['path']
                locks.append(Path(row['path']))
        production=read(BASE/'commands.json')['production']
        locks.extend([Path(production[0]),Path(production[0]).with_name('ar.exe'),
                      next(Path(x) for x in production if x.endswith('.obj'))])
    else:
        for n,h in read(BUILD/'report.json')['executables'].items():
            assert sha(BUILD/n)==h,n
            locks.append(BUILD/n)
        locks.append(BUILD/'report.json')
    return list(dict.fromkeys(locks))


def main(stage):
    BASE.mkdir(exist_ok=True)
    guard.BASE=guard.CASE=BASE
    guard.DOC=DOC
    size=guard.size_bytes
    guard.size_bytes=lambda p: size(BASE) if p.name.startswith('process_') else 0
    locks=locks_for(stage)
    snap=BASE/('executed_'+stage)
    snap.mkdir(exist_ok=False)
    for source in [Path(__file__).resolve(),DOC/'PLAN.md',*(NATIVE/n for n in GAS_HEADERS+('staticAirFactoryV261.cpp',))]:
        shutil.copyfile(source,snap/source.name)
    if stage in ('prepare','build'):
        argv=[sys.executable,str(Path(__file__).resolve()),'task_'+stage]
        limit,kind=(60,'offline') if stage=='prepare' else (240,'build')
        cwd=ROOT
    elif stage in VALID or stage in INVALID:
        exe='SU2_CFD.exe' if stage in INVALID else 'inspect.exe'
        argv=[str(BUILD/exe),'case.cfg']
        limit,kind,cwd=12,'control',BASE/stage
    else:
        raise ValueError(stage)
    ok=supervise(guard,BASE,stage,argv,limit,160,locks,cwd,kind)
    if stage in INVALID:
        receipt=read(BASE/('process_'+stage)/'report.json')
        log=(BASE/('process_'+stage)/'stdout.txt').read_text(encoding='utf8',errors='replace')
        return (not ok and receipt['exitCode']!=0 and receipt['stopReason'] is None
                and receipt['inputsStable'] and not receipt['liveGroupAfterCleanup']
                and INVALID[stage] in log and 'RAH caloric initialization:' not in log)
    return ok


if __name__=='__main__':
    stage=sys.argv[1]
    if stage=='task_prepare':
        prepare()
    elif stage=='task_build':
        build()
    else:
        raise SystemExit(0 if main(stage) else 1)
