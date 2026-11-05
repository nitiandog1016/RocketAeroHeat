"""固定变比热内嵌引擎，先实际源语法检查，再一次全量构建，最后恢复共享源码。"""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import subprocess
import zipfile

from validateIdealAir import run

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'tools/runtime/su2-thermal-build'
SOURCE=BASE/'source'
BUILD=BASE/'build-v117'
OUT=ROOT/'outputs/v154_validation/build01'
PATCH=ROOT/'outputs/v154_validation/patch02'
FLAVOR='rah-ideal-air-v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf8')


def main():
    assert not OUT.exists(); OUT.mkdir(parents=True)
    proof=json.loads((ROOT/'outputs/v154_validation/adapter02/report.json').read_text())
    assert proof['passed']
    archive=BASE/'su2.zip'
    assert sha(archive)=='c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    old=ROOT/'tools/runtime/su2/8.6.0-rah-axial-wall-v1'
    previous=json.loads((old/'runtime.json').read_text())
    permitted={str(Path(p).relative_to(SOURCE)).replace('\\','/'):h
        for p,h in previous['buildSourceSHA256'].items() if str(p).startswith(str(SOURCE))}
    originals={}; differences=[]
    with zipfile.ZipFile(str(archive)) as bundle:
        for name in bundle.namelist():
            if name.endswith('/'): continue
            relative='/'.join(name.split('/')[1:]); raw=bundle.read(name); path=SOURCE/relative
            assert path.is_file(),relative
            if path.read_bytes()!=raw:
                assert sha(path)==permitted.get(relative),relative
                differences.append(dict(path=relative,archiveSHA256=hashlib.sha256(raw).hexdigest(),currentSHA256=sha(path)))
            originals[relative]=raw
    dependencies={}
    for name in ('eigen','mel'):
        record=json.loads((BASE/(name+'_download.json')).read_text())
        assert sha(BASE/(name+'.zip'))==record['sha256']
        with zipfile.ZipFile(str(BASE/(name+'.zip'))) as bundle:
            count=0
            for entry in bundle.namelist():
                if entry.endswith('/'): continue
                relative='/'.join(entry.split('/')[1:])
                assert (SOURCE/'externals'/name/relative).read_bytes()==bundle.read(entry)
                count+=1
        dependencies[name]=dict(filesVerified=count,archiveSHA256=record['sha256'])
    save(OUT/'upstream_audit.json',dict(filesVerified=len(originals),differences=differences,dependencies=dependencies))
    patch=json.loads((PATCH/'patch_manifest.json').read_text())
    temporary={SOURCE/name:(SOURCE/name).read_bytes() for name in patch}
    reset=['Common/src/geometry/CGeometry.cpp','SU2_CFD/src/python_wrapper_structure.cpp']
    temporary.update({SOURCE/name:(SOURCE/name).read_bytes() for name in reset})
    headers=['CSU2IdealAir.hpp','SU2EnergyState.hpp','IdealAir.hpp']
    destinations=[SOURCE/'SU2_CFD/include/fluid'/name for name in headers]
    assert not any(path.exists() for path in destinations)
    for name in headers:
        assert sha(ROOT/'third_party/ideal_air'/name)==sha(ROOT/'outputs/v154_validation/adapter02/executed_sources/third_party/ideal_air'/name)
    graph=BUILD/'build.ninja'
    raw_graph=graph.read_bytes()
    options={item['name']:item['value'] for item in json.loads((BUILD/'meson-info/intro-buildoptions.json').read_text())}
    assert options['with-mpi']=='disabled' and options['with-omp'] and options['cpu-arch']=='x86-64'
    assert not any(options[name] for name in ('enable-coolprop','enable-cgns','enable-tecio'))
    compiler=BASE/'toolchain/w64devkit/bin/g++.exe'
    ninja=Path('D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe')
    save(OUT/'cache_before.json',{str(p.relative_to(BUILD)):dict(sha256=sha(p),bytes=p.stat().st_size)
        for p in BUILD.rglob('*') if p.is_file() and p.suffix in ('.obj','.a','.exe')})
    try:
        for name,identity in patch.items():
            assert sha(PATCH/name)==identity['patchedSHA256']
            assert hashlib.sha256(originals[name]).hexdigest()==identity['originalSHA256']
            (SOURCE/name).write_bytes((PATCH/name).read_bytes())
        for name in reset: (SOURCE/name).write_bytes(originals[name])
        for name,destination in zip(headers,destinations):
            destination.write_bytes((ROOT/'third_party/ideal_air'/name).read_bytes())
        lines=raw_graph.decode().replace('\r\n','\n').splitlines(True)
        for index,line in enumerate(lines):
            if re.match(r'build .*: cpp_COMPILER ',line):
                assert ' |' not in line
                lines[index]=line.rstrip()+' | RAH_FORCE_IDEAL_AIR_V154\n'
        custom=BUILD/'ideal-air-v154.ninja'; assert not custom.exists()
        custom.write_text(''.join(lines)+'\nbuild RAH_FORCE_IDEAL_AIR_V154: phony\n',encoding='utf8')
        shutil.copyfile(str(custom),str(OUT/custom.name))
        inputs=[Path(__file__).resolve(),ROOT/'scripts/prepareSU2IdealAirPatch.py',ROOT/'scripts/validateIdealAir.py',
            ROOT/'docs/validation/v154/PLAN.md',compiler,ninja,graph,custom]+[ROOT/'third_party/ideal_air'/name for name in headers]
        lock={str(path):sha(path) for path in inputs}
        tree={str(path):sha(path) for path in SOURCE.rglob('*') if path.is_file()}
        save(OUT/'source_tree.json',tree); save(OUT/'build_sources.json',lock)
        # 实际待编译翻译单元先进行语法检查，避免全量构建才发现调用接口错误。
        flags=['-std=c++17','-O3','-fopenmp','-ffast-math','-fno-finite-math-only','-fsyntax-only',
            '-I'+str(BUILD/'Common/src'),'-I'+str(SOURCE/'externals/eigen'),
            '-I'+str(SOURCE/'externals/mel'),'-I'+str(SOURCE/'externals/CLI11')]
        for name in ('Common/src/CConfig.cpp','SU2_CFD/src/solvers/CEulerSolver.cpp','SU2_CFD/src/solvers/CNSSolver.cpp'):
            run(OUT,'syntax_'+Path(name).stem,[str(compiler)]+flags+[str(SOURCE/name)],120)
        command=[str(ninja),'-C',str(BUILD),'-f',custom.name,'-j','2','SU2_CFD/src/SU2_CFD.exe']
        env=dict(os.environ,PATH=str(compiler.parent)+os.pathsep+os.environ['PATH'])
        planned=subprocess.check_output(command+['-n'],env=env).decode()
        (OUT/'planned.log').write_text(planned)
        assert 'Regenerating' not in planned and 'Compiling' in planned
        run(OUT,'compile',command,720)
        assert all(sha(Path(path))==expected for path,expected in tree.items()),'Build source changed'
        assert all(sha(Path(path))==expected for path,expected in lock.items()),'Build input changed'
        actual=(OUT/'compile.log').read_text(errors='replace')
        objects=[path.rstrip('\r') for path in re.findall(r'Compiling C\+\+ object (.*)',planned)]
        assert objects and all('Compiling C++ object '+path in actual for path in objects)
        members={}
        for name in ('Common/src/libSU2Common.a','SU2_CFD/src/libSU2core.a'):
            actual_members=subprocess.check_output([str(compiler.parent/'ar.exe'),'t',str(BUILD/name)],cwd=str(BUILD)).decode().splitlines()
            expected=next(line for line in raw_graph.decode().splitlines() if line.startswith('build '+name+': ')).split(': ',1)[1].split()[1:]
            canonical=lambda path:str((BUILD/path).resolve()).lower()
            assert sorted(map(canonical,actual_members))==sorted(map(canonical,expected))
            members[name]=len(actual_members)
        executable=BUILD/'SU2_CFD/src/SU2_CFD.exe'
        imports=subprocess.check_output([str(compiler.parent/'objdump.exe'),'-p',str(executable)]).decode()
        dlls=re.findall(r'DLL Name:\s*(\S+)',imports)
        assert set(name.lower() for name in dlls)=={'kernel32.dll','msvcrt.dll'}
        destination=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
        assert not destination.exists(); (destination/'bin').mkdir(parents=True)
        shutil.copyfile(str(executable),str(destination/'bin/SU2_CFD.exe'))
        for item in previous['files']:
            if item['path'].endswith(('.exe','.diff')): continue
            path=old/item['path']; assert sha(path)==item['sha256']
            target=destination/item['path']; target.parent.mkdir(parents=True,exist_ok=True)
            os.link(str(path),str(target))
        shutil.copyfile(str(PATCH/'ideal_air.diff'),str(destination/'ideal_air.diff'))
        shutil.copyfile(str(ROOT/'third_party/ideal_air/LICENSE.CoolProp'),str(destination/'LICENSE.CoolProp'))
        runtime=dict(flavor=FLAVOR,version='8.6.0+'+FLAVOR,status='compiled-tests-pending',
            executable=str(destination/'bin/SU2_CFD.exe'),sourceCommit=previous['sourceCommit'],archiveSHA256=sha(archive),
            thermodynamicModel='Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1',temperatureRangeK=[60,2000],
            gasConstant=8.31451/.02896546,physicsQualified=False,systemDLLs=dlls,
            allCFDObjectsRebuilt=len(objects),archiveMembers=members,buildSourceSHA256=lock,
            sourceTreeSHA256=sha(OUT/'source_tree.json'),buildReceipt=json.loads((OUT/'compile_process.json').read_text()),
            restrictions=['direct SI steady single-zone NS/SST','general Roe','resolved wall','static PT Riemann inlet and pressure outlet','research interface until validated'],
            files=[dict(path=path.relative_to(destination).as_posix(),bytes=path.stat().st_size,sha256=sha(path))
                for path in sorted(destination.rglob('*')) if path.is_file()])
        save(destination/'runtime.json',runtime); save(OUT/'runtime_installed.json',runtime)
        print(json.dumps(dict(installed=str(destination),compiledObjects=len(objects),physicsQualified=False)),flush=True)
    finally:
        for path,raw in temporary.items():
            target=OUT/'executed_sources'/path.relative_to(SOURCE)
            target.parent.mkdir(parents=True,exist_ok=True); target.write_bytes(path.read_bytes())
            path.write_bytes(raw)
        for path in destinations:
            if path.exists():
                target=OUT/'executed_sources'/path.relative_to(SOURCE)
                target.parent.mkdir(parents=True,exist_ok=True); target.write_bytes(path.read_bytes())
                path.resolve().relative_to(SOURCE.resolve())
                path.unlink()
        restored={str(path):sha(path)==hashlib.sha256(raw).hexdigest() for path,raw in temporary.items()}
        save(OUT/'source_restore.json',restored); assert all(restored.values())


if __name__=='__main__':
    main()
