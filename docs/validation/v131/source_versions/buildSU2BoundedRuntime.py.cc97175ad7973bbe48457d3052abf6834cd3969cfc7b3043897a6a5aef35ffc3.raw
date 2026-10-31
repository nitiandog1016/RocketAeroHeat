"""v131受控有界输运构建：全量重建原CFD目标，恢复源文件，不覆盖任何已安装运行时。"""
from pathlib import Path
import hashlib,json,os,re,shutil,subprocess,zipfile
import buildSU2ThermalRuntime as supervisor
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'tools/runtime/su2-thermal-build'
SOURCE=BASE/'source';BUILD=BASE/'build-v117';OUT=ROOT/'outputs/v131_validation/build_final'
FLAVOR='rah-bounded-conservative-v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def main():
    assert not OUT.exists();OUT.mkdir(parents=True);supervisor.OUT=OUT
    archive=BASE/'su2.zip'
    assert sha(archive)=='c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    old=ROOT/'tools/runtime/su2/8.6.0-rah-axial-wall-v1'
    previous=json.loads((old/'runtime.json').read_text())
    permitted={str(Path(p).relative_to(SOURCE)).replace('\\','/'):h for p,h in previous['buildSourceSHA256'].items() if str(p).startswith(str(SOURCE))}
    originals={};differences=[]
    # 模块一：上游全部归档文件和两个实际头文件依赖逐字核对。
    with zipfile.ZipFile(str(archive)) as z:
        for n in z.namelist():
            if n.endswith('/'):continue
            rel='/'.join(n.split('/')[1:]);raw=z.read(n);p=SOURCE/rel
            assert p.is_file(),rel
            if p.read_bytes()!=raw:
                assert sha(p)==permitted.get(rel),rel
                differences.append(dict(path=rel,archiveSHA256=hashlib.sha256(raw).hexdigest(),currentSHA256=sha(p)))
            originals[rel]=raw
    dependency={}
    for name in ['eigen','mel']:
        dep_archive=BASE/(name+'.zip');record=json.loads((BASE/(name+'_download.json')).read_text())
        assert sha(dep_archive)==record['sha256']
        count=0
        with zipfile.ZipFile(str(dep_archive)) as z:
            for n in z.namelist():
                if n.endswith('/'):continue
                rel='/'.join(n.split('/')[1:])
                assert (SOURCE/'externals'/name/rel).read_bytes()==z.read(n),rel
                count+=1
        dependency[name]=dict(archiveSHA256=sha(dep_archive),filesVerified=count,provenance=record)
    write(OUT/'upstream_audit.json',dict(filesVerified=len(originals),differences=differences,dependencies=dependency))
    # 模块二：剔除旧热边界试验；仅保留三项包含路径修复及本次四行Jacobian补丁。
    temporary={}
    for rel in ['Common/src/geometry/CGeometry.cpp','SU2_CFD/src/python_wrapper_structure.cpp',
                'SU2_CFD/include/solvers/CScalarSolver.hpp','SU2_CFD/include/solvers/CScalarSolver.inl']:
        p=SOURCE/rel;temporary[p]=p.read_bytes()
    proof=json.loads((ROOT/'outputs/v131_validation/preflight/report.json').read_text())
    assert proof['passed']
    graph=BUILD/'build.ninja';graphraw=graph.read_bytes()
    options=json.loads((BUILD/'meson-info/intro-buildoptions.json').read_text())
    opt={v['name']:v['value'] for v in options}
    assert opt['with-mpi']=='disabled' and opt['with-omp'] and opt['cpu-arch']=='x86-64'
    assert not any(opt[n] for n in ['enable-coolprop','enable-cgns','enable-tecio'])
    cc=BASE/'toolchain/w64devkit/bin';ninja=Path('D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe')
    env=dict(os.environ,PATH=str(cc)+os.pathsep+str(ninja.parent)+os.pathsep+os.environ['PATH'])
    cache_before={str(p.relative_to(BUILD)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BUILD.rglob('*') if p.is_file()}
    write(OUT/'cache_before.json',cache_before)
    try:
        for rel in ['Common/src/geometry/CGeometry.cpp','SU2_CFD/src/python_wrapper_structure.cpp']:
            (SOURCE/rel).write_bytes(originals[rel])
        for rel in ['SU2_CFD/include/solvers/CScalarSolver.hpp','SU2_CFD/include/solvers/CScalarSolver.inl']:
            p=SOURCE/rel;candidate=ROOT/'outputs/v131_validation/preflight'/(p.name+'.patched')
            assert p.read_bytes()==originals[rel]
            p.write_bytes(candidate.read_bytes())
        # 原图只追加无文件的强制依赖，所有CFD链接输入重编，避免旧FrozenEnergy对象混入。
        lines=graphraw.decode().replace('\r\n','\n').splitlines(True);objects=[]
        for i,line in enumerate(lines):
            if re.match(r'build .*: cpp_COMPILER ',line):
                assert ' |' not in line
                objects.append(line.split(': cpp_COMPILER')[0][6:])
                lines[i]=line.rstrip()+' | RAH_FORCE_V131\n'
        custom=BUILD/'bounded-v131-final.ninja';assert not custom.exists()
        custom.write_text(''.join(lines)+'\nbuild RAH_FORCE_V131: phony\n',encoding='utf-8')
        shutil.copyfile(str(custom),str(OUT/'bounded-v131.ninja'))
        tree={str(p):sha(p) for p in SOURCE.rglob('*') if p.is_file()}
        write(OUT/'source_tree.json',tree)
        locked={str(p):sha(p) for p in [Path(__file__).resolve(),ninja,cc/'g++.exe',custom,graph,
            ROOT/'third_party/su2_bounded/bounded_conservative_jacobian.diff',ROOT/'docs/validation/v131/JACOBIAN_PREFLIGHT.md']}
        write(OUT/'build_sources.json',locked)
        command=[str(ninja),'-C',str(BUILD),'-f',custom.name,'-j','2','SU2_CFD/src/SU2_CFD.exe']
        planned=subprocess.check_output(command+['-n'],env=env).decode()
        (OUT/'planned.log').write_text(planned,encoding='utf-8')
        assert 'Regenerating' not in planned and 'Compiling' in planned
        supervisor.run('compile',command,SOURCE,env,720)
        assert all(sha(Path(p))==h for p,h in tree.items()),'构建期间源码变化'
        assert all(sha(Path(p))==h for p,h in locked.items()),'构建输入变化'
        actual=(OUT/'compile.log').read_text(errors='replace')
        expected_objects=re.findall(r'Compiling C\+\+ object (.*)',planned)
        assert expected_objects and all('Compiling C++ object '+p in actual for p in expected_objects)
        # 薄归档的成员与原图逐项匹配，不只依据链接成功判断来源。
        members={}
        for lib in ['Common/src/libSU2Common.a','SU2_CFD/src/libSU2core.a']:
            native=subprocess.check_output([str(cc/'ar.exe'),'t',str(BUILD/lib)],cwd=str(BUILD),env=env).decode().splitlines()
            expected=next(l for l in graphraw.decode().splitlines() if l.startswith('build '+lib+': ')).split(': ',1)[1].split()[1:]
            canonical=lambda p:str(Path(p).resolve()).lower()
            actual_paths=[canonical(BUILD/p) for p in native]
            expected_paths=[canonical(BUILD/p) for p in expected]
            assert sorted(actual_paths)==sorted(expected_paths),(lib,len(native),len(expected))
            members[lib]=dict(memberCount=len(native),exactMembers=True)
        exe=BUILD/'SU2_CFD/src/SU2_CFD.exe'
        imports=subprocess.check_output([str(cc/'objdump.exe'),'-p',str(exe)],env=env).decode()
        dlls=re.findall(r'DLL Name:\s*(\S+)',imports)
        assert set(n.lower() for n in dlls)=={'kernel32.dll','msvcrt.dll'}
        dest=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR);assert not dest.exists();(dest/'bin').mkdir(parents=True)
        shutil.copyfile(str(exe),str(dest/'bin/SU2_CFD.exe'))
        for entry in previous['files']:
            if entry['path'].endswith('.exe') or entry['path'].endswith('.diff'):continue
            p=old/entry['path'];assert sha(p)==entry['sha256'];q=dest/entry['path'];q.parent.mkdir(parents=True,exist_ok=True);os.link(str(p),str(q))
        patch=ROOT/'third_party/su2_bounded/bounded_conservative_jacobian.diff'
        shutil.copyfile(str(patch),str(dest/patch.name))
        runtime=dict(flavor=FLAVOR,version='8.6.0+'+FLAVOR,status='compiled-tests-pending',
            executable=str(dest/'bin/SU2_CFD.exe'),archiveSHA256=sha(archive),sourceCommit=previous['sourceCommit'],
            license=previous['license'],physicsQualified=False,systemDLLs=dlls,buildReceipt=json.loads((OUT/'compile_process.json').read_text()),
            sourceTreeSHA256=sha(OUT/'source_tree.json'),buildSourceSHA256=locked,
            correction='conservative bounded scalar divergence Jacobian divided by density',
            patchSHA256=sha(patch),priorThermalAndFrozenHooksRemoved=True,allCFDObjectsRebuilt=True,archiveMembers=members,
            restrictions=['steady single-zone ideal-gas SST-2003m','SI dimensional','constant isothermal wall',
                'MGLEVEL=0','explicit controlled study only','no HEOS or SA qualification'],
            files=[dict(path=p.relative_to(dest).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(dest.rglob('*')) if p.is_file()])
        write(dest/'runtime.json',runtime);write(OUT/'runtime_installed.json',runtime)
        print(json.dumps(dict(installed=str(dest),allObjectsRebuilt=len(expected_objects),physicsQualified=False)),flush=True)
    finally:
        for p,raw in temporary.items():
            (OUT/(p.name+'.used')).write_bytes(p.read_bytes());p.write_bytes(raw)
        restored={str(p):sha(p)==hashlib.sha256(raw).hexdigest() for p,raw in temporary.items()}
        write(OUT/'source_restore.json',restored);assert all(restored.values())
        after={str(p.relative_to(BUILD)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in BUILD.rglob('*') if p.is_file()}
        write(OUT/'cache_after.json',after)
        write(OUT/'cache_delta.json',dict(lengthBefore=sum(r['bytes'] for r in cache_before.values()),
            lengthAfter=sum(r['bytes'] for r in after.values()),changed=[p for p in after if after[p]!=cache_before.get(p)],
            sourceRestored=True,installedHistoricalRuntimesNotModified=True))
if __name__=='__main__':main()

