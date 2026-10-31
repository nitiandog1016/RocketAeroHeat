"""固定SU2源码有界诊断/恢复候选构建；恢复已留归档，不下载或覆盖正式运行时。"""
from pathlib import Path
import hashlib,json,os,sys,subprocess,zipfile,difflib,shutil
import buildSU2ThermalRuntime as build
from prepareSU2ThermalBuild import extract
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'tools/runtime/su2-thermal-build'
SOURCE=BASE/'source';OUT=ROOT/'outputs/v117_validation';FLAVOR='rah-frozen-energy-v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(parents=True,exist_ok=True);build.OUT=OUT
    compiler=BASE/'compiler.7z.exe'
    assert sha(compiler)=='54c2bb9962e8a30c8270e34c5efc3d5fee579f53087fa7d435b95fbba2af7b2b'
    bindir=BASE/'toolchain/w64devkit/bin'
    if not (bindir/'g++.exe').exists():
        (BASE/'toolchain').mkdir(exist_ok=True)
        build.run('extract_compiler',['D:/Bandizip/bz.exe','x','-y','-o:'+str(BASE/'toolchain'),str(compiler)],BASE,os.environ,120)
    assert sha(BASE/'meson.zip')=='4e2b6d43a766cc7f8ea2922339b358f9f325cc8e3dd406debb16a83382e9b532'
    extract('meson',BASE/'meson')
    ninja=Path('D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe')
    env=dict(os.environ,PATH=str(bindir)+os.pathsep+str(ninja.parent)+os.pathsep+os.environ['PATH'],
             CC=str(bindir/'gcc.exe'),CXX=str(bindir/'g++.exe'))
    build.run('unit_compile',[str(bindir/'g++.exe'),'-std=c++17','-O2','-Wall','-Wextra',
        str(ROOT/'tests/testSU2FrozenEnergyRecovery.cpp'),'-o',str(OUT/'unit.exe')],ROOT,env,60)
    build.run('unit_run',[str(OUT/'unit.exe')],ROOT,env,30)
    src=SOURCE/'SU2_CFD/src/solvers/CEulerSolver.cpp'
    assert sha(BASE/'su2.zip')=='c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    with zipfile.ZipFile(str(BASE/'su2.zip')) as z:
        original=z.read(next(n for n in z.namelist() if n.endswith('/SU2_CFD/src/solvers/CEulerSolver.cpp'))).decode('utf-8')
    assert hashlib.sha256(original.encode()).hexdigest()=='a02047d187aeb4dcd5dffadbce9c21243ace41df79f806e58091b99b3465dc20'
    headers=[ROOT/'third_party/su2_thermal'/n for n in ['RAHFrozenEnergyRecovery.hpp','RAHFrozenEnergyDiagnostics.hpp']]
    hook=ROOT/'third_party/su2_thermal/FrozenEnergyHook.cpp.inc'
    begin=original.index('void CEulerSolver::ComputeUnderRelaxationFactor(')
    end=original.index('\nvoid CEulerSolver::IdentifySolutionOutliers',begin)
    block=original[begin:end]
    setup='''
  const char* rah_mode=std::getenv("RAH_FROZEN_ENERGY_MODE");
  const bool rah_enabled=rah_mode && *rah_mode;
  const bool rah_apply=rah_enabled && std::string(rah_mode)=="apply";
  if (rah_enabled && !rah_apply && std::string(rah_mode)!="diagnostic")
    SU2_MPI::Error("Unknown frozen-energy mode",CURRENT_FUNCTION);
  unsigned long rah_checked=0,rah_frozen=0,rah_eligible=0,rah_applied=0;
'''
    block=block.replace('  SU2_OMP_FOR_STAT(omp_chunk_size)',setup+'\n  SU2_OMP_FOR_STAT(omp_chunk_size)',1)
    block=block.replace('    nodes->SetUnderRelaxation',hook.read_text(encoding='utf-8')+'\n    nodes->SetUnderRelaxation',1)
    block=block.replace('  END_SU2_OMP_FOR\n','''  END_SU2_OMP_FOR
  if (rah_enabled) {
    SU2_OMP_CRITICAL {
      rah_frozen_log().count(config->GetInnerIter(),rah_checked,rah_frozen,rah_eligible,rah_applied);
    } END_SU2_OMP_CRITICAL
  }
''',1)
    prefix='''#include <cstdlib>
#include "RAHFrozenEnergyDiagnostics.hpp"
static rah::FrozenEnergyLog& rah_frozen_log() { static rah::FrozenEnergyLog logger; return logger; }
'''
    patched=prefix+original[:begin]+block+original[end:]
    assert src.read_bytes() in (original.encode(),patched.encode()),'未记录源码改动，拒绝覆盖'
    src.write_bytes(patched.encode())
    for h in headers:(src.parent/h.name).write_bytes(h.read_bytes())
    diff=''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),fromfile='a/SU2_CFD/src/solvers/CEulerSolver.cpp',tofile='b/SU2_CFD/src/solvers/CEulerSolver.cpp'))
    (OUT/'source_patch.diff').write_bytes(diff.encode())
    # 原v115的热边界及包含路径修复必须逐项与已交付运行时绑定一致。
    previous=json.loads((ROOT/'tools/runtime/su2/8.6.0-rah-axial-wall-v1/runtime.json').read_text(encoding='utf-8'))
    locked={str(src):sha(src),str(Path(__file__).resolve()):sha(Path(__file__)),str(hook):sha(hook)}
    for p,h in previous['buildSourceSHA256'].items():
        path=Path(p)
        if str(path).startswith(str(SOURCE)):
            assert sha(path)==h;locked[str(path)]=h
    for p in headers+[bindir/'g++.exe',ninja,ROOT/'tests/testSU2FrozenEnergyRecovery.cpp']:
        locked[str(p)]=sha(p)
    # 全部构建输入清单，防止仅记录补丁而遗漏被编译的上游文件。
    source_files={str(p):sha(p) for p in SOURCE.rglob('*') if p.is_file()}
    (OUT/'source_tree.json').write_text(json.dumps(source_files,indent=2)+'\n')
    locked[str(OUT/'source_tree.json')]=sha(OUT/'source_tree.json')
    (OUT/'build_sources.json').write_text(json.dumps(locked,indent=2)+'\n')
    builddir=BASE/'build-v117'
    build.run('configure',[sys.executable,str(BASE/'meson/meson.py'),'setup',str(builddir),str(SOURCE),
        '-Dwith-mpi=disabled','-Dwith-omp=true','-Denable-cgns=false','-Denable-tecio=false',
        '-Denable-pywrapper=false','-Denable-tests=false','-Dcpu-arch=x86-64','--buildtype=release'],SOURCE,env,120)
    build.run('compile',[str(ninja),'-C',str(builddir),'-j','2','SU2_CFD/src/SU2_CFD.exe'],SOURCE,env,720)
    assert all(sha(Path(p))==h for p,h in locked.items())
    assert all(sha(Path(p))==h for p,h in source_files.items())
    exe=builddir/'SU2_CFD/src/SU2_CFD.exe'
    imports=subprocess.check_output([str(bindir/'objdump.exe'),'-p',str(exe)]).decode()
    import re
    dlls=re.findall(r'DLL Name:\s*(\S+)',imports)
    assert set(map(str.lower,dlls))<=set(['kernel32.dll','msvcrt.dll'])
    dest=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR);assert not dest.exists();(dest/'bin').mkdir(parents=True)
    shutil.copyfile(str(exe),str(dest/'bin/SU2_CFD.exe'))
    old=ROOT/'tools/runtime/su2/8.6.0-rah-axial-wall-v1'
    for p in old.iterdir():
        if p.is_file() and p.name not in ['runtime.json','source_patch.diff']:shutil.copyfile(str(p),str(dest/p.name))
    shutil.copyfile(str(old/'source_patch.diff'),str(dest/'thermal_source_patch.diff'))
    shutil.copyfile(str(OUT/'source_patch.diff'),str(dest/'source_patch.diff'))
    state=dict(flavor=FLAVOR,version='8.6.0+'+FLAVOR,status='experimental-not-physics-qualified',
        executable=str(dest/'bin/SU2_CFD.exe'),archiveSHA256=sha(BASE/'su2.zip'),sourceCommit=previous['sourceCommit'],
        systemDLLs=dlls,buildSourceSHA256=locked,buildReceipt=json.loads((OUT/'compile_process.json').read_text()),
        compilerVersion=subprocess.check_output([str(bindir/'g++.exe'),'--version']).decode(),
        license=previous['license'],physicsQualified=False,
        restrictions=['explicit diagnostic/apply only','steady single-zone ideal-gas NS','SI dimensional','MGLEVEL=0',
                      'serial-MPI OpenMP','adiabatic or constant isothermal','experimental relaxation, no qualification'],
        files=[dict(path=p.relative_to(dest).as_posix(),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(dest.rglob('*')) if p.is_file()])
    (dest/'runtime.json').write_text(json.dumps(state,indent=2)+'\n')
    (OUT/'runtime_installed.json').write_bytes((dest/'runtime.json').read_bytes())
    print(json.dumps(dict(installed=str(dest),executableSHA256=sha(dest/'bin/SU2_CFD.exe'))),flush=True)
if __name__=='__main__':main()
