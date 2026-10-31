"""固定源码原生壁温适配构建；两线程/720秒，失败保留真实终态及可续构建状态。"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,time,zipfile,difflib

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'tools/runtime/su2-thermal-build'
SOURCE=BASE/'source';OUT=ROOT/'outputs/v115_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(name,command,cwd,env,seconds):
    receipt=OUT/(name+'_process.json');assert not receipt.exists(),'已有实际记录，不盲目重跑'
    state=dict(command=list(map(str,command)),cwd=str(cwd),status='starting',maximumSeconds=seconds)
    start=time.monotonic();timed=False
    with (OUT/(name+'.log')).open('wb') as log:
        proc=subprocess.Popen(command,cwd=str(cwd),env=env,stdout=log,stderr=subprocess.STDOUT,
                              creationflags=subprocess.CREATE_NO_WINDOW)
        state.update(status='running',processID=proc.pid)
        receipt.write_text(json.dumps(state,indent=2))
        try:code=proc.wait(timeout=seconds)
        except subprocess.TimeoutExpired:
            timed=True
            subprocess.run(['taskkill','/PID',str(proc.pid),'/T','/F'],stdout=log,stderr=log)
            code=proc.wait()
    state.update(status='terminal',exitCode=code,timedOut=timed,elapsedSeconds=time.monotonic()-start)
    receipt.write_text(json.dumps(state,indent=2));print(json.dumps(state),flush=True)
    assert code==0 and not timed,'构建步骤未完成，检查真实日志并保留对象，不冒充可用引擎'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--attempt',default='01');args=parser.parse_args()
    assert args.attempt.isdigit();OUT.mkdir(parents=True,exist_ok=True)
    source=SOURCE/'Common/src/geometry/CGeometry.cpp'
    header=ROOT/'third_party/su2_thermal/RAHAxialWallTemperature.hpp'
    body=ROOT/'third_party/su2_thermal/SetCustomBoundary.cpp.inc'
    with zipfile.ZipFile(str(BASE/'su2.zip')) as z:
        name=next(n for n in z.namelist() if n.endswith('/Common/src/geometry/CGeometry.cpp'))
        original=z.read(name).decode('utf-8')
        wrapper_name=next(n for n in z.namelist() if n.endswith('/SU2_CFD/src/python_wrapper_structure.cpp'))
        wrapper_original=z.read(wrapper_name).decode('utf-8')
    begin=original.index('void CGeometry::SetCustomBoundary(CConfig* config) {')
    end=original.index('\nvoid CGeometry::UpdateCustomBoundaryConditions',begin)
    patched=original[:begin]+body.read_text(encoding='utf-8')+'\n'+original[end:]
    patched=patched.replace('#include <unordered_set>','#include <unordered_set>\n#include "RAHAxialWallTemperature.hpp"',1)
    assert source.read_bytes() in (original.encode(),patched.encode()),'上游源存在未记录修改'
    source.write_bytes(patched.encode());(source.parent/header.name).write_bytes(header.read_bytes())
    diff=''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),
         fromfile='a/Common/src/geometry/CGeometry.cpp',tofile='b/Common/src/geometry/CGeometry.cpp'))
    (OUT/('source_patch'+args.attempt+'.diff')).write_bytes(diff.encode())
    # 无自定义标记时上游不会调用SetCustomBoundary，必须在其调用点拒绝遗漏。
    marker='void CDriver::PreprocessPythonInterface(CConfig** config, CGeometry**** geometry, CSolver***** solver) {'
    guard='''
  const char* rah_file = std::getenv("RAH_WALL_TEMPERATURE_FILE");
  if (rah_file && *rah_file && (nZone != 1 || config[0]->GetnMarker_PyCustom() != 1))
    SU2_MPI::Error("RAH wall table requires exactly one custom wall marker in one zone", CURRENT_FUNCTION);
'''
    wrapper_new=wrapper_original.replace(marker,marker+guard,1)
    wrapper_new='#include <cstdlib>\n'+wrapper_new
    wrapper=SOURCE/'SU2_CFD/src/python_wrapper_structure.cpp'
    assert wrapper.read_bytes() in (wrapper_original.encode(),wrapper_new.encode())
    wrapper.write_bytes(wrapper_new.encode())
    # 上游8.6.0此FSI文件的相对包含少一级；仅修复编译路径，不修改FSI算法。
    fsi=SOURCE/'SU2_CFD/src/interfaces/fsi/CFlowTractionInterface.cpp'
    with zipfile.ZipFile(str(BASE/'su2.zip')) as z:
        fsi_name=next(n for n in z.namelist() if n.endswith('/SU2_CFD/src/interfaces/fsi/CFlowTractionInterface.cpp'))
        fsi_original=z.read(fsi_name).decode('utf-8')
    fsi_new=fsi_original.replace('"../../../Common/include/interface_interpolation/CInterpolator.hpp"',
                                 '"../../../../Common/include/interface_interpolation/CInterpolator.hpp"',1)
    assert fsi_new != fsi_original
    assert fsi.read_bytes() in (fsi_original.encode(),fsi_new.encode())
    fsi.write_bytes(fsi_new.encode())
    diff+=''.join(difflib.unified_diff(fsi_original.splitlines(True),fsi_new.splitlines(True),
                 fromfile='a/SU2_CFD/src/interfaces/fsi/CFlowTractionInterface.cpp',
                 tofile='b/SU2_CFD/src/interfaces/fsi/CFlowTractionInterface.cpp'))
    include_fixes=[]
    with zipfile.ZipFile(str(BASE/'su2.zip')) as z:
        for path,old in [
            ('SU2_CFD/include/interfaces/cfd/CMixingPlaneInterface.hpp','../../../Common/include/containers/C2DContainer.hpp'),
            ('SU2_CFD/src/interfaces/cfd/CMixingPlaneInterface.cpp','../../../Common/include/interface_interpolation/CInterpolator.hpp')]:
            p=SOURCE/path;original_text=z.read(next(n for n in z.namelist() if n.endswith('/'+path))).decode('utf-8')
            fixed=original_text.replace('"'+old+'"','"../'+old+'"',1)
            assert fixed!=original_text and p.read_bytes() in (original_text.encode(),fixed.encode())
            assert (p.parent/('../'+old)).is_file()
            p.write_bytes(fixed.encode());include_fixes.append(p)
            diff+=''.join(difflib.unified_diff(original_text.splitlines(True),fixed.splitlines(True),
                                             fromfile='a/'+path,tofile='b/'+path))
    diff+=''.join(difflib.unified_diff(wrapper_original.splitlines(True),wrapper_new.splitlines(True),
                 fromfile='a/SU2_CFD/src/python_wrapper_structure.cpp',tofile='b/SU2_CFD/src/python_wrapper_structure.cpp'))
    (OUT/('source_patch'+args.attempt+'.diff')).write_bytes(diff.encode())
    for name in ('mel','eigen'):
        dst=SOURCE/'externals'/name
        for p in (BASE/name).rglob('*'):
            if not p.is_file() or p.name=='.rah_extracted.json':continue
            q=dst/p.relative_to(BASE/name);q.parent.mkdir(parents=True,exist_ok=True)
            if q.exists():assert q.read_bytes()==p.read_bytes()
            else:shutil.copyfile(str(p),str(q))
    binpath=BASE/'toolchain/w64devkit/bin'
    ninja=Path('D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe')
    env=dict(os.environ,PATH=str(binpath)+os.pathsep+str(ninja.parent)+os.pathsep+os.environ['PATH'],
             CC=str(binpath/'gcc.exe'),CXX=str(binpath/'g++.exe'))
    lock={str(p.resolve()):sha(p) for p in [Path(__file__),header,body,source,wrapper,fsi,binpath/'g++.exe',ninja]+include_fixes}
    (OUT/('build_sources'+args.attempt+'.json')).write_bytes((json.dumps(lock,indent=2)+'\n').encode())
    # 官方预配置明确禁用未需要的可选依赖下载；本次只编译正常CFD目标。
    if not (SOURCE/'su2preconfig.timestamp').exists():
        command=[sys.executable,'preconfigure.py','--no-codi','--no-medi','--no-opdi','--no-mpp',
                 '--no-coolprop','--no-mel','--no-fado','--no-mlpcpp','--no-eigen']
        run('preconfigure'+args.attempt,command,SOURCE,env,60)
    build=BASE/'build'
    if not (build/'build.ninja').exists():
        cmd=[sys.executable,str(BASE/'meson/meson.py'),'setup',str(build),str(SOURCE),
             '-Dwith-mpi=disabled','-Dwith-omp=true','-Denable-cgns=false','-Denable-tecio=false',
             '-Denable-pywrapper=false','-Denable-tests=false','-Dcpu-arch=x86-64',
             '--buildtype=release','--prefix='+str(BASE/'installed')]
        run('configure'+args.attempt,cmd,SOURCE,env,120)
    run('compile'+args.attempt,[str(ninja),'-C',str(build),'-j','2','SU2_CFD/src/SU2_CFD.exe'],SOURCE,env,720)
    assert all(sha(ROOT/p)==h for p,h in lock.items()),'构建中源码变化'
    print('Native thermal executable compiled; numerical and compatibility tests still required.',flush=True)

if __name__=='__main__':main()
