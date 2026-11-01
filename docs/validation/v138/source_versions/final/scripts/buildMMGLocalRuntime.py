"""固定MMG静态构建及局部入口；仅调整Windows构建脚本调用，不修改重网格内核。"""
from pathlib import Path
import difflib,json,os,time,tarfile
from prepareCGALRemeshingRuntime import ROOT,sha
from prepareMMGLocalRuntime import DEST,EVID
from buildCGALRemeshingRuntime import run

def main():
    start=time.monotonic();source=DEST/'source';build=DEST/'build'
    meta=json.loads((DEST/'source.json').read_text());assert meta['status']=='source-prepared'
    bindir=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin'
    # 准备失败发生在配置之前时，允许复原唯一已知构建补丁再按同源执行。
    with tarfile.open(str(DEST/'mmg-5.8.0.tar.gz')) as archive:
        original=archive.extractfile('mmg-5.8.0/CMakeLists.txt').read()
    f=source/'CMakeLists.txt';expected=original.decode().replace('COMMAND ./git_log_mmg.sh ',
        'COMMAND "'+(bindir/'sh.exe').as_posix()+'" ./git_log_mmg.sh ')
    assert f.read_text() in [original.decode(),expected]
    f.write_bytes(original)
    for n,r in meta['files'].items():assert sha(DEST/n)==r['sha256']
    cmake=ROOT/'tools/runtime/su2-coolprop-build/cmake/bin/cmake.exe';ninja=Path('D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe')
    os.environ['PATH']=str(bindir)+os.pathsep+os.environ['PATH']
    # 显式调用已存在sh；源tar内无.git，原脚本如实生成No git commit found。
    f=source/'CMakeLists.txt';old=f.read_text();needle='COMMAND ./git_log_mmg.sh '
    assert old.count(needle)==1
    new=old.replace(needle,'COMMAND "'+(bindir/'sh.exe').as_posix()+'" ./git_log_mmg.sh ')
    with f.open('w',newline='\n') as out:out.write(new)
    (EVID/'windows_build.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='upstream/CMakeLists.txt',tofile='local/CMakeLists.txt')))
    src=ROOT/'third_party/mmg/local_tetra_remesh.cpp'
    paths=list(source.rglob('*'))+[Path(__file__).resolve(),src,cmake,ninja,bindir/'gcc.exe',bindir/'g++.exe',bindir/'sh.exe',ROOT/'scripts/buildCGALRemeshingRuntime.py']
    locked={p.resolve().as_posix():sha(p) for p in paths if p.is_file()}
    configure=[str(cmake),'-S',str(source),'-B',str(build),'-G','Ninja','-DCMAKE_MAKE_PROGRAM='+str(ninja),
        '-DCMAKE_C_COMPILER='+str(bindir/'gcc.exe'),'-DCMAKE_CXX_COMPILER='+str(bindir/'g++.exe'),
        '-DCMAKE_BUILD_TYPE=Release','-DBUILD=MMG3D','-DLIBMMG3D_STATIC=ON','-DLIBMMG3D_SHARED=OFF',
        '-DBUILD_TESTING=OFF','-DUSE_SCOTCH=OFF','-DUSE_VTK=OFF','-DUSE_ELAS=OFF','-DCMAKE_C_FLAGS=-fcommon']
    receipts=[]
    for name,command,budget in [('configure',configure,60),('compile',[str(cmake),'--build',str(build),'--target','libmmg3d_a','--parallel','2'],240)]:
        r=run(command,EVID/name,budget,locked);receipts.append(r)
        print(name,r['exitCode'],r['seconds'],flush=True)
        assert r['exitCode']==0 and all(r['sourceStableAfter'].values()),name
    libs=list((build/'lib').glob('*.a'));assert len(libs)==1,libs
    exe=DEST/'local_tetra_remesh.exe'
    r=run([str(bindir/'g++.exe'),'-std=c++17','-O2','-static','-I'+str(build/'include'),str(src),str(libs[0]),'-o',str(exe)],EVID/'wrapper',60,locked)
    receipts.append(r);print('wrapper',r['exitCode'],r['seconds'],flush=True)
    assert r['exitCode']==0 and all(r['sourceStableAfter'].values())
    files={p.relative_to(DEST).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in DEST.rglob('*') if p.is_file()}
    record=dict(version='5.8.0',status='built-not-physics-qualified',executable=str(exe),executableSHA256=sha(exe),
        archiveSHA256=meta['archiveSHA256'],sourceSHA256=locked,processes=receipts,files=files,seconds=time.monotonic()-start,
        fileLengthBytes=sum(v['bytes'] for v in files.values()),physicsQualified=False)
    assert record['seconds']<360 and record['fileLengthBytes']<350*1024**2
    (DEST/'runtime.json').write_text(json.dumps(record,indent=2)+'\n');(EVID/'runtime.json').write_bytes((DEST/'runtime.json').read_bytes())
    print('built',record['seconds'],record['fileLengthBytes'],flush=True)
if __name__=='__main__':main()
