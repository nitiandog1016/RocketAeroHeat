"""安装已真实构建的独立壁温运行时；检查DLL、许可证、构建源与成品身份。"""
from pathlib import Path
import argparse,hashlib,json,re,shutil,subprocess
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'tools/runtime/su2-thermal-build'
OUT=ROOT/'outputs/v115_validation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--attempt',required=True);args=p.parse_args()
    receipt=json.loads((OUT/('compile'+args.attempt+'_process.json')).read_text())
    assert receipt['status']=='terminal' and receipt['exitCode']==0 and not receipt['timedOut']
    locked=json.loads((OUT/('build_sources'+args.attempt+'.json')).read_text())
    assert all(sha(Path(p))==h for p,h in locked.items())
    exe=BASE/'build/SU2_CFD/src/SU2_CFD.exe';binpath=BASE/'toolchain/w64devkit/bin'
    imports=subprocess.check_output([str(binpath/'objdump.exe'),'-p',str(exe)]).decode()
    dlls=re.findall(r'DLL Name:\s*(\S+)',imports)
    assert set(map(str.lower,dlls))<=set(['kernel32.dll','msvcrt.dll']),dlls
    dest=ROOT/'tools/runtime/su2/8.6.0-rah-axial-wall-v1';assert not dest.exists()
    (dest/'bin').mkdir(parents=True)
    shutil.copyfile(str(exe),str(dest/'bin/SU2_CFD.exe'))
    for name in ['COPYING','LICENSE.md','AUTHORS.md']:
        shutil.copyfile(str(BASE/'source'/name),str(dest/name))
    shutil.copyfile(str(BASE/'toolchain/w64devkit/COPYING.MinGW-w64-runtime.txt'),str(dest/'COPYING.MinGW-w64-runtime.txt'))
    for name in ['COPYING.MPL2','COPYING.README']:
        shutil.copyfile(str(BASE/'source/externals/eigen'/name),str(dest/('Eigen_'+name)))
    for name in ['mel','CLI11']:
        shutil.copyfile(str(BASE/'source/externals'/name/'LICENSE'),str(dest/(name+'_LICENSE')))
    # 静态GCC运行库的许可例外随源码获取位置一并保留；运行时只在本地，不分发二进制。
    shutil.copyfile(str(OUT/('source_patch'+args.attempt+'.diff')),str(dest/'source_patch.diff'))
    manifest=dict(flavor='rah-axial-wall-v1',version='8.6.0+rah-axial-wall-v1',
        status='compiled-not-physics-qualified',executable=str(dest/'bin/SU2_CFD.exe'),
        archiveSHA256=sha(BASE/'su2.zip'),sourceCommit='7048fd35f56d4bf52e586429d6746a36a10255f8',
        systemDLLs=dlls,buildReceipt=receipt,buildSourceSHA256=locked,
        compilerVersion=subprocess.check_output([str(binpath/'g++.exe'),'--version']).decode(),
        compilerArchiveSHA256=sha(BASE/'compiler.7z.exe'),
        mesonNinjaSHA256=sha(binpath/'ninja.exe'),
        license='LGPL-2.1-or-later; local binary, source changes and rebuild recipe published',
        restrictions=['SI','one-zone','fixed-mesh','MGLEVEL=0','serial-MPI','custom isothermal wall'],
        physicsQualified=False,files=[dict(path=p.relative_to(dest).as_posix(),sha256=sha(p),bytes=p.stat().st_size)
        for p in sorted(dest.rglob('*')) if p.is_file()])
    assert sum(q['bytes'] for q in manifest['files'])<150*1024**2
    (dest/'runtime.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    (OUT/'runtime_installed.json').write_bytes((dest/'runtime.json').read_bytes())
    print(json.dumps(dict(installed=str(dest),bytes=sum(q['bytes'] for q in manifest['files']),dlls=dlls)))
if __name__=='__main__':main()
