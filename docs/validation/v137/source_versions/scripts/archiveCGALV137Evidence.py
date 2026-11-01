"""复核本轮已有结果并绑定源码/许可证/存储；不再次执行重网格或流场。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from scipy.io import loadmat
from prepareCGALRemeshingRuntime import ROOT,DEST,EVID,sha

def main():
    out=ROOT/'docs/validation/v137';versions=out/'source_versions';versions.mkdir(exist_ok=True)
    test=json.loads((out/'local_tests/report.json').read_text());real=json.loads((out/'real_admissibility/report.json').read_text())
    assert test['passed'] and not real['permittedDirectImport'] and not real['realRemeshingExecuted']
    for r in [test,real]:
        locked=r.get('lockedSourceSHA256',r.get('sourceSHA256'))
        assert all(sha(ROOT/n)==h for n,h in locked.items())
    for item in test['cases']:
        if 'filesSHA256' in item:
            for n,h in item['filesSHA256'].items():assert sha(out/'local_tests'/item['case']/n)==h
    # 体积按有向边界积分再核对，与原四面体体积和凸包分别比较。
    from cgalTetrahedralContract import faces_and_neighbors
    d=loadmat(str(ROOT/'outputs/v136_validation/explicit/hybrid_volume.mat'));P=d['vertices'];T=d['tetrahedra'].astype(int)-1
    checks=[]
    for name,ids,expected in [('full',np.arange(len(T)),real['fullExterior'])]+[
        ('component'+str(r['component']),np.array(r['sourceTetRowsOneBased'])-1,r) for r in real['components']]:
        F,_=faces_and_neighbors(T[ids]);ref=P[np.unique(F)].mean(axis=0);q=P[F]-ref
        volume=float(np.sum(np.einsum('ij,ij->i',q[:,0],np.cross(q[:,1],q[:,2])))/6)
        error=abs(volume-expected['tetVolume'])/expected['tetVolume'];assert error<2e-10
        assert expected['convexHullVolume']-volume>2e-10*expected['convexHullVolume']
        checks.append(dict(scope=name,boundaryVolume=volume,cellVsBoundaryRelative=error))
    ga=json.loads((ROOT/'outputs/v136_validation/explicit/geometry_audit.json').read_text())
    assert abs(real['fullExterior']['convexHullVolume']/ga['outerVolume']-1)<2e-12
    sources=['scripts/prepareCGALRemeshingRuntime.py','scripts/buildCGALRemeshingRuntime.py',
        'scripts/cgalTetrahedralContract.py','scripts/reviewCGALRocketAdmissibility.py',
        'tests/testCGALLocalRemeshing.py','third_party/cgal/local_tetra_remesh.cpp',
        'scripts/reviewSU2ExplicitQuality.py','scripts/auditSU2HybridTopology.py','scripts/su2PyramidGeometry.py',
        'scripts/archiveCGALV137Evidence.py']
    for n in sources:
        target=versions/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(str(ROOT/n),str(target))
    licenses=out/'runtime/licenses';licenses.mkdir(exist_ok=True)
    for n in ['LICENSE','LICENSE.GPL','LICENSE.LGPL','LICENSE.BSL']:
        shutil.copyfile(str(DEST/'cgal'/n),str(licenses/('CGAL-'+n)))
    shutil.copyfile(str(DEST/'boost/LICENSE_1_0.txt'),str(licenses/'BOOST-LICENSE_1_0.txt'))
    for n in ['IO/File_medit.h','tetrahedral_remeshing.h','version.h']:
        shutil.copyfile(str(DEST/'cgal/include/CGAL'/n),str(out/'official'/('CGAL-6.0.1-'+Path(n).name)))
    runtime=json.loads((DEST/'runtime.json').read_text());exe=Path(runtime['executable']);assert sha(exe)==runtime['executableSHA256']
    dependencies=shutil.which('objdump') or str(ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/objdump.exe')
    import subprocess
    imports=subprocess.check_output([dependencies,'-p',str(exe)],creationflags=subprocess.CREATE_NO_WINDOW).decode()
    dlls=[line.strip().split(':',1)[1].strip() for line in imports.splitlines() if 'DLL Name:' in line]
    assert set(dlls)<=set(['KERNEL32.dll','msvcrt.dll','USER32.dll','ADVAPI32.dll','SHELL32.dll','WS2_32.dll'])
    failed=(ROOT/'scripts/reviewCGALRocketAdmissibility.py').read_text(encoding='utf8')
    failed=failed.replace("out.mkdir(parents=True,exist_ok=True);assert not (out/'report.json').exists()","out.mkdir(parents=True,exist_ok=False)")
    failed=failed.replace("source/'geometry_audit.json',Path(__file__).resolve(),","source/'geometry_audit.json',Path(__file__),")
    (versions/'review_before_relative_path_fix.py').write_text(failed,encoding='utf8')
    preparation=dict(status='failed-before-numeric-review',exitCode=1,
        error="ValueError: relative scripts/reviewCGALRocketAdmissibility.py not under absolute E:/RocketAeroHeat",
        corrected='resolve __file__; allow existing empty evidence directory, refuse report overwrite',
        historicalSource='source_versions/review_before_relative_path_fix.py',
        limitation='source reconstructed from the two recorded edits; no pre-failure SHA was captured',
        numericalCandidateExecuted=False)
    (out/'preparation_failure.json').write_text(json.dumps(preparation,indent=2)+'\n')
    runtimeFiles=[p for p in DEST.rglob('*') if p.is_file()];rtbytes=sum(p.stat().st_size for p in runtimeFiles)
    assert rtbytes<=650*1024**2
    report=dict(passed=True,independentVolumeChecks=checks,systemDLLs=dlls,
        nativeBuildExit=runtime['buildReceipt']['exitCode'],nativeBuildSeconds=runtime['buildReceipt']['seconds'],
        filesStable=True,physicsQualified=False,newQualifiedArea=0,
        runtimeFileCount=len(runtimeFiles),runtimeFileLengthBytes=rtbytes,runtimeUploaded=False,
        originalV136MATUnchanged=True,realRemeshingExecuted=False,noSU2OrMATLAB=True,
        deletedUniqueDataFiles=0,storageScope='file lengths, not physical allocated space; no whole outputs rescan',
        currentSourceSHA256={n:sha(ROOT/n) for n in sources})
    (out/'independent_review.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
