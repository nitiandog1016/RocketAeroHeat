"""模块用途。三对象有界构建：恢复v164流动对象，只改变SST共同松弛入口。"""
import difflib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from prepareSU2SSTEnergyPatch import ROOT,SOURCE,BUILD,PARENT,FILE,HEADER,FLAVOR,overlay
from buildSU2IdealAirFirstSample import arguments
from su2IdealAirInterface import SST_ENERGY
from su2ThermalRuntime import read,sha
from validateIdealAir import run,save

OUT=ROOT/'outputs/v167_validation/build01'
RESTORE=('SU2_CFD/src/variables/CNSVariable.cpp','SU2_CFD/src/solvers/CNSSolver.cpp')


def main():
    OUT.mkdir(parents=True,exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    fixture=read(ROOT/'outputs/v167_validation/fixture02/report.json')
    assert fixture['passed'] and all(sha(Path(n))==v for n,v in fixture['sourceSHA256'].items())
    before={str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}
    parent=read(PARENT/'objects_after.json')
    assert len(before)==298 and before==read(ROOT/'outputs/v165_validation/build01/objects_after.json')
    save(OUT/'objects_before.json',before)
    old=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-admissible-v1'
    previous=read(old/'runtime.json')
    for entry in previous['files']:assert sha(old/entry['path'])==entry['sha256']
    entries=read(BUILD/'compile_commands.json')
    selected=[next(row for row in entries if row['file'].replace('\\','/').endswith('/'+name)) for name in (*RESTORE,FILE)]
    objects={str(Path(row['output'])) for row in selected}
    assert objects.issubset(before) and len(objects)==3
    commands={'compile%d'%i:arguments(row['command']) for i,row in enumerate(selected)}
    commands['link']=read(PARENT/'commands.json')['link']
    save(OUT/'commands.json',commands)
    worker=[sys.executable,str(Path(__file__).resolve()),'--worker']
    with overlay(OUT) as snapshots:
        for i in range(3):run(OUT,'compile%d'%i,worker+['compile%d'%i],60)
        assert sum(read(OUT/('compile%d_process.json'%i))['elapsedSeconds'] for i in range(3))<=120
        unchanged={n:sha(BUILD/n)==v==parent[n] for n,v in before.items() if n not in objects}
        restored={str(Path(row['output'])):sha(BUILD/row['output'])==parent[str(Path(row['output']))] for row in selected[:2]}
        save(OUT/'unchanged_objects.json',unchanged);save(OUT/'restored_objects.json',restored)
        assert len(unchanged)==295 and all(unchanged.values()) and all(restored.values())
        # 既有归档为thin archive，直接引用对象；核对成员，不能只把新obj放在目录里。
        archive=BUILD/'SU2_CFD/src/libSU2core.a'
        assert archive.read_bytes()[:8]==b'!<thin>\n'
        ar=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/ar.exe'
        members=subprocess.check_output([str(ar),'t',str(archive)],cwd=str(BUILD)).decode().splitlines()
        assert all(any(Path(m).name==Path(row['output']).name for m in members) for row in selected)
        save(OUT/'archive_members.json',dict(thin=True,members=members))
        run(OUT,'link',worker+['link'],30)
        save(OUT/'objects_after.json',{n:sha(BUILD/n) for n in before})
    target=ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
    (target/'bin').mkdir(parents=True,exist_ok=False)
    shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'),str(target/'bin/SU2_CFD.exe'))
    for entry in previous['files']:
        if entry['path'].endswith(('.exe','.diff')):continue
        dest=target/entry['path'];dest.parent.mkdir(parents=True,exist_ok=True)
        os.link(str(old/entry['path']),str(dest))
    delta=''
    for name in (FILE,HEADER):
        prior=(SOURCE/name).read_text('utf8') if (SOURCE/name).exists() else ''
        delta+=''.join(difflib.unified_diff(prior.splitlines(True),snapshots[SOURCE/name].decode('utf8').splitlines(True),
                                         fromfile='a/'+name,tofile='b/'+name))
    (OUT/'sst_energy.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((old/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime=dict(previous)
    runtime.update(flavor=FLAVOR,version='8.6.0+'+FLAVOR,executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-tests-pending',parentRuntimeSHA256=sha(old/'runtime.json'),
        numericalBaselineFlavor=FLAVOR,numericalParentFlavor=previous['flavor'],
        samePhysicsAndNumerics=False,sstEnergy=dict(SST_ENERGY),
        incrementalObjectsRebuilt=3,inheritedVerifiedObjects=295,restoredParentObjects=2,allCFDObjectsRebuilt=0,
        sourceTreeSHA256=sha(OUT/'source_tree.json'),buildReceipt=read(OUT/'link_process.json'),
        buildSourceSHA256={str(p):sha(p) for p in (Path(__file__).resolve(),OUT/'executed_builder.py',
            ROOT/'scripts/prepareSU2SSTEnergyPatch.py',ROOT/'third_party/ideal_air/SSTEnergyUpdate.hpp',
            ROOT/'docs/validation/v167/PLAN.md',OUT/'commands.json',OUT/'objects_before.json',OUT/'source_tree.json')},
        files=[dict(path=p.relative_to(target).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    save(target/'runtime.json',runtime);save(OUT/'runtime_installed.json',runtime)
    print(json.dumps(dict(installed=str(target),objectsRebuilt=3)))


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--worker':
        raise SystemExit(subprocess.call(read(OUT/'commands.json')[sys.argv[2]],cwd=str(BUILD)))
    main()
