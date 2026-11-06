"""模块用途。继承v164准确对象，仅重编两个必要对象，发布独立限步数值包。"""
import difflib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from prepareSU2ConservativeStepPatch import ROOT,SOURCE,BUILD,PARENT,FILE,SOLVER,HEADER,FLAVOR,overlay
from buildSU2IdealAirFirstSample import arguments
from su2IdealAirInterface import CONSERVATIVE_STEP
from su2ThermalRuntime import read,sha
from validateIdealAir import run,save

OUT=ROOT/'outputs/v165_validation/build01'


def main():
    OUT.mkdir(parents=True,exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    assert read(ROOT/'outputs/v165_validation/fixture02/report.json')['passed']
    before={str(p.relative_to(BUILD)):sha(p) for p in BUILD.rglob('*.obj')}
    assert len(before)==298 and before==read(PARENT/'objects_after.json')
    save(OUT/'objects_before.json',before)
    old=ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-admissible-v1'
    previous=read(old/'runtime.json')
    for entry in previous['files']: assert sha(old/entry['path'])==entry['sha256']
    entries=read(BUILD/'compile_commands.json')
    selected=[next(row for row in entries if row['file'].replace('\\','/').endswith('/'+name)) for name in (FILE,SOLVER)]
    objects={str(Path(row['output'])) for row in selected}
    assert objects.issubset(before) and len(objects)==2
    commands={'compile%d'%i:arguments(row['command']) for i,row in enumerate(selected)}
    commands['link']=read(PARENT/'commands.json')['link']
    save(OUT/'commands.json',commands)
    worker=[sys.executable,str(Path(__file__).resolve()),'--worker']
    with overlay(OUT) as snapshots:
        for i in range(2):run(OUT,'compile%d'%i,worker+['compile%d'%i],60)
        assert sum(read(OUT/('compile%d_process.json'%i))['elapsedSeconds'] for i in range(2))<=120
        unchanged={n:sha(BUILD/n)==v for n,v in before.items() if n not in objects}
        assert len(unchanged)==296 and all(unchanged.values())
        save(OUT/'unchanged_objects.json',unchanged)
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
    for name in (FILE,SOLVER,HEADER):
        parent=PARENT/'executed_sources'/name
        prior=parent.read_text('utf8') if parent.exists() else ''
        delta+=''.join(difflib.unified_diff(prior.splitlines(True),snapshots[SOURCE/name].decode('utf8').splitlines(True),
                                         fromfile='a/'+name,tofile='b/'+name))
    (OUT/'conservative_step.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((old/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime=dict(previous);runtime.pop('admissibility')
    runtime.update(flavor=FLAVOR,version='8.6.0+'+FLAVOR,executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-tests-pending',parentRuntimeSHA256=sha(old/'runtime.json'),
        numericalBaselineFlavor=FLAVOR,numericalParentFlavor=previous['flavor'],
        samePhysicsAndNumerics=False,conservativeStep=dict(CONSERVATIVE_STEP),
        incrementalObjectsRebuilt=2,inheritedVerifiedObjects=296,allCFDObjectsRebuilt=0,
        sourceTreeSHA256=sha(OUT/'source_tree.json'),buildReceipt=read(OUT/'link_process.json'),
        buildSourceSHA256={str(p):sha(p) for p in (Path(__file__).resolve(),OUT/'executed_builder.py',
            ROOT/'scripts/prepareSU2ConservativeStepPatch.py',ROOT/'third_party/ideal_air/ConservativeStepTrace.hpp',
            ROOT/'third_party/ideal_air/ConservativeStep.inc',ROOT/'docs/validation/v165/PLAN.md',
            OUT/'commands.json',OUT/'objects_before.json',OUT/'source_tree.json')},
        files=[dict(path=p.relative_to(target).as_posix(),bytes=p.stat().st_size,sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    save(target/'runtime.json',runtime);save(OUT/'runtime_installed.json',runtime)
    print(json.dumps(dict(installed=str(target),objectsRebuilt=2)))


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--worker':
        raise SystemExit(subprocess.call(read(OUT/'commands.json')[sys.argv[2]],cwd=str(BUILD)))
    main()
