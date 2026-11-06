"""模块用途。继承v163准确对象，仅重编状态回退对象并发布独立数值身份。"""
import difflib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from prepareSU2AdmissiblePatch import ROOT, SOURCE, BUILD, PARENT, FILE, HEADER, FLAVOR, overlay
from buildSU2IdealAirFirstSample import arguments
from su2IdealAirInterface import ADMISSIBILITY
from su2ThermalRuntime import read, sha
from validateIdealAir import run, save

OUT = ROOT/'outputs/v164_validation/build01'


def main():
    """模块一。准确继承审核；共享源临时覆盖后必恢复，旧运行包保持原文件。"""
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    assert read(ROOT/'outputs/v164_validation/fixture02/report.json')['passed']
    before = {str(p.relative_to(BUILD)): sha(p) for p in BUILD.rglob('*.obj')}
    assert len(before) == 298 and before == read(PARENT/'objects_after.json')
    save(OUT/'objects_before.json', before)
    old = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-observed-v1'
    previous = read(old/'runtime.json')
    for entry in previous['files']:
        assert sha(old/entry['path']) == entry['sha256']
    entries = read(BUILD/'compile_commands.json')
    selected = next(row for row in entries if row['file'].replace('\\', '/').endswith('/'+FILE))
    obj = str(Path(selected['output']))
    assert obj in before
    commands = dict(compile=arguments(selected['command']), link=read(PARENT/'commands.json')['link'])
    save(OUT/'commands.json', commands)
    worker = [sys.executable, str(Path(__file__).resolve()), '--worker']
    with overlay(OUT) as snapshots:
        run(OUT, 'compile', worker+['compile'], 60)
        unchanged = {name: sha(BUILD/name) == digest for name, digest in before.items() if name != obj}
        assert len(unchanged) == 297 and all(unchanged.values())
        save(OUT/'unchanged_objects.json', unchanged)
        run(OUT, 'link', worker+['link'], 30)
        save(OUT/'objects_after.json', {name: sha(BUILD/name) for name in before})
    # 模块二。只新增具名包，运行库采用已有逐字文件；数值补丁与父补丁一并保留。
    target = ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
    (target/'bin').mkdir(parents=True, exist_ok=False)
    shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'), str(target/'bin/SU2_CFD.exe'))
    for entry in previous['files']:
        if entry['path'].endswith(('.exe', '.diff')):
            continue
        destination = target/entry['path']
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.link(str(old/entry['path']), str(destination))
    delta = ''
    for name in (FILE, HEADER):
        parent = PARENT/'executed_sources'/name
        old_text = parent.read_text('utf8') if parent.exists() else ''
        delta += ''.join(difflib.unified_diff(old_text.splitlines(True),
            snapshots[SOURCE/name].decode('utf8').splitlines(True), fromfile='a/'+name, tofile='b/'+name))
    (OUT/'admissibility.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((old/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime = dict(previous)
    # 本轮修改实际候选态处理，不能继续沿用“数值相同”的父包结论。
    runtime.update(flavor=FLAVOR, version='8.6.0+'+FLAVOR, executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-tests-pending', parentRuntimeSHA256=sha(old/'runtime.json'),
        numericalBaselineFlavor=FLAVOR, numericalParentFlavor=previous['flavor'],
        samePhysicsAndNumerics=False, admissibility=dict(ADMISSIBILITY),
        incrementalObjectsRebuilt=1, inheritedVerifiedObjects=297, allCFDObjectsRebuilt=0,
        sourceTreeSHA256=sha(OUT/'source_tree.json'), buildReceipt=read(OUT/'link_process.json'),
        buildSourceSHA256={str(p): sha(p) for p in (Path(__file__).resolve(), OUT/'executed_builder.py',
            ROOT/'scripts/prepareSU2AdmissiblePatch.py', ROOT/'third_party/ideal_air/AdmissibilityTrace.hpp',
            ROOT/'docs/validation/v164/PLAN.md', OUT/'commands.json', OUT/'objects_before.json', OUT/'source_tree.json')},
        files=[dict(path=p.relative_to(target).as_posix(), bytes=p.stat().st_size, sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    save(target/'runtime.json', runtime)
    save(OUT/'runtime_installed.json', runtime)
    print(json.dumps(dict(installed=str(target), objectsRebuilt=1)))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        raise SystemExit(subprocess.call(read(OUT/'commands.json')[sys.argv[2]], cwd=str(BUILD)))
    main()
