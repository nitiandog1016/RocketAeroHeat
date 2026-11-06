"""模块用途。恢复准确三维基线，增加终止记录和高精度输出，限三个对象并恢复共享源。"""
import difflib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from buildSU2IdealAirRuntime import ROOT, SOURCE, BUILD, BASE, sha, save
from buildSU2IdealAirFirstSample import arguments
from prepareSU2IdealAirPatch import replace
from validateIdealAir import run

OUT = ROOT/'outputs/v163_validation/build01'
PARENT = ROOT/'outputs/v156_validation/build02'
FLAVOR = 'rah-ideal-air-observed-v1'
FILES = ('SU2_CFD/src/solvers/CEulerSolver.cpp', 'SU2_CFD/src/SU2_CFD.cpp',
         'SU2_CFD/src/output/filewriter/CSU2FileWriter.cpp')


def main():
    """先核对继承对象和源码，再受控编译，发布带独立身份的运行包。"""
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    assert json.loads((ROOT/'outputs/v163_validation/fixture01/report.json').read_text())['passed']
    baseline = json.loads((ROOT/'outputs/v159_validation/build01/objects_before.json').read_text())
    before = {str(p.relative_to(BUILD)): sha(p) for p in BUILD.rglob('*.obj')}
    euler_object = str(Path('SU2_CFD/src/libSU2core.a.p/solvers_CEulerSolver.cpp.obj'))
    assert len(before) == len(baseline) == 298 and set(before) == set(baseline)
    assert [name for name in before if before[name] != baseline[name]] == [euler_object]
    save(OUT/'objects_before.json', before)
    old = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-phase-v2'
    previous = json.loads((old/'runtime.json').read_text())
    for item in previous['files']:
        assert sha(old/item['path']) == item['sha256']
    old_tree = json.loads((PARENT/'source_tree.json').read_text())
    snapshots = {SOURCE/p.relative_to(PARENT/'executed_sources'): p.read_bytes()
                 for p in (PARENT/'executed_sources').rglob('*') if p.is_file()}
    for name in FILES:
        path = SOURCE/name
        if path not in snapshots:
            snapshots[path] = path.read_bytes()
        assert sha_bytes(snapshots[path]) == old_tree[str(path)]
    inherited = dict(snapshots)
    main_file = SOURCE/FILES[1]
    main_text = replace(snapshots[main_file].decode('utf8'), '#include "../include/SU2_CFD.hpp"',
        '#include "../include/SU2_CFD.hpp"\n#include "../include/NativeTermination.hpp"')
    main_text = replace(main_text, 'int main(int argc, char *argv[]) {',
        'int main(int argc, char *argv[]) {\n  // 仅观察异常并保留原终止行为。\n  rah_termination::install();')
    snapshots[main_file] = main_text.encode('utf8')
    writer = SOURCE/FILES[2]
    text = replace(snapshots[writer].decode('utf8'),
        '#include "../../../include/output/filewriter/CSU2FileWriter.hpp"',
        '#include "../../../include/output/filewriter/CSU2FileWriter.hpp"\n#include <limits>')
    text = replace(text, '  data.precision(15);',
        '  // scientific的precision统计小数位，17位有效数字采用max_digits10减一。\n'
        '  data.precision(std::numeric_limits<double>::max_digits10 - 1);')
    snapshots[writer] = text.encode('utf8')
    header = SOURCE/'SU2_CFD/include/NativeTermination.hpp'
    assert not header.exists()
    snapshots[header] = (ROOT/'third_party/su2_thermal/NativeTermination.hpp').read_bytes()
    temporary = {path: path.read_bytes() if path.exists() else None for path in snapshots}
    entries = json.loads((BUILD/'compile_commands.json').read_text())
    selected = [next(row for row in entries if row['file'].replace('\\', '/').endswith('/'+name)) for name in FILES]
    commands = {'compile'+str(i): arguments(row['command']) for i, row in enumerate(selected)}
    ninja = 'D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe'
    link = subprocess.check_output([ninja, '-C', str(BUILD), '-t', 'commands',
                                   'SU2_CFD/src/SU2_CFD.exe']).decode().splitlines()[-1]
    commands['link'] = arguments(link)
    save(OUT/'commands.json', commands)
    try:
        for path, raw in snapshots.items():
            path.write_bytes(raw)
        changed_sources = {str(main_file), str(writer)}
        assert all(sha(Path(name)) == digest for name, digest in old_tree.items() if name not in changed_sources)
        tree = {name: sha(Path(name)) for name in old_tree}
        tree[str(header)] = sha(header)
        save(OUT/'source_tree.json', tree)
        started = time.monotonic()
        for i in range(3):
            remaining = 120-(time.monotonic()-started)
            assert remaining > 0
            run(OUT, 'compile'+str(i), [sys.executable, str(Path(__file__).resolve()), '--worker', 'compile'+str(i)], remaining)
        assert sha(BUILD/euler_object) == baseline[euler_object], '三维基线对象未逐字恢复'
        rebuilt = {str(Path(row['output'])) for row in selected}
        unchanged = {name: sha(BUILD/name) == digest for name, digest in before.items() if name not in rebuilt}
        assert len(unchanged) == 295 and all(unchanged.values())
        save(OUT/'unchanged_objects.json', unchanged)
        run(OUT, 'link', [sys.executable, str(Path(__file__).resolve()), '--worker', 'link'], 30)
        save(OUT/'objects_after.json', {name: sha(BUILD/name) for name in before})
        assert all(sha(Path(name)) == digest for name, digest in tree.items())
    finally:
        restored = {}
        for path, raw in temporary.items():
            snapshot = OUT/'executed_sources'/path.relative_to(SOURCE)
            snapshot.parent.mkdir(parents=True, exist_ok=True)
            snapshot.write_bytes(path.read_bytes())
            if raw is None:
                path.resolve().relative_to(SOURCE.resolve())
                path.unlink()
                restored[str(path)] = not path.exists()
            else:
                path.write_bytes(raw)
                restored[str(path)] = path.read_bytes() == raw
        save(OUT/'source_restore.json', restored)
        assert all(restored.values())
    target = ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
    (target/'bin').mkdir(parents=True, exist_ok=False)
    shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'), str(target/'bin/SU2_CFD.exe'))
    for item in previous['files']:
        if item['path'].endswith(('.exe', '.diff')):
            continue
        destination = target/item['path']
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.link(str(old/item['path']), str(destination))
    delta = ''
    for path in (main_file, writer, header):
        name = path.relative_to(SOURCE).as_posix()
        delta += ''.join(difflib.unified_diff(inherited.get(path, b'').decode('utf8').splitlines(True),
            snapshots[path].decode('utf8').splitlines(True), fromfile='a/'+name, tofile='b/'+name))
    (OUT/'observation.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((old/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime = dict(previous)
    runtime.update(flavor=FLAVOR, version='8.6.0+'+FLAVOR, executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-tests-pending', parentRuntimeSHA256=sha(old/'runtime.json'),
        numericalBaselineFlavor=previous['flavor'], asciiSignificantDigits=17,
        terminationObservation='cpp-current-exception-windows-stack-delegate-v1',
        incrementalObjectsRebuilt=3, inheritedVerifiedObjects=295, allCFDObjectsRebuilt=0,
        restoredEulerObjectSHA256=sha(BUILD/euler_object), sourceTreeSHA256=sha(OUT/'source_tree.json'),
        buildReceipt=json.loads((OUT/'link_process.json').read_text()),
        buildSourceSHA256={str(p): sha(p) for p in (Path(__file__).resolve(), OUT/'executed_builder.py',
            ROOT/'docs/validation/v163/PLAN.md', ROOT/'third_party/su2_thermal/NativeTermination.hpp',
            OUT/'commands.json', OUT/'objects_before.json', OUT/'source_tree.json')},
        files=[dict(path=p.relative_to(target).as_posix(), bytes=p.stat().st_size, sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    save(target/'runtime.json', runtime)
    save(OUT/'runtime_installed.json', runtime)
    print(json.dumps(dict(installed=str(target), objectsRebuilt=3)))


def sha_bytes(raw):
    import hashlib
    return hashlib.sha256(raw).hexdigest()


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        command = json.loads((OUT/'commands.json').read_text())[sys.argv[2]]
        raise SystemExit(subprocess.call(command, cwd=str(BUILD)))
    main()
