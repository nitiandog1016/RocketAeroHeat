"""只纠正实际流动writer的首帧条件，核对既有对象后单文件编译并恢复共享源。"""
from pathlib import Path
import difflib
import json
import os
import shlex
import shutil
import subprocess
import sys

from buildSU2IdealAirRuntime import ROOT, SOURCE, BUILD, BASE, sha, save
from prepareSU2IdealAirPatch import replace
from validateIdealAir import run

OUT = ROOT/'outputs/v156_validation/build02'
PREVIOUS = ROOT/'outputs/v156_validation/build01'
FLAVOR = 'rah-ideal-air-phase-v2'


def arguments(command):
    return [item[1:-1] if item.startswith('"') and item.endswith('"') else item
            for item in shlex.split(command, posix=False)]


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    assert json.loads((PREVIOUS/'compile_process.json').read_text())['exitCode'] == 0
    old = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-phase-v1'
    previous = json.loads((old/'runtime.json').read_text())
    for entry in previous['files']:
        assert sha(old/entry['path']) == entry['sha256']
    assert sha(BUILD/'SU2_CFD/src/SU2_CFD.exe') == sha(old/'bin/SU2_CFD.exe')
    snapshots = {SOURCE/path.relative_to(PREVIOUS/'executed_sources'): path.read_bytes()
                 for path in (PREVIOUS/'executed_sources').rglob('*') if path.is_file()}
    changed = SOURCE/'SU2_CFD/src/output/CFlowOutput.cpp'
    original = changed.read_bytes()
    text = replace(original.decode('utf8'),
        '    return ((Iter > 0) && Iter % config->GetVolumeOutputFrequency(iFile) == 0) || force_writing;',
        '    // 新变比热包保存第一次完成更新，稳态局部0对应第1次更新。\n'
        '    return ((Iter > 0 || config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP) &&\n'
        '            Iter % config->GetVolumeOutputFrequency(iFile) == 0) || force_writing;')
    snapshots[changed] = text.encode('utf8')
    temporary = {path: path.read_bytes() if path.exists() else None for path in snapshots}
    before = {str(path.relative_to(BUILD)): sha(path) for path in BUILD.rglob('*.obj')}
    save(OUT/'objects_before.json', before)
    try:
        for path, raw in snapshots.items():
            path.write_bytes(raw)
        old_tree = json.loads((PREVIOUS/'source_tree.json').read_text())
        assert all(sha(Path(name)) == digest for name, digest in old_tree.items() if Path(name) != changed)
        save(OUT/'source_tree.json', {name: sha(Path(name)) for name in old_tree})
        compile_entry = next(row for row in json.loads((BUILD/'compile_commands.json').read_text())
            if row['file'].replace('\\', '/').endswith('/src/output/CFlowOutput.cpp'))
        ninja = 'D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe'
        link = subprocess.check_output([ninja, '-C', str(BUILD), '-t', 'commands', 'SU2_CFD/src/SU2_CFD.exe']).decode().splitlines()[-1]
        save(OUT/'commands.json', dict(compile=arguments(compile_entry['command']), link=arguments(link)))
        # 父监督器继续管理真实子进程树，工作目录仅在固定构建目录内切换。
        for action in ('compile', 'link'):
            run(OUT, action, [sys.executable, str(Path(__file__).resolve()), '--worker', action], 60)
        object_name = compile_entry['output'].replace('\\', '/')
        unchanged = {name: sha(BUILD/name) == digest for name, digest in before.items() if name.replace('\\', '/') != object_name}
        assert len(unchanged) == 297 and all(unchanged.values())
        save(OUT/'unchanged_objects.json', unchanged)
        assert sha(BUILD/object_name) != before[object_name]
        target = ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
        (target/'bin').mkdir(parents=True, exist_ok=False)
        shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'), str(target/'bin/SU2_CFD.exe'))
        for entry in previous['files']:
            if entry['path'].endswith(('.exe', '.diff')):
                continue
            destination = target/entry['path']
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.link(str(old/entry['path']), str(destination))
        delta = ''.join(difflib.unified_diff(original.decode().splitlines(True), text.splitlines(True),
            fromfile='a/SU2_CFD/src/output/CFlowOutput.cpp', tofile='b/SU2_CFD/src/output/CFlowOutput.cpp'))
        (target/'ideal_air.diff').write_bytes((old/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
        (OUT/'writer.diff').write_bytes(delta.encode('utf8'))
        runtime = dict(previous)
        runtime.update(flavor=FLAVOR, version='8.6.0+'+FLAVOR, executable=str(target/'bin/SU2_CFD.exe'),
            parentRuntimeSHA256=sha(old/'runtime.json'), incrementalObjectsRebuilt=1,
            allCFDObjectsRebuilt=0, inheritedVerifiedObjects=297,
            sourceTreeSHA256=sha(OUT/'source_tree.json'), buildReceipt=json.loads((OUT/'link_process.json').read_text()),
            buildSourceSHA256={str(path): sha(path) for path in (Path(__file__).resolve(),
                ROOT/'docs/validation/v156/RECOVERY_PLAN.md', OUT/'commands.json', OUT/'objects_before.json')},
            files=[dict(path=path.relative_to(target).as_posix(), bytes=path.stat().st_size, sha256=sha(path))
                   for path in sorted(target.rglob('*')) if path.is_file()])
        save(target/'runtime.json', runtime)
        save(OUT/'runtime_installed.json', runtime)
        assert all(sha(Path(name)) == digest for name, digest in json.loads((OUT/'source_tree.json').read_text()).items())
        print(json.dumps(dict(installed=str(target), changedObjects=1, inheritedObjects=297)))
    finally:
        restored = {}
        for path, raw in temporary.items():
            target = OUT/'executed_sources'/path.relative_to(SOURCE)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(path.read_bytes())
            if raw is None:
                path.resolve().relative_to(SOURCE.resolve())
                path.unlink()
                restored[str(path)] = not path.exists()
            else:
                path.write_bytes(raw)
                restored[str(path)] = path.read_bytes() == raw
        save(OUT/'source_restore.json', restored)
        assert all(restored.values())


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        command = json.loads((OUT/'commands.json').read_text())[sys.argv[2]]
        raise SystemExit(subprocess.call(command, cwd=str(BUILD)))
    main()
