"""补齐变比热全状态入口和外推出口，只编译一个翻译单元并恢复共享源。"""
import difflib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from buildSU2IdealAirRuntime import ROOT, SOURCE, BUILD, sha, save
from buildSU2IdealAirFirstSample import arguments
from prepareSU2IdealAirPatch import replace
from validateIdealAir import run

OUT = ROOT/'outputs/v159_validation/build01'
PARENT = ROOT/'outputs/v156_validation/build02'
FLAVOR = 'rah-ideal-air-axisymmetric-v1'
FILE = 'SU2_CFD/src/solvers/CEulerSolver.cpp'


def patch(text):
    """模块一。修改限定于两处已有边界函数和显式新物性分支。"""
    start = text.index('void CEulerSolver::BC_Supersonic_Inlet(')
    end = text.index('\nvoid CEulerSolver::', start+1)
    inlet = text[start:end]
    anchor = '    su2double Energy = Pressure / (Density * Gamma_Minus_One) + 0.5 * Velocity2;'
    inlet = replace(inlet, anchor, anchor+'''
    // 变比热入口以同一PT状态构造热内能及一般EOS导数，SST k由下方原逻辑加入。
    const bool variable_air = config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP;
    su2double Secondary_inlet[2] = {0.0, 0.0};
    if (variable_air) {
      GetFluidModel()->SetTDState_PT(Pressure, Temperature);
      Energy = GetFluidModel()->GetStaticEnergy() + 0.5 * Velocity2;
      Secondary_inlet[0] = GetFluidModel()->GetdPdrho_e();
      Secondary_inlet[1] = GetFluidModel()->GetdPde_rho();
    }''')
    anchor = '    conv_numerics->SetPrimitive(V_domain, V_inlet);'
    inlet = replace(inlet, anchor, anchor+'''
    if (variable_air)
      conv_numerics->SetSecondary(nodes->GetSecondary(iPoint), Secondary_inlet);''')
    text = text[:start]+inlet+text[end:]
    start = text.index('void CEulerSolver::BC_Supersonic_Outlet(')
    end = text.index('\nvoid CEulerSolver::', start+1)
    outlet = text[start:end]
    anchor = '      conv_numerics->SetPrimitive(V_domain, V_outlet);'
    outlet = replace(outlet, anchor, anchor+'''
      // 外推边界两侧使用当前内部状态的一般EOS导数。
      if (config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP)
        conv_numerics->SetSecondary(nodes->GetSecondary(iPoint), nodes->GetSecondary(iPoint));''')
    anchor = '       visc_numerics->SetPrimitive(V_domain, V_domain);'
    outlet = replace(outlet, anchor, anchor+'''
       if (config->GetKind_FluidModel() == IDEAL_AIR_VARIABLE_CP)
         visc_numerics->SetSecondary(nodes->GetSecondary(iPoint), nodes->GetSecondary(iPoint));''')
    return text[:start]+outlet+text[end:]


def main():
    """模块二。核对父运行包和准确源树，增量编译/链接，独立发布及恢复。"""
    OUT.mkdir(parents=True, exist_ok=False)
    (OUT/'executed_builder.py').write_bytes(Path(__file__).read_bytes())
    old = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-phase-v2'
    previous = json.loads((old/'runtime.json').read_text())
    for item in previous['files']:
        assert sha(old/item['path']) == item['sha256']
    assert sha(BUILD/'SU2_CFD/src/SU2_CFD.exe') == sha(old/'bin/SU2_CFD.exe')
    old_tree = json.loads((PARENT/'source_tree.json').read_text())
    snapshots = {SOURCE/p.relative_to(PARENT/'executed_sources'): p.read_bytes()
                 for p in (PARENT/'executed_sources').rglob('*') if p.is_file()}
    changed = SOURCE/FILE
    inherited = snapshots[changed]
    snapshots[changed] = patch(inherited.decode('utf8')).encode('utf8')
    temporary = {p: p.read_bytes() if p.exists() else None for p in snapshots}
    before = {str(p.relative_to(BUILD)): sha(p) for p in BUILD.rglob('*.obj')}
    assert len(before) == 298
    save(OUT/'objects_before.json', before)
    try:
        for p, raw in snapshots.items():
            p.write_bytes(raw)
        assert all(sha(Path(name)) == digest for name, digest in old_tree.items() if Path(name) != changed)
        save(OUT/'source_tree.json', {name: sha(Path(name)) for name in old_tree})
        entry = next(row for row in json.loads((BUILD/'compile_commands.json').read_text())
                     if row['file'].replace('\\', '/').endswith('/src/solvers/CEulerSolver.cpp'))
        ninja = 'D:/matlab/toolbox/shared/coder/ninja/win64/ninja.exe'
        link = subprocess.check_output([ninja, '-C', str(BUILD), '-t', 'commands',
                                        'SU2_CFD/src/SU2_CFD.exe']).decode().splitlines()[-1]
        save(OUT/'commands.json', dict(compile=arguments(entry['command']), link=arguments(link)))
        run(OUT, 'compile', [sys.executable, str(Path(__file__).resolve()), '--worker', 'compile'], 150)
        run(OUT, 'link', [sys.executable, str(Path(__file__).resolve()), '--worker', 'link'], 30)
        object_name = str(Path(entry['output']))
        unchanged = {name: sha(BUILD/name) == digest for name, digest in before.items() if name != object_name}
        assert len(unchanged) == 297 and all(unchanged.values())
        assert sha(BUILD/object_name) != before[object_name]
        save(OUT/'unchanged_objects.json', unchanged)
        assert all(sha(Path(name)) == digest for name, digest in json.loads((OUT/'source_tree.json').read_text()).items())
    finally:
        restored = {}
        for p, raw in temporary.items():
            snapshot = OUT/'executed_sources'/p.relative_to(SOURCE)
            snapshot.parent.mkdir(parents=True, exist_ok=True)
            snapshot.write_bytes(p.read_bytes())
            if raw is None:
                p.resolve().relative_to(SOURCE.resolve())
                p.unlink()
                restored[str(p)] = not p.exists()
            else:
                p.write_bytes(raw)
                restored[str(p)] = p.read_bytes() == raw
        save(OUT/'source_restore.json', restored)
        assert all(restored.values())
    target = ROOT/'tools/runtime/su2'/('8.6.0-'+FLAVOR)
    (target/'bin').mkdir(parents=True, exist_ok=False)
    shutil.copyfile(str(BUILD/'SU2_CFD/src/SU2_CFD.exe'), str(target/'bin/SU2_CFD.exe'))
    for item in previous['files']:
        if item['path'].endswith(('.exe', '.diff')):
            continue
        dest = target/item['path']
        dest.parent.mkdir(parents=True, exist_ok=True)
        os.link(str(old/item['path']), str(dest))
    delta = ''.join(difflib.unified_diff(inherited.decode('utf8').splitlines(True),
        snapshots[changed].decode('utf8').splitlines(True), fromfile='a/'+FILE, tofile='b/'+FILE))
    (OUT/'boundary.diff').write_bytes(delta.encode('utf8'))
    (target/'ideal_air.diff').write_bytes((old/'ideal_air.diff').read_bytes()+delta.encode('utf8'))
    runtime = dict(previous)
    runtime.update(flavor=FLAVOR, version='8.6.0+'+FLAVOR, executable=str(target/'bin/SU2_CFD.exe'),
        status='compiled-tests-pending', parentRuntimeSHA256=sha(old/'runtime.json'),
        incrementalObjectsRebuilt=1, inheritedVerifiedObjects=297, allCFDObjectsRebuilt=0,
        sourceTreeSHA256=sha(OUT/'source_tree.json'), buildReceipt=json.loads((OUT/'link_process.json').read_text()),
        buildSourceSHA256={str(p): sha(p) for p in (Path(__file__).resolve(), OUT/'executed_builder.py',
            ROOT/'docs/validation/v159/PLAN.md', OUT/'commands.json', OUT/'objects_before.json')},
        files=[dict(path=p.relative_to(target).as_posix(), bytes=p.stat().st_size, sha256=sha(p))
               for p in sorted(target.rglob('*')) if p.is_file()])
    save(target/'runtime.json', runtime)
    save(OUT/'runtime_installed.json', runtime)
    print(json.dumps(dict(installed=str(target), objectsRebuilt=1)))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        command = json.loads((OUT/'commands.json').read_text())[sys.argv[2]]
        raise SystemExit(subprocess.call(command, cwd=str(BUILD)))
    main()
