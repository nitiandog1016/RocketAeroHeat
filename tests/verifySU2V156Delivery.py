"""v156只读交付复核。重用实际流场和独立物性输出，核对最终源码与准确执行源。"""
from pathlib import Path
import copy
import csv
import json
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from validateIdealAir import save, sha, relative
from su2ThermalRuntime import read, resolve
from su2IdealAirCheckpoints import inventory
from prepareSU2IdealAirRestart import prepare


def main():
    output = ROOT/'outputs/v156_validation'
    out = output/'review01'
    out.mkdir(exist_ok=False)
    native = read(output/'native_report.json')
    assert native['passed'] and all(native['inputsStable'].values())
    assert read(output/'terminal01/report.json')['passed']
    summaries = {}
    max_energy = 0.
    for name in ('full', 'sparse', 'restart'):
        folder = output/'cases01'/name
        receipt = read(folder/'engine_process.json')
        assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0 and all(receipt['inputsStable'].values())
        for path, digest in receipt['outputSHA256'].items():
            assert sha(folder/path) == digest, path
        checkpoints = inventory(folder, receipt)
        if name == 'full':
            assert checkpoints['missingScheduledLocalIterations'] == [0]
        else:
            assert checkpoints['normalScheduleComplete']
        summaries[name] = dict(seconds=receipt['elapsedSeconds'], exitCode=0,
            savedLocalIterations=[entry['localIteration'] for entry in checkpoints['complete']],
            updates=checkpoints['seriesUpdatesAtLastFiniteHistory'],
            scheduledComplete=checkpoints['normalScheduleComplete'])
        audits = [folder/'ideal_air_audit'] if name != 'sparse' else sorted(folder.glob('audit_*'))
        for audit in audits:
            report = read(audit/'report.json')
            assert report['interfacePassed'] and not report['physicsQualified']
            for path, digest in report['inputSHA256'].items():
                # 原生计算及原始审核锁定时均稳定。已更新审核源码由精确快照承接。
                current = Path(path)
                if current.exists() and sha(current) == digest:
                    continue
                snapshots = [output/'executed_sources'/current.relative_to(ROOT),
                             output/'executed_sources02'/current.relative_to(ROOT)]
                assert any(p.exists() and sha(p) == digest for p in snapshots), path
            step = report['checkpointLocalIteration']
            with (folder/('restart_flow_%06d.csv' % step)).open() as stream:
                keys = next(csv.reader(stream))
            data = np.loadtxt(str(folder/('restart_flow_%06d.csv' % step)), delimiter=',', skiprows=1)
            fields = {key: data[:, index] for index, key in enumerate(keys)}
            oracle = np.genfromtxt(str(audit/'oracle.csv'), delimiter=',', names=True)
            velocity = np.column_stack([fields['Momentum_'+axis]/fields['Density'] for axis in ('x', 'y', 'z')])
            energy = fields['Energy']/fields['Density']-.5*np.sum(velocity**2, axis=1)-fields['Flow_Primitive_K']
            error = relative(energy, oracle['e'])
            max_energy = max(max_energy, error)
            assert error < 2e-10 and relative(fields['Pressure'], oracle['p']) < 2e-11
    # 从最终代码的新准备入口读取既有初场和实际源帧，原生运行预算为0。
    initial = out/'initial_contract'
    restarted = out/'restart_contract'
    prepare(ROOT/'outputs/v155_validation/cases01/small', initial, 6, 2, 32*1024**2,
            snapshot=ROOT/'outputs/v155_validation/executed_sources')
    prepare(output/'cases01/sparse', restarted, 3, 2, 32*1024**2, checkpoint=2)
    sys.path.insert(0, str(ROOT/'tests'))
    from testSU2IdealAirCheckpoints import reject_contracts
    # 独立最终合同负例使用新目录，不覆盖早先7项负例。
    import testSU2IdealAirCheckpoints
    testSU2IdealAirCheckpoints.OUT = out
    reject_contracts(initial)
    assert not (initial/'engine_process.json').exists() and not (restarted/'engine_process.json').exists()
    patch = read(output/'phase_patch/patch_manifest.json')
    for name, identity in patch.items():
        assert sha(output/'phase_patch'/name) == identity['patchedSHA256']
        assert (output/'phase_patch'/name).read_bytes() == (output/'build02/executed_sources'/name).read_bytes()
    runtime = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-phase-v2'
    manifest = read(runtime/'runtime.json')
    assert manifest['incrementalObjectsRebuilt'] == 1 and manifest['inheritedVerifiedObjects'] == 297
    for entry in manifest['files']:
        assert sha(runtime/entry['path']) == entry['sha256']
    objdump = ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/objdump.exe'
    imports = subprocess.check_output([str(objdump), '-p', manifest['executable']]).decode()
    dlls = re.findall(r'DLL Name:\s*(\S+)', imports)
    assert set(name.lower() for name in dlls) == {'kernel32.dll', 'msvcrt.dll'}
    for build in ('build01', 'build02'):
        assert all(read(output/build/'source_restore.json').values())
    save(out/'report.json', dict(passed=True, cases=summaries, independentPhaseEnergyError=max_energy,
        finalSourcePreflightPassed=True, finalPatchedFilesMatchExecuted=10, systemDLLs=dlls,
        additionalNativeRuns=0, physicsQualified=False, newQualifiedArea=0))
    print(json.dumps(dict(passed=True, cases=summaries, independentPhaseEnergyError=max_energy)))


if __name__ == '__main__':
    main()
