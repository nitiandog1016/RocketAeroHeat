"""按v155事前计划验证正常MATLAB导出、新初场和唯一小型三维原生读取。"""
from pathlib import Path
import csv
import json
import os
import shutil
import sys

import numpy as np

from validateIdealAir import ROOT, run, save, sha, write_rows, read_rows, relative
from su2ThermalRuntime import resolve, read
from su2SSTInterface import state
from su2IdealAirInterface import R, settings
from prepareSU2SSTCase import set_value

OUT = ROOT/'outputs/v155_validation'


def independent_initial(folder):
    """模块一。逐节点检查写出量，独立CoolProp抽样覆盖壁、近壁、远场和内部状态。"""
    init = read(folder/'ideal_air_initialization.json')
    manifest = read(folder/'case_manifest.json')
    data, keys, dimension, rho, velocity, energy = state(folder/'initial_flow.csv', True)
    assert dimension == 3
    with np.load(str(folder/'ideal_air_wall_distance.npz')) as saved:
        detail = dict(saved)
    eta = detail['distance']/init['thicknessM']
    blend = eta**2*(3-2*eta)
    u_inf = np.asarray(manifest['expectedVelocity'])
    model = manifest['sst']
    k = data[:, keys.index('Turb_Kin_Energy')]
    assert relative(velocity, blend[:, None]*u_inf) < 2e-14
    assert relative(k, blend*model['kInfinity']) < 2e-14
    wall, far = detail['wallIDs'], detail['farIDs']
    assert np.all(velocity[wall] == 0) and np.all(k[wall] == 0)
    assert relative(velocity[far], np.broadcast_to(u_inf, velocity[far].shape)) < 2e-14
    c = settings((folder/'case.cfg').read_text())
    temperature = float(c['FREESTREAM_PRESSURE'])/(R*rho)
    assert max(abs(temperature[wall]-init['givenWallTemperatureK'])) < 2e-12
    assert max(abs(temperature[far]-float(c['FREESTREAM_TEMPERATURE']))) < 2e-12
    # 几何坐标与原网格的全量相等在正常导出及启动前验证器中核对。
    selected = np.unique(np.r_[np.linspace(0, len(data)-1, min(64, len(data)), dtype=int),
        wall[np.linspace(0, len(wall)-1, min(12, len(wall)), dtype=int)],
        far[np.linspace(0, len(far)-1, min(12, len(far)), dtype=int)],
        np.argsort(blend)[np.linspace(0, len(data)-1, min(32, len(data)), dtype=int)]])
    audit = folder/'initial_independent_audit'
    audit.mkdir(exist_ok=False)
    helper = ROOT/'outputs/v153_validation/native03/idealAirProbe.exe'
    assert sha(helper) == read(helper.parent/'runtime.json')['sha256']
    write_rows(audit/'oracle_input.csv', [('rhoT', float(rho[i]), float(temperature[i])) for i in selected])
    run(audit, 'oracle', [str(helper), '--oracle', str(audit/'oracle_input.csv'), str(audit/'oracle.csv')], 30)
    oracle = read_rows(audit/'oracle.csv')
    assert np.all(oracle['accepted'] == 1)
    errors = dict(thermalEnergy=relative(energy[selected], oracle['e']),
                  pressure=relative(np.full(len(selected), float(c['FREESTREAM_PRESSURE'])), oracle['p']))
    assert max(errors.values()) < 2e-10
    np.savetxt(str(audit/'node_ids.csv'), selected, fmt='%d', header='PointID', comments='')
    report = dict(passed=True, nodes=len(data), originalFaces=init['originalFaces'],
        independentSampleCount=len(selected), errors=errors,
        temperatureRangeK=[float(min(temperature)), float(max(temperature))],
        expectedVelocity=manifest['expectedVelocity'], inputMach=manifest['flight']['Mach'],
        preparationSeconds=read(folder/'ideal_air_preparation.json')['elapsedSeconds'],
        initialStateSHA256=sha(folder/'initial_flow.csv'), physicsQualified=False, newQualifiedArea=0)
    save(audit/'report.json', report)
    return report


def negative_contracts(source):
    """模块二。错误来源、准备失败和未适配保存方式必须在创建原生进程前拒绝。"""
    output = OUT/'contracts01'
    output.mkdir(exist_ok=False)
    rejected = []
    for name in ('incomplete', 'wrongMach', 'wrongVelocity', 'wrongK', 'sparseOutput', 'changedState'):
        folder = output/name
        folder.mkdir()
        for filename in ('mesh.su2', 'case.cfg', 'case_manifest.json', 'sst_contract.json',
                         'ideal_air_contract.json', 'ideal_air_initialization.json',
                         'ideal_air_preparation.json', 'ideal_air_wall_distance.npz', 'initial_flow.csv'):
            if filename == 'mesh.su2':
                os.link(str(source/filename), str(folder/filename))
            else:
                shutil.copyfile(str(source/filename), str(folder/filename))
        manifest = read(folder/'case_manifest.json')
        contract = manifest['idealAir']
        text = (folder/'case.cfg').read_text()
        if name == 'incomplete':
            save(folder/'ideal_air_preparation.json', dict(status='preparing', exitCode=None))
        elif name == 'wrongMach':
            text = set_value(text, 'MACH_NUMBER', '4')
        elif name == 'wrongVelocity':
            manifest['expectedVelocity'][0] *= .8
        elif name == 'wrongK':
            manifest['sst']['kInfinity'] *= 2
            save(folder/'sst_contract.json', manifest['sst'])
            manifest['sstContractSHA256'] = sha(folder/'sst_contract.json')
        elif name == 'sparseOutput':
            text = set_value(text, 'OUTPUT_WRT_FREQ', '2')
        else:
            with (folder/'initial_flow.csv').open('a') as stream:
                stream.write('\n')
        (folder/'case.cfg').write_text(text)
        contract['configurationSHA256'] = sha(folder/'case.cfg')
        save(folder/'ideal_air_contract.json', contract)
        manifest['idealAirContractSHA256'] = sha(folder/'ideal_air_contract.json')
        save(folder/'case_manifest.json', manifest)
        try:
            resolve(folder, text)
        except AssertionError as error:
            rejected.append(dict(name=name, error=str(error)))
        assert not (folder/'engine_process.json').exists()
    assert len(rejected) == 6
    for relative in ('outputs/v154_validation/cases01/laminar', 'outputs/v154_validation/cases01/sst',
                     'outputs/v145_validation/fixture02/case'):
        folder = ROOT/relative
        resolve(folder, (folder/'case.cfg').read_text())
    report = dict(passed=True, rejected=rejected, legacyReadContracts=3, nativeRuns=0)
    save(output/'report.json', report)
    return report


def main():
    """模块三。已有回执阻止重复任务，任何失败保留终态并停止后续数值阶段。"""
    assert read(OUT/'unit01/report.json')['passed']
    supervisor = OUT/'supervision'
    supervisor.mkdir(exist_ok=False)
    sources = [Path(__file__).resolve(), ROOT/'scripts/prepareSU2IdealAirCase.py',
        ROOT/'scripts/prepareSU2IdealAirWallDistance.py', ROOT/'scripts/idealAirProperties.py',
        ROOT/'scripts/su2IdealAirInterface.py', ROOT/'scripts/su2ThermalRuntime.py',
        ROOT/'scripts/runEmbeddedSU2.py', ROOT/'scripts/auditSU2IdealAirState.py',
        ROOT/'tests/testSU2IdealAirInitialState.py', ROOT/'tests/testSU2IdealAirInitialExport.m',
        ROOT/'docs/validation/v155/PLAN.md']
    locks = {str(path): sha(path) for path in sources}
    snapshot = OUT/'executed_sources'
    for path in sources:
        target = snapshot/path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
    save(OUT/'source_lock.json', locks)
    expression = "addpath(pwd); addpath('tests'); testSU2IdealAirInitialExport('outputs/v155_validation/cases01');"
    command = [sys.executable, '-X', 'utf8', str(ROOT/'scripts/runMatlabResearch.py'), 'export01', expression,
               '--folder', 'outputs/v155_validation/matlab']
    for base in ('near_wall_fixture', 'mesh'):
        for name in ('mesh_input.mat', 'hybrid_volume.mat', 'mesh_process.json', 'mesh_report.json'):
            command += ['--input', 'outputs/v145_validation/'+base+'/'+name]
    run(supervisor, 'matlab_export', command, 180)
    summaries = {}
    for name in ('small', 'real'):
        folder = OUT/'cases01'/name
        command = [sys.executable, '-X', 'utf8', str(ROOT/'scripts/prepareSU2IdealAirCase.py'),
                   '--folder', str(folder), '--sst', '--wall-distance-m', '.001',
                   '--provenance', 'v155 explicit fully turbulent wall-distance initial guess; transition unknown']
        run(supervisor, name+'_preparation', command, 180)
        summaries[name] = independent_initial(folder)
        if name == 'small':
            negative_contracts(folder)
            save(folder/'engine_request.json', dict(workDirectory=str(folder), configFile='case.cfg',
                meshFile='mesh.su2', threads=2, maximumSeconds=60, dryRun=False,
                sourceFiles=[path.relative_to(ROOT).as_posix() for path in sources]))
            run(supervisor, 'small_native', [sys.executable, '-X', 'utf8', str(ROOT/'scripts/runEmbeddedSU2.py'),
                '--request', str(folder/'engine_request.json')], 70)
            run(supervisor, 'small_audit', [sys.executable, '-X', 'utf8', str(ROOT/'scripts/auditSU2IdealAirState.py'),
                '--folder', str(folder)], 60)
    stable = {path: sha(Path(path)) == digest for path, digest in locks.items()}
    assert all(stable.values())
    total = sum(path.stat().st_size for path in OUT.rglob('*') if path.is_file())
    assert total < 160*1024**2, total
    save(OUT/'report.json', dict(passed=True, initialStates=summaries, sourceSHA256=locks,
        inputsStable=stable, nativeSmallRuns=1, nativeFullDomainRuns=0,
        fileLengthBytes=total, physicalAllocatedBytes=None, physicsQualified=False, newQualifiedArea=0))
    print(json.dumps(dict(passed=True, initialStates=summaries, fileLengthBytes=total)))


if __name__ == '__main__':
    main()
