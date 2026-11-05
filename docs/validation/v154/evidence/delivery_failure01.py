"""独立核对v154导出、真实终态、准确执行源和小型公开证据。"""
from pathlib import Path
from datetime import datetime
import csv
import hashlib
import json
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from validateIdealAir import save, sha
from su2ThermalRuntime import read

OUT = ROOT/'outputs/v154_validation'
DEST = ROOT/'docs/validation/v154'


def main():
    evidence = DEST/'evidence'
    evidence.mkdir(exist_ok=False)
    copied = {}

    def include(relative):
        source = OUT/relative
        target = evidence/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(str(source), str(target))
        assert sha(source) == sha(target)
        copied[target.relative_to(ROOT).as_posix()] = sha(target)

    errors = {}
    runs = {}
    for mode, directory in (('laminar', 'ideal_air_audit02'), ('sst', 'ideal_air_audit')):
        folder = OUT/'cases01'/mode
        receipt = read(folder/'engine_process.json')
        assert receipt['status'] == 'terminal' and receipt['exitCode'] == 0 and all(receipt['inputsStable'].values())
        for path, expected in receipt['inputSHA256'].items():
            assert sha(Path(path)) == expected, path
        for name, expected in receipt['outputSHA256'].items():
            assert sha(folder/name) == expected, name
        audit = read(folder/directory/'report.json')
        assert audit['interfacePassed'] and not audit['originalResidualThresholdPassed']
        # 使用Python CSV独立读取MATLAB导出，核对全部九个原始字段及温度/方法角色。
        reference = np.genfromtxt(str(folder/directory/'wall.csv'), names=True, delimiter=',')
        with (OUT/'ui02'/(mode+'.csv')).open(encoding='utf-8-sig', newline='') as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == len(reference) == 65
        error = 0.
        for field in reference.dtype.names:
            actual = np.array([float(row[field]) for row in rows])
            error = max(error, float(max(abs(actual-reference[field])/np.maximum(1., abs(reference[field])))))
        assert error < 2e-12
        assert all(float(row['given_Tw_K']) == 300 and np.isnan(float(row['solvedWallTemperature_K'])) for row in rows)
        assert all(row['PhysicsQualified'] == '0' and row['SavedLocalIteration'] == '9' for row in rows)
        assert all('OTH' in row['Method'] and 'Roe' in row['Method'] for row in rows)
        if mode == 'sst':
            assert all('SST-2003m' in row['Method'] for row in rows)
        errors[mode] = error
        runs[mode] = dict(elapsedSeconds=receipt['elapsedSeconds'], exitCode=0, nodes=audit['nodes'],
                          localIteration=audit['checkpointLocalIteration'], errors=audit['errors'])
        for name in ('case.cfg', 'case_manifest.json', 'ideal_air_contract.json', 'engine_request.json',
                     'engine_process.json', 'engine_stdout.log', 'history.csv', 'ideal_air_review.json'):
            include(Path('cases01')/mode/name)
        if mode == 'sst':
            include(Path('cases01')/mode/'sst_contract.json')
        for name in ('report.json', 'wall.csv', 'input_locks.json', 'executed_audit.py'):
            include(Path('cases01')/mode/directory/name)
        for name in ('oracle_process.json', 'current_process.json'):
            include(Path('cases01')/mode/'ideal_air_audit'/name)
        include(Path('ui02')/(mode+'.csv'))
        include(Path('ui02')/(mode+'.csv.json'))
    for version in ('adapter01', 'adapter02'):
        for path in (OUT/version).rglob('*'):
            if path.is_file() and path.suffix in ('.json', '.csv', '.hpp', '.cpp', '.py', '.md', '.log'):
                include(path.relative_to(OUT))
    for name in ('compile_process.json', 'source_restore.json', 'build_sources.json', 'source_tree.json',
                 'runtime_installed.json', 'upstream_audit.json', 'syntax_CConfig_process.json',
                 'syntax_CEulerSolver_process.json', 'syntax_CNSSolver_process.json'):
        include(Path('build01')/name)
    assert all(read(OUT/'build01/source_restore.json').values())
    for path, expected in read(OUT/'build01/build_sources.json').items():
        assert sha(Path(path)) == expected, path
    for path in (OUT/'build01/executed_sources').rglob('*'):
        if path.is_file():
            relative = path.relative_to(OUT/'build01/executed_sources')
            source = ROOT/'tools/runtime/su2-thermal-build/source'/relative
            assert sha(path) == read(OUT/'build01/source_tree.json')[str(source)]
    for relative in ('patch01/prepareSU2IdealAirPatch.py', 'patch02/ideal_air.diff', 'patch02/patch_manifest.json',
                     'contracts01/failed_test.py', 'contracts02/report.json', 'ui02/capture_done.json',
                     'cases01/laminar/ideal_air_audit/executed_audit.py'):
        include(Path(relative))
    matlab = {}
    for name in ('ui01', 'ui02', 'legacy01'):
        receipt = read(OUT/(name+'_process.json'))
        assert receipt['status'] == 'terminal' and all(receipt['sourceStableAfter'].values())
        expected_code = 1 if name == 'ui01' else 0
        assert receipt['exitCode'] == expected_code
        elapsed = (datetime.fromisoformat(receipt['endedUTC'])-datetime.fromisoformat(receipt['startedUTC'])).total_seconds()
        matlab[name] = dict(exitCode=expected_code, elapsedSeconds=elapsed, sourceCount=len(receipt['sourceSHA256']))
        include(Path(name+'_process.json'))
        include(Path(name+'_stdout.log'))
        if name != 'ui01':
            include(Path(name)/'report.json')
    # 既存OpenCode改动仅保存本地准确副本，公开提交继续隔离这些工作区改动。
    local = OUT/'ui_existing_sources'
    local.mkdir(exist_ok=False)
    for name in ('web/viewer.js', 'web/section.js', 'tests/verifyWebSectionExport.py'):
        path = ROOT/name
        assert sha(path) == read(OUT/'ui02_process.json')['sourceSHA256'][name]
        target = local/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(str(path), str(target))
    files = [path for path in OUT.rglob('*') if path.is_file()]
    length = sum(path.stat().st_size for path in files)
    case_length = sum(path.stat().st_size for path in (OUT/'cases01').rglob('*') if path.is_file())
    assert case_length < 80*1024**2
    runtime = ROOT/'tools/runtime/su2/8.6.0-rah-ideal-air-v1'
    runtime_length = sum(path.stat().st_size for path in runtime.rglob('*') if path.is_file())
    report = dict(passed=True, nativeRuns=runs, MATLAB=matlab, exportMaximumRelativeError=errors,
        outputFileLengthBytes=length, caseFileLengthBytes=case_length, runtimeFileLengthBytes=runtime_length,
        physicalAllocatedBytes=None, independentDataDeleted=False, packageCreated=False,
        sourceSHA256={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}, evidenceSHA256=copied,
        physicsQualified=False, newQualifiedArea=0)
    save(DEST/'independent_delivery.json', report)
    print(json.dumps({key: report[key] for key in ('passed', 'exportMaximumRelativeError',
                     'outputFileLengthBytes', 'caseFileLengthBytes', 'runtimeFileLengthBytes')}))


if __name__ == '__main__':
    main()
