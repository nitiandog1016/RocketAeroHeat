"""独立读原CSV、原壁线段和工作台导出，核对v151交付及准确执行源。"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/v151_validation'
HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def raw_csv(path):
    with path.open() as stream:
        keys = [s.strip() for s in next(csv.reader(stream))]
    data = np.loadtxt(str(path), delimiter=',', skiprows=1)
    return {key: data[:, i] for i, key in enumerate(keys)}


def copy_exact(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    data = source.read_bytes()
    if target.exists():
        assert target.read_bytes() == data, str(target)
    else:
        target.write_bytes(data)
    return sha(target)


def verify():
    native = OUT/'acc_sbli01'
    receipt = read(native/'engine_process.json')
    binding = read(native/'initialization_run_binding.json')
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == binding['exitCode'] == 0
    assert all(receipt['inputsStable'].values()) and all(binding['sourceStableAfter'].values())
    for path, expected in dict(receipt['sourceSHA256'], **binding['sourceSHA256']).items():
        assert sha(path) == expected, path
    for name, expected in receipt['outputSHA256'].items():
        assert sha(native/name) == expected, name
    checkpoints = sorted(native.glob('restart_flow_*.csv'))
    assert [int(p.stem.split('_')[-1]) for p in checkpoints] == list(range(100, 1200, 100))+[1199]
    assert sha(native/'restart_flow.csv') == sha(checkpoints[-1])

    # 独立读取SU2原边界，原坐标和wall标记直接决定导出的空间采样。
    with (native/'mesh.su2').open() as stream:
        assert stream.readline().strip() == 'NDIME= 2'
        for _ in range(int(stream.readline().split('=')[1])):
            stream.readline()
        points = np.array([list(map(float, stream.readline().split()))
                           for _ in range(int(stream.readline().split('=')[1]))])
        markers = {}
        for _ in range(int(stream.readline().split('=')[1])):
            name = stream.readline().split('=')[1].strip()
            rows = [list(map(int, stream.readline().split()))
                    for _ in range(int(stream.readline().split('=')[1]))]
            markers[name] = np.unique(np.array(rows)[:, 1:])
    ids = markers['wall']
    initial = raw_csv(native/'initial_flow.csv')
    rho = initial['Density']
    u = np.column_stack((initial['Momentum_x'], initial['Momentum_y']))/rho[:, None]
    thermal = initial['Energy']/rho-np.sum(u*u, axis=1)/2-initial['Turb_Kin_Energy']
    temperature = thermal*.4/287.058
    wall_error = float(max(abs(temperature[ids]-311)))
    assert wall_error < 1e-11 and np.all(u[ids] == 0) and min(temperature) > 0
    assert np.array_equal(points[:, :2], np.column_stack((initial['x'], initial['y'])))
    condition = read(native/'study_manifest.json')['benchmark']['conditions']
    pressure_error = float(max(abs(rho*287.058*temperature/condition['pressurePa']-1)))
    assert pressure_error < 2e-12
    latest = raw_csv(checkpoints[-1])
    report = read(OUT/'review01/report.json')
    assert report['savedLocalIteration'] == report['lastRecordedIteration'] == 1199
    assert report['sourceLocalIteration'] is None and not report['validationPassed']
    assert report['convergenceAcceptanceImplemented'] and len(report['acceptance']['failed']) == 21
    assert not report['acceptance']['pending'] and not report['physicsQualified']
    residual_error = 0.
    for name, key in [('Density', 'rms[Rho]'), ('Energy', 'rms[RhoE]')]:
        raw = np.log10(np.sqrt(np.mean(latest['Residual_'+name]**2)))
        residual_error = max(residual_error, abs(raw-report['finalRecordedResiduals'][key]))
    assert residual_error < 1e-8

    x = latest['x'][ids]
    arc = np.where(x >= 0, -x, -x/np.cos(np.deg2rad(20)))
    order = np.argsort(arc)
    with (OUT/'workbench02/benchmark_review.json.stations.csv').open(encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    maximum_error = 0.
    for row in rows:
        field = 'Pressure' if row['Quantity'] == 'pressure' else 'Heat_Flux'
        expected = np.interp(float(row['S_m']), arc[order], latest[field][ids][order])
        maximum_error = max(maximum_error, abs(expected-float(row['Predicted']))/max(abs(expected), 1.))
        assert row['Unit'] == ('Pa' if field == 'Pressure' else 'W/m2')
        assert row['Method'] == report['model'] and int(row['SavedIteration']) == 1199
        assert row['ValidationPassed'] in ('0', 'false')
    assert len(rows) == 51 and maximum_error < 2e-12

    runs = {}
    for name in ('workbench01', 'workbench02'):
        run = read(OUT/(name+'_process.json'))
        sup = read(OUT/(name+'_supervisor/process.json'))
        assert run['status'] == sup['status'] == 'terminal'
        assert run['exitCode'] == sup['exitCode'] == 0
        assert all(run['sourceStableAfter'].values()) and all(sup['sourceStableAfter'].values())
        for relative, expected in run['sourceSHA256'].items():
            current = ROOT/relative
            if sha(current) != expected:
                old = HERE/'ui_before_scientific_format'/current.name
                assert old.is_file() and sha(old) == expected, relative
        assert sha(OUT/name/'benchmark_review.json') == sha(OUT/'review01/report.json')
        assert sha(OUT/name/'benchmark_review.json.png') == sha(OUT/'review01/convergence.png')
        runs[name] = dict(seconds=sup['seconds'], exitCode=0, stableSources=len(run['sourceSHA256']),
                          captureMarkerBeforeDeadline=read(OUT/name/'report.json')['nativeScreenshotObserved'])
    selected = ['acc_sbli01/'+name for name in (
        'case.cfg', 'engine.log', 'engine_process.json', 'study_manifest.json',
        'initialization_manifest.json', 'initialization_run_binding.json', 'history.csv')]
    selected += ['review01/report.json', 'review01/convergence.png', 'checkpoint_trend.json', 'run_bounded.py']
    for name in runs:
        selected += [name+'/report.json', name+'/capture_done.json', name+'_process.json',
                     name+'_request.json', name+'_stdout.log', name+'_supervisor/process.json',
                     name+'_supervisor/stdout.log']
    selected += ['workbench02/benchmark_review.json.stations.csv']
    archived = {name: copy_exact(OUT/name, HERE/'evidence'/name) for name in selected}
    sources = ['scripts/prepareSU2AxisymmetricInitialState.py', 'scripts/su2BenchmarkAcceptance.py',
               'scripts/auditSU2AxisymmetricStudy.py', 'scripts/su2AxisymmetricAnalysis.py',
               'createSU2BenchmarkReview.m', 'readSU2BenchmarkReview.m',
               'tests/testSU2AxisymmetricInitialization.py', 'tests/testSU2BenchmarkAcceptanceReview.m']
    source_hashes = {name: copy_exact(ROOT/name, HERE/'executed_sources'/name) for name in sources}
    return dict(passed=True, nativeRuns=1, updates=1200, nativeSeconds=receipt['elapsedSeconds'],
                savedFields=len(checkpoints), initializationWallTemperatureErrorK=wall_error,
                initializationPressureRelativeError=pressure_error, rawResidualLogError=residual_error,
                stationExportMaximumRelativeError=maximum_error, matlab=runs,
                nativeWindowObservedInTool=True, latestCaptureMarkerLate=True,
                newQualifiedArea=0, outputPathFileLengthsBytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),
                allocatedDiskBytesMeasured=False, evidenceSHA256=archived, executionSourceSHA256=source_hashes)


if __name__ == '__main__':
    result = verify()
    (HERE/'delivery_verification.json').write_bytes((json.dumps(result, indent=2)+'\n').encode('utf8'))
    print(json.dumps({key: result[key] for key in ('passed', 'savedFields', 'stationExportMaximumRelativeError',
                                                 'rawResidualLogError', 'outputPathFileLengthsBytes')}))
