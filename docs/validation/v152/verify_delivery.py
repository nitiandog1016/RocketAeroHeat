"""直接读原网格、CSV及工作台导出，独立核对续算身份、测站和热积分。"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'outputs/v152_validation'
HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def raw_csv(path):
    with Path(path).open(encoding='utf-8-sig') as stream:
        labels = [key.strip() for key in next(csv.reader(stream, skipinitialspace=True))]
    data = np.loadtxt(str(path), delimiter=',', skiprows=1)
    return {name: data[:, i] for i, name in enumerate(labels)}


def verify():
    native = OUT/'acc_sbli01'
    receipt = read(native/'engine_process.json')
    manifest = read(native/'study_manifest.json')
    binding = read(native/'initialization_run_binding.json')
    report = read(OUT/'review01/report.json')
    assert receipt['status'] == 'terminal'
    assert all(receipt['inputsStable'].values()) and all(binding['sourceStableAfter'].values())
    locks = dict(receipt['sourceSHA256'], **binding['sourceSHA256'])
    for path, expected in locks.items():
        assert sha(path) == expected, path
    for name, expected in receipt['outputSHA256'].items():
        assert sha(native/name) == expected, name
    line = manifest['continuation']
    assert line == report['continuation'] and line['sourceSavedLocalIteration'] == 1199
    assert sha(native/'initial_flow.csv') == line['sourceStateSHA256'] == sha(line['sourceCheckpoint'])
    source = Path(line['sourceDirectory'])
    assert sha(native/'mesh.su2') == sha(source/'mesh.su2')
    history = raw_csv(native/'history.csv')
    last = int(history['Inner_Iter'][-1])
    assert report['lastRecordedIteration'] == last
    assert report['completedUpdatesInSeries'] == 1200+last+1 <= 12000
    assert np.array_equal(history['Inner_Iter'], np.arange(last+1))
    volume_files = sorted(native.glob('restart_flow_*.csv'))
    surface_files = sorted(native.glob('surface_flow_*.vtk'))
    assert len(surface_files) > len(volume_files)
    frame = raw_csv(volume_files[-1])
    assert int(volume_files[-1].stem.split('_')[-1]) == report['savedLocalIteration']
    with (native/'mesh.su2').open() as stream:
        assert stream.readline().strip() == 'NDIME= 2'
        count = int(stream.readline().split('=')[1])
        for _ in range(count):
            stream.readline()
        n = int(stream.readline().split('=')[1])
        points = np.array([list(map(float, stream.readline().split())) for _ in range(n)])
        markers = {}
        for _ in range(int(stream.readline().split('=')[1])):
            name = stream.readline().split('=')[1].strip()
            count = int(stream.readline().split('=')[1])
            markers[name] = np.array([list(map(int, stream.readline().split()))[1:] for _ in range(count)])
    wall = np.unique(markers['wall'])
    wall = wall[np.argsort(-points[wall, 0])]
    s = -points[wall, 0]
    s[s >= 0] /= np.cos(20*np.pi/180)
    predicted = {}
    for quantity, key in (('pressure', 'Pressure'), ('heat', 'Heat_Flux')):
        stations = report['currentFieldAssessment']['stations'][quantity]
        predicted[quantity] = np.interp(stations['S'], s, frame[key][wall])
        np.testing.assert_allclose(predicted[quantity], stations['predicted'], rtol=1e-13, atol=1e-10)
    # 线性q、r的乘积为二次多项式，两点Gauss对每条原线段精确积分。
    segments = markers['wall']
    length = np.linalg.norm(points[segments[:, 1], :2]-points[segments[:, 0], :2], axis=1)
    power = 0.
    for xi in (.5-.5/np.sqrt(3), .5+.5/np.sqrt(3)):
        q = (1-xi)*frame['Heat_Flux'][segments[:, 0]]+xi*frame['Heat_Flux'][segments[:, 1]]
        radius = (1-xi)*points[segments[:, 0], 1]+xi*points[segments[:, 1], 1]
        power += np.pi*np.sum(length*q*radius)
    power_error = abs(power/report['currentFieldAssessment']['revolvedHeatPowerW']-1)
    assert power_error < 2e-13
    saved_step = report['savedLocalIteration']
    rhs_error = max(abs(np.log10(np.sqrt(np.mean(frame['Residual_'+name]**2)))-history[key][saved_step])
                    for name, key in (('Density', 'rms[Rho]'), ('Energy', 'rms[RhoE]')))
    assert rhs_error < 2e-8
    steady = report['temporalStability']
    assert steady['sampling'] == 'native surface every 100 local indices'
    indices = [int(path.stem.split('_')[-1]) for path in surface_files]
    if steady['available']:
        selected = [index for index in indices if index >= steady['firstIteration']]
        assert selected[-1] == last and len(selected) == steady['frames']
        assert min(np.diff(selected)) > 0 and max(np.diff(selected)) <= 100
        assert last-selected[0] >= 1000
    # 工作台正式入口生成的CSV，以原始高精度体场再核对一遍。
    ui = OUT/'workbench01'
    ui_process = read(OUT/'ui_process.json')
    assert ui_process['status'] == 'terminal' and ui_process['exitCode'] == 0
    assert all(ui_process['sourceStableAfter'].values())
    assert sha(ui/'benchmark_review.json') == sha(OUT/'review01/report.json')
    assert sha(ui/'benchmark_review.json.png') == sha(OUT/'review01/convergence.png')
    with (ui/'benchmark_review.json.stations.csv').open(encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 51
    errors = []
    for row in rows:
        key = 'Pressure' if row['Quantity'] == 'pressure' else 'Heat_Flux'
        value = np.interp(float(row['S_m']), s, frame[key][wall])
        errors.append(abs(float(row['Predicted'])-value)/max(abs(value), 1))
        assert int(row['SourceLocalIteration']) == 1199
        assert int(row['SavedIteration']) == report['savedLocalIteration']
        assert int(row['CompletedUpdatesInSeries']) == report['completedUpdatesInSeries']
    assert max(errors) < 1e-12
    return dict(passed=True, nativeExitCode=receipt['exitCode'],
                lockedSources=len(locks), volumeFrames=len(volume_files), surfaceFrames=len(surface_files),
                completedUpdatesInSeries=report['completedUpdatesInSeries'],
                stationExportRelativeError=max(errors), powerRelativeError=power_error,
                nativeRHSLogError=rhs_error, matlabSourceCount=len(ui_process['sourceSHA256']),
                newQualifiedArea=0, outputPathFileLengths=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),
                diskAllocationMeasured=False)


if __name__ == '__main__':
    result = verify()
    (HERE/'delivery_verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
    print(json.dumps(result))
