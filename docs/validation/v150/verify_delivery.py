"""只读复核v150真实终态和导出，并归档本轮小型证据。"""
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'outputs/v150_validation'
ARCHIVE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def copy_exact(source, target):
    """小证据逐字归档，已有证据只有字节一致时才允许复用。"""
    data = source.read_bytes()
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        assert target.read_bytes() == data, str(target)
    else:
        target.write_bytes(data)
    return sha(target)


def verify():
    native = OUT / 'acc_sbli01'
    receipt = read(native / 'engine_process.json')
    assert receipt['status'] == 'terminal' and receipt['exitCode'] == 1
    assert receipt['stopReason'] is None and all(receipt['inputsStable'].values())
    for path, expected in receipt['sourceSHA256'].items():
        assert sha(path) == expected, path
    for name, expected in receipt['outputSHA256'].items():
        assert sha(native / name) == expected, name
    assert 'NaN detected' in (native / 'engine.log').read_text()
    assert not list(native.glob('restart_flow*.csv'))
    assert not list(native.glob('surface_flow*'))
    with (native / 'history.csv').open() as stream:
        rows = [{key.strip(): value for key, value in row.items()}
                for row in csv.DictReader(stream, skipinitialspace=True)]
    iterations = [int(row['Inner_Iter']) for row in rows]
    assert iterations == list(range(326))
    latest = read(OUT / 'review02/report.json')
    original = read(OUT / 'review01/report.json')
    assert latest['lastRecordedIteration'] == 325
    assert latest['sourceLocalIteration'] == 1
    assert latest['savedLocalIteration'] is None and latest['currentFieldAssessment'] is None
    assert not latest['hasSavedNewField'] and not latest['validationPassed']
    assert not latest['convergenceAcceptanceImplemented']
    assert len(latest['missingCurrentChecks']) == 8
    assert latest['sourceInitialFieldAssessment'] == original['sourceInitialFieldAssessment']
    for field, value in latest['finalRecordedResiduals'].items():
        assert float(rows[-1][field]) == value
    for review in (original, latest):
        assert sha(review['historyPlot']['path']) == review['historyPlot']['sha256']
    assert sha(OUT / 'workbench02/benchmark_review.json') == sha(OUT / 'review01/report.json')
    assert sha(OUT / 'workbench02/benchmark_review.json.png') == sha(OUT / 'review01/convergence.png')
    ui = read(OUT / 'workbench02/report.json')
    for key in ('mainWorkbenchEntry', 'actualFailureShown', 'inputStatePreserved',
                'fakeSavedFrameRejected', 'exportExact', 'returnedToModelStep'):
        assert ui[key], key
    assert not ui['nativeScreenshotObserved']
    assert not (OUT / 'layout_done.json').exists()

    # 每次实际MATLAB执行的源锁，按准确旧源恢复；接续后的元数据修订独立列出。
    changed = {}
    runs = {}
    for name, code in (('workbench01', 1), ('workbench02', 0), ('layout01', 0)):
        run = read(OUT / (name + '_process.json'))
        supervisor = read(OUT / (name + '_supervisor/process.json'))
        assert run['status'] == supervisor['status'] == 'terminal'
        assert run['exitCode'] == supervisor['exitCode'] == code
        assert all(run['sourceStableAfter'].values())
        assert all(supervisor['sourceStableAfter'].values())
        for rel, expected in run['sourceSHA256'].items():
            path = ROOT / rel
            if sha(path) != expected:
                candidates = [ARCHIVE / 'executed_sources' / rel]
                if rel == 'tests/testSU2BenchmarkReview.m':
                    candidates.append(ARCHIVE / 'failures/testSU2BenchmarkReview_v1.m')
                matches = [p for p in candidates if p.is_file() and sha(p) == expected]
                assert matches, rel
                changed[name + '/' + rel] = str(matches[0].relative_to(ROOT))
        runs[name] = dict(exitCode=code, seconds=supervisor['seconds'],
                          sourcesAndInputs=len(run['sourceSHA256']), stableDuringRun=True)

    selected = [
        'acc_sbli01/case.cfg', 'acc_sbli01/engine_process.json',
        'acc_sbli01/study_manifest.json', 'acc_sbli01/engine.log', 'acc_sbli01/history.csv',
        'review01/report.json', 'review01/convergence.png', 'review02/report.json',
        'workbench02/report.json', 'run_bounded.py']
    for name in runs:
        selected += [name + '_process.json', name + '_request.json', name + '_stdout.log',
                     name + '_supervisor/process.json', name + '_supervisor/stdout.log']
    archived = {rel: copy_exact(OUT / rel, ARCHIVE / 'evidence' / rel) for rel in selected}
    lengths = sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
    assert lengths < 120 * 1024**2
    return dict(passed=True, nativeRuns=1, nativeExitCode=1, lastRecordedIteration=325,
                newSavedFields=0, numericalRerun=False, matlabRuns=runs,
                exactArchivedExecutionSources=changed, sourceFieldAssessmentUnchanged=True,
                exportJSONAndPNGExact=True, nativeScreenshotObserved=False,
                convergenceAcceptanceImplemented=False, newQualifiedArea=0,
                outputPathFileLengthsBytes=lengths, includesHardlinkedInputs=True,
                physicalAllocatedBytesMeasured=False, archivedSHA256=archived)


if __name__ == '__main__':
    result = verify()
    (ARCHIVE / 'delivery_verification.json').write_bytes(
        (json.dumps(result, indent=2) + '\n').encode('utf-8'))
    print(json.dumps({k: result[k] for k in (
        'passed', 'nativeRuns', 'nativeExitCode', 'lastRecordedIteration', 'newSavedFields',
        'nativeScreenshotObserved', 'outputPathFileLengthsBytes')}))
