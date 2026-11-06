"""用途：只重画既有筛选报告的温度散点图，不重跑CFD或审核，不覆盖旧证据。"""
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from screenNACAConstantCpDevelopment import read, save, sha, history_table, plot_history_temperature

start = time.monotonic()
out = ROOT / 'outputs/v174_validation/plot_final01'
out.mkdir(parents=True, exist_ok=False)
locks = {}
report = dict(renderCompleted=False, newCFDTasks=0, newAudits=0, newScreeningRuns=0,
    unsavedIntermediateTemperaturesKnown=False, temperatureInterpolationApplied=False)


def bind(path, expected=None):
    path = Path(path).resolve()
    value = sha(path)
    assert expected is None or value == expected
    locks[str(path)] = value


try:
    proof_path = ROOT / 'outputs/v174_validation/screen_real01/report.json'
    bind(proof_path)
    proof = read(proof_path)
    assert proof['reviewCompleted'] and proof['candidateAccepted'] is False
    old = proof_path.parent
    bind(old / proof['figure'], proof['figureSHA256'])
    bind(old / 'executed_screenNACAConstantCpDevelopment.py',
         proof['inputSHA256'][str(ROOT / 'scripts/screenNACAConstantCpDevelopment.py')])
    history = ROOT / 'outputs/v174_validation/real01/history.csv'
    bind(history, proof['inputSHA256'][str(history)])
    source = ROOT / 'scripts/screenNACAConstantCpDevelopment.py'
    for path in (Path(__file__).resolve(), source):
        bind(path)
        (out / ('executed_' + path.name)).write_bytes(path.read_bytes())
    rows, rejected, contiguous = history_table(history)
    assert contiguous and not rejected and len(rows) == proof['historyRows']
    plot_history_temperature(rows, proof['temperatureFrames'], out / 'development_history_temperature.png')
    report.update(renderCompleted=True, sourceScreenReport=str(proof_path),
        sourceScreenCandidateAccepted=False, physicalQualificationStillFalse=True,
        completeSavedLocalIterations=proof['completeSavedLocalIterations'],
        figure='development_history_temperature.png',
        figureSHA256=sha(out / 'development_history_temperature.png'))
except Exception as error:
    report.update(error=type(error).__name__ + ': ' + str(error), traceback=traceback.format_exc())
    raise
finally:
    stable = {name: Path(name).is_file() and sha(Path(name)) == digest for name, digest in locks.items()}
    report.update(inputSHA256=locks, inputsStable=stable, originalScreenEvidenceUnchanged=all(stable.values()),
        elapsedSeconds=time.monotonic() - start)
    save(out / 'report.json', report)
assert report['renderCompleted'] and all(stable.values())
print(dict(renderCompleted=True, inputsStable=len(stable), elapsedSeconds=report['elapsedSeconds']))
