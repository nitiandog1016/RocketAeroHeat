"""交付时只读核对：Run14新完成门槛兼容、预算和源码时点；不启动原生求解。"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'outputs/v180_validation'
sys.path.insert(0, str(ROOT / 'scripts'))
from auditEilmerCompletion import audit_completion


def read(name):
    return json.loads((BASE / name).read_text(encoding='utf-8-sig'))


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    start = time.monotonic()
    dest = BASE / 'closeout01.json'
    assert not dest.exists()
    # 新完成审核也用于Run14；只读已有日志/末载荷，避免无必要重跑控制。
    control = ROOT / 'outputs/v179_validation/control01'
    completion = audit_completion(control, control.parent / 'control01_run/process00/stdout.txt')
    assert completion['finalStep'] == 831 and completion['stopReason'] == 'maximum-wall-clock'
    matlab = read('eilmer_development01_process.json')
    assert matlab['status'] == 'terminal' and matlab['exitCode'] == 0
    assert len(matlab['sourceSHA256']) == 1885 and all(matlab['sourceStableAfter'].values())
    oldrun = read('real01_run/report.json')
    assert oldrun['passed'] and read('real01_inspection.json')['passed'] is False
    assert all(oldrun['commands'][0]['inputStable'].values())
    raw = (BASE / 'readEilmerNative-at-real-run.py.raw').read_bytes()
    assert digest(raw) == oldrun['commands'][0]['inputSHA256'][str(ROOT / 'scripts/readEilmerNative.py')]
    # 对MATLAB时点的完成审核源按准确SHA恢复；不把最终修正版冒充已执行旧版本。
    current = (ROOT / 'scripts/auditEilmerCompletion.py').read_bytes()
    old = current.replace(b"'absolute-global-residual-target','relative-global-residual-target'", b"'relative-residual'")
    expected = matlab['sourceSHA256']['scripts/auditEilmerCompletion.py']
    restored = digest(old) == expected
    if restored:
        target = BASE / 'auditEilmerCompletion-at-matlab.py.raw'
        assert not target.exists()
        target.write_bytes(old)
    result = dict(passed=True, run14Completion=completion, newNativeRuns=0, sourceLockCount=1885,
                  realRunInputCount=len(oldrun['commands'][0]['inputStable']),
                  originalWrongPassedReceiptPreserved=True, realRunReaderSHA256=digest(raw),
                  matlabAuditSourceRestored=restored, matlabAuditSourceExpectedSHA256=expected,
                  finalAuditSourceSHA256=digest(current), newQualifiedAreaM2=0,
                  elapsedSeconds=time.monotonic()-start)
    dest.write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
    # 分项列示；父监督已含MATLAB内部准备，不重复计入“独立总耗时”。
    rows = []
    for name, scope in [
        ('memory02/report.json', 'group-memory-test'), ('cancel01/report.json', 'cancel-timeout-test'),
        ('profile01.json', 'configuration-test'), ('small01_prepare/report.json', 'small-preparation-entry'),
        ('small01_run/report.json', 'small-run-entry'), ('small01_phase.json', 'small-phase-readonly'),
        ('real01_prepare/report.json', 'real-preparation-entry'), ('real01_run/report.json', 'real-failed-run-entry'),
        ('real01_inspection.json', 'real-failure-readonly'), ('completion01/report.json', 'completion-check-initial'),
        ('completion02/report.json', 'completion-check-final'), ('legacy_interface01/report.json', 'legacy-readonly'),
        ('matlab_supervision.json', 'matlab-parent-including-preparation'), ('memory_route/report.json', 'source-review-script'),
        ('closeout01.json', 'new-guard-run14-readonly')]:
        r = read(name)
        rows.append(dict(receipt='evidence/'+name, scope=scope, elapsedSeconds=r['elapsedSeconds']))
    failures = {name: read(name)['elapsedSeconds'] for name in ['memory01/normal/report.json','memory01/limit/report.json']}
    budget = dict(measuredEntries=rows, measuredEntrySumSeconds=sum(x['elapsedSeconds'] for x in rows),
                  initialMemoryTestProcessSeconds=failures, initialMemoryTestFullWallClockUnknown=True,
                  bootstrapFailedToolSeconds=read('bootstrap01_failure.json')['toolWallSeconds'],
                  matlabInternalTestSeconds=read('matlab01/report.json')['elapsedSeconds'],
                  realProcessSupervisorSeconds=oldrun['commands'][0]['elapsedSeconds'],
                  limits=dict(planCombinedProgramSeconds=2490, realEntrySeconds=1800, realParentSeconds=1860, matlabParentSeconds=120),
                  deviation='Early small preparation output limit was 1024MiB instead of planned 256MiB; actual about 30KB. No retrospective change.',
                  scope='Measured entry sum is not complete session time; nested MATLAB/process times shown separately. Development, shell coordination and source review reading time are not claimed as measured total.',
                  nativeRealAttempts=1, repeatedRealAttempts=0, engineBuilds=0, newMainUIRuns=0)
    (ROOT / 'docs/validation/v180/budget_receipt.json').write_text(json.dumps(budget, indent=2)+'\n', encoding='utf8')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
