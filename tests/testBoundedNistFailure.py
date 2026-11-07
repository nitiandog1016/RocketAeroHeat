"""新映射实际失败的拒绝/身份核验，旧表域与模型不改。"""
from pathlib import Path
import json,sys,subprocess,shutil,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareOpenFoamControl import sha,save
from readOpenFoamNistRocketState import verify_identity


def main():
    start=time.monotonic();base=ROOT/'outputs/v203_validation';out=base/'failure_tests01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_test.py')
    run=json.loads((base/'rocket_run01/report.json').read_text());admission=json.loads((base/'rocket_run01/admission.json').read_text());audit=json.loads((base/'rocket_audit01/report.json').read_text())
    assert audit['passed'] and not audit['candidatePassed'] and audit['rejectedByPressureBeforeRootSolve']
    assert not run['saved'] and not run['passed'] and not audit['priorV202TemperatureFailureRootCauseEstablished']
    rows=[]
    try:verify_identity(base/'rocket01',run,admission)
    except ValueError:rows.append('v202_reader_rejects_new_inverse_identity')
    else:raise AssertionError('旧读取器不能隐式接收新映射')
    for name in ['runBoundedNistRocket.py','checkBoundedNistDispatch.py']:
        p=subprocess.run([sys.executable,str(ROOT/'scripts'/name)],capture_output=True,timeout=10)
        assert p.returncode!=0 and b'ValueError' in p.stderr and b'ModuleNotFoundError' not in p.stderr
        rows.append('closed_'+name)
    protected=json.loads((ROOT/'outputs/v176_validation/protected_initial.json').read_text(encoding='utf-8'))
    assert len(protected)==29 and all(sha(ROOT/n)==h for n,h in protected.items())
    save(out/'report.json',dict(passed=True,checks=rows,contractRejections=8,contractRejectionEvidence='outputs/v203_validation/bounded01/report.json',
        exactProtectedFiles=29,CFDRuns=0,elapsedSeconds=time.monotonic()-start))
    print(json.dumps(dict(passed=True,checks=len(rows),protected=29,elapsedSeconds=time.monotonic()-start)))


if __name__=='__main__':main()
