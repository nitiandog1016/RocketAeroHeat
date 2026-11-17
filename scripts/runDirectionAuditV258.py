"""模块职责：一次方向审核及新原面诊断出口的累计60秒预算与准确输入监督。"""
import sys,shutil
from pathlib import Path
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v258_validation';DOC=ROOT/'docs/validation/v258';CLOSED=True


def main(stage):
    if CLOSED:raise SystemExit('v258 direction audit/export closed; follow NEXT_DECISION')
    assert stage in ['audit','export'];guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    script=ROOT/'scripts'/('auditDirectionConsistencyV258.py' if stage=='audit' else 'exportDirectionConsistencyV258.py')
    prev=ROOT/'outputs/v257_validation/rocket01'
    locks=[Path(__file__).resolve(),script,DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',
        ROOT/'scripts/auditFullResidualV256.py',ROOT/'scripts/auditKinsolCoupledV257.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',BASE/'process_control/report.json',ROOT/'outputs/v250_validation/rocket01/original_faces.txt']
    locks += [prev/n for n in ['initial_unknowns.bin','linear_refused_trial.bin','initial_state.bin','initial_full_residual.bin','scaling.csv']]
    locks += [p for p in (BASE/'rocket01').glob('*') if p.is_file()]
    if stage=='export':locks += [BASE/'audit01/report.json',BASE/'audit01/surface_values.npz',BASE/'audit01/direction_consistency.png',
        ROOT/'outputs/v257_validation/workbench01/index.html',ROOT/'outputs/v257_validation/workbench01/coupled-data.js',
        ROOT/'outputs/v252_validation/workbench02/case-data.js',ROOT/'web/kinsol-coupled-v257.js',ROOT/'web/direction-consistency-v258.js',
        ROOT/'tests/testDirectionViewerV258.js']
    used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='audit');assert used<60
    snap=BASE/('source_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix in ['.py','.js','.md','.html','.json'] and p.stat().st_size<2*1024**2:
            q=snap/(sha(p)+p.suffix)
            if not q.exists():shutil.copyfile(p,q)
    return supervise(guard,BASE,stage,[sys.executable,str(script)],60-used,512,locks,ROOT,'audit')


if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
