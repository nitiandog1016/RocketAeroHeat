"""模块职责：审核与原面科研出口的准确输入、累计120秒预算及终态监督。"""
from pathlib import Path
import sys,shutil
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v257_validation';DOC=ROOT/'docs/validation/v257';CLOSED=False


def main(stage):
    if CLOSED:raise SystemExit('v257 fixed coupled audit/export closed; follow NEXT_DECISION')
    assert stage in ['audit','export'];guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    script=ROOT/'scripts'/('auditKinsolCoupledV257.py' if stage=='audit' else 'exportKinsolCoupledV257.py')
    locks=[Path(__file__).resolve(),script,DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',
        ROOT/'scripts/auditFullResidualV256.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',BASE/'process_control/report.json',BASE/'process_control/stdout.txt',
        ROOT/'outputs/v250_validation/rocket01/original_faces.txt',ROOT/'outputs/v256_validation/rocket02/initial_state.bin',
        ROOT/'outputs/v256_validation/rocket02/base_full_residual.bin']
    locks += [p for p in (BASE/'rocket01').glob('*') if p.is_file()]
    if stage=='export':locks += [BASE/'audit01/report.json',BASE/'audit01/surface_values.npz',ROOT/'web/kinsol-coupled-v257.js',
        ROOT/'tests/testKinsolViewerV257.js',ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz',
        ROOT/'outputs/v252_validation/workbench02/index.html',ROOT/'outputs/v252_validation/workbench02/case-data.js',
        ROOT/'web/viewer.js',ROOT/'web/section.js',ROOT/'web/export.js']
    used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']=='audit');assert used<120
    snap=BASE/('source_'+stage);snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix in ['.py','.md','.js','.html','.json','.txt'] and p.stat().st_size<2*1024**2:
            target=snap/(sha(p)+p.suffix)
            if not target.exists():shutil.copyfile(p,target)
    return supervise(guard,BASE,stage,[sys.executable,str(script)],120-used,1024,locks,ROOT,'audit')


if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
