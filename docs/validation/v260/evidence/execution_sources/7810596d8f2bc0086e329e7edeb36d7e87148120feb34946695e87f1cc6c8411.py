"""模块职责：量热helper独立审核的准确输入与60秒预算；不再次查询原生接口。"""
import sys,shutil
from pathlib import Path
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,sha
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v260_validation';DOC=ROOT/'docs/validation/v260';CLOSED=False


def main():
    if CLOSED:raise SystemExit('v260 caloric audit closed; follow NEXT_DECISION')
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;measure=guard.size_bytes;guard.size_bytes=lambda p:measure(BASE) if p==BASE or p.name.startswith('process_') else 0
    script=ROOT/'scripts/auditCaloricInverseV260.py'
    locks=[Path(__file__).resolve(),script,DOC/'PLAN.md',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',
        ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json',
        ROOT/'outputs/v220_validation/source/SU2_CFD/include/fluid/PhaseAirTable.hpp',BASE/'process_query/report.json']
    locks += [p for p in (BASE/'case01').glob('*') if p.is_file()]
    snap=BASE/'source_audit';snap.mkdir(exist_ok=False)
    for p in locks:
        if p.suffix in ['.py','.md']:shutil.copyfile(p,snap/(sha(p)+p.suffix))
    return supervise(guard,BASE,'audit',[sys.executable,str(script)],60,64,locks,ROOT,'audit')


if __name__=='__main__':raise SystemExit(0 if main() else 1)
