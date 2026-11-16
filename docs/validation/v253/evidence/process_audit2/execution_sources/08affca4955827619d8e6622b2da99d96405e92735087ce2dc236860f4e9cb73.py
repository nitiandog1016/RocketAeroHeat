"""模块职责：具名真实入口任务的固定监督；不重启旧MATLAB/CFD分支。"""
from pathlib import Path
import sys
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read
from terminalCapacityV249 import supervise
BASE=ROOT/'outputs/v253_validation';DOC=ROOT/'docs/validation/v253';CLOSED=False
INPUT=[ROOT/'outputs/v222_validation/geometry01/cache.npz',ROOT/'outputs/v228_validation/prepare03/region.json',
       ROOT/'outputs/v226_validation/matlab02/rocket_metadata.json']+[ROOT/'outputs/v252_validation/matlab01'/(k+'.json') for k in ('base','fine','history')]

def main(stage):
    if CLOSED:raise SystemExit('v253 entropy inlet task closed; follow NEXT_DECISION')
    assert stage in ('produce','produce2','produce3','audit','audit2','export')
    BASE.mkdir(exist_ok=True);guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
    original=guard.size_bytes;guard.size_bytes=lambda p:original(BASE) if p==BASE or p.name.startswith('process_') else 0
    script=ROOT/'scripts'/({'produce':'produceEntropyInletV253.py','produce2':'produceEntropyInletV253.py','produce3':'produceEntropyInletV253.py','audit':'auditEntropyInletV253.py','audit2':'auditEntropyInletV253.py','export':'exportEntropyInletV253.py'}[stage])
    locks=INPUT+[Path(__file__).resolve(),script,DOC/'PLAN.md',ROOT/'scripts/phaseAwareEntropyInlet.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'scripts/phaseAwareViscousAir.py',
        ROOT/'data/thermo/phase_aware_dry_air_v1.json',ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py']
    if stage=='produce3':locks += [BASE/'inlet01/shock_mass_source.json',BASE/'inlet01/controls.json']
    if stage in ('audit','audit2','export'):locks += [BASE/'inlet01/shock_mass_source.json',BASE/'inlet02/inlets.json',BASE/'inlet02/inlet_profiles.npz',BASE/'inlet02/report.json']
    kind='offline' if stage.startswith('produce') else 'audit'
    used=sum(read(p)['elapsedSeconds'] for p in BASE.glob('process_*/report.json') if read(p)['budgetKind']==kind)
    return supervise(guard,BASE,stage,[sys.executable,str(script),stage],min(120,(120 if kind=='offline' else 180)-used),128,locks,ROOT,kind)

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
