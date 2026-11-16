"""模块职责：原MATLAB终态缺失时只读恢复；固定180秒独立审核/科研出口预算。"""
from pathlib import Path
import sys
import runPhaseNonsimilarV252 as run
from runSU2RocketPhaseV222 import ROOT,read
from terminalCapacityV249 import supervise
CLOSED=False

def main(stage):
    if CLOSED:raise SystemExit('v252 audit tasks closed; follow NEXT_DECISION')
    assert stage in ('recover','audit','export')
    g=run.setup();base=run.BASE
    script=ROOT/'scripts'/('exportPhaseNonsimilarV252.py' if stage=='export' else 'auditPhaseNonsimilarV252.py')
    locks=[Path(__file__).resolve(),script,ROOT/'scripts/runPhaseNonsimilarV252.py',
           ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
           ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',run.DOC/'PLAN.md']
    locks+=list((base/'matlab01').glob('*.json'))+[ROOT/n for n in run.OWN+run.LIB]
    used=sum(read(p)['elapsedSeconds'] for p in base.glob('process_*/report.json') if read(p)['budgetKind']=='audit')
    return supervise(g,base,stage,[sys.executable,str(script),stage],min(120,180-used),256,locks,ROOT,'audit')

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
