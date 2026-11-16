"""模块职责：原MATLAB终态缺失时只读恢复；固定180秒独立审核/科研出口预算。"""
from pathlib import Path
import sys
import runPhaseNonsimilarV252 as run
from runSU2RocketPhaseV222 import ROOT,read
from terminalCapacityV249 import supervise
CLOSED=True

def main(stage):
    if CLOSED:raise SystemExit('v252 audit tasks closed; follow NEXT_DECISION')
    assert stage in ('recover','audit','export','export2','guards')
    g=run.setup();base=run.BASE
    script=ROOT/'scripts'/('exportPhaseNonsimilarV252.py' if stage.startswith('export') else 'auditPhaseNonsimilarV252.py')
    locks=[Path(__file__).resolve(),script,ROOT/'scripts/runPhaseNonsimilarV252.py',
           ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',
           ROOT/'scripts/runSU2RocketPhaseV222.py',ROOT/'scripts/terminalCapacityV249.py',run.DOC/'PLAN.md']
    locks+=list((base/'matlab01').glob('*.json'))+[ROOT/n for n in run.OWN+run.LIB]
    if stage=='guards':
        import shutil
        out=base/'guards01';out.mkdir(exist_ok=False)
        for name in run.OWN+run.LIB+['tests/testPhaseNonsimilarGuardsV252.m']:
            dest=out/'executed_source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,dest)
        locks+=run.INPUT+[ROOT/'tests/testPhaseNonsimilarGuardsV252.m']
        expression="restoredefaultpath;maxNumCompThreads(1);set(groot,'defaultFigureVisible','off');addpath(pwd);addpath(fullfile(pwd,'tests'));testPhaseNonsimilarGuardsV252('"+out.as_posix()+"');"
        return supervise(g,base,stage,['D:/matlab/bin/matlab.exe','-batch',expression],60,256,locks,ROOT,'matlab')
    if stage.startswith('export'):
        locks += [ROOT/n for n in ('web/index.html','web/viewer.js','web/section.js','web/export.js',
            'web/phase-profiles-v252.js','web/phase-profiles-v252.html','tests/testPhaseNonsimilarViewerV252.js')]
    used=sum(read(p)['elapsedSeconds'] for p in base.glob('process_*/report.json') if read(p)['budgetKind']=='audit')
    return supervise(g,base,stage,[sys.executable,str(script),stage],min(120,180-used),256,locks,ROOT,'audit')

if __name__=='__main__':raise SystemExit(0 if main(sys.argv[1]) else 1)
