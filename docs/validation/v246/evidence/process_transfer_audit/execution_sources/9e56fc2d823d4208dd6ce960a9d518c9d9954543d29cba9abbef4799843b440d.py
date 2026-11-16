"""模块职责：只读实际层间体积积分，核对真实限制/延拓的质量和物理能量身份。"""
from pathlib import Path
import csv,sys
import numpy as np
import runFASRocketControlV246 as run
from runSU2RocketPhaseV222 import ROOT,read,save
BASE=run.BASE;DOC=run.DOC;CLOSED=False


def review():
    out=BASE/'transfer_audit01';out.mkdir(exist_ok=False)
    with (BASE/'control_observe01/levels.csv').open() as f:rows=list(csv.DictReader(f))
    index={(int(r['update']),r['stage'],int(r['level'])):r for r in rows};assert len(index)==len(rows)==100
    names=['volumeRho','volumeMx','volumeMy','volumeMz','volumeEphysical'];errs={}
    def values(step,label,level):return np.array([float(index[(step,label,level)][n]) for n in names])
    def check(label,a,b):
        error=float(max(abs(a-b)/np.maximum(1,abs(b))));assert error<2e-12,(label,error);errs[label]=max(errs.get(label,0),error)
    for step in range(1,6):
        for label in ('restrict_in','restrict_out','coarse_preprocessed','forcing_in','forcing_out','correction_in','correction_out','prolong_in','prolong_out','flow_after_time'):
            for level in (0,1):assert (step,label,level) in index
        fine=values(step,'restrict_in',0);coarse=values(step,'restrict_out',1)
        check('restrictedMassPhysicalEnergy',coarse[[0,4]],fine[[0,4]])
        check('restrictionLeavesFineState',values(step,'restrict_out',0),fine)
        check('forcingLeavesFineState',values(step,'forcing_in',0),values(step,'forcing_out',0))
        check('forcingLeavesCoarseState',values(step,'forcing_in',1),values(step,'forcing_out',1))
        check('correctionConstructionLeavesFineState',values(step,'correction_in',0),values(step,'correction_out',0))
        check('correctionConstructionLeavesCoarseState',values(step,'correction_in',1),values(step,'correction_out',1))
        before=values(step,'prolong_in',0);coarse=values(step,'prolong_in',1);after=values(step,'prolong_out',0)
        check('prolongedMassPhysicalEnergy',after[[0,4]],(.5*before+.5*coarse)[[0,4]])
    report=dict(passed=True,completeCycles=5,realLayerRows=100,volumeMassEnergyIdentitiesPassed=True,
        maxScaledErrors=errs,physicalEnergyOffsetJkg=300000,finalCycleResidualNotSaved=True,
        steadyBalanceVerified=False,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'report.json',report);print(report)


def main():
    if CLOSED:raise SystemExit('v246 transfer review closed; follow NEXT_DECISION')
    g=run.prep.guard;g.BASE=BASE;g.CASE=BASE;g.DOC=DOC;size=g.size_bytes;g.size_bytes=lambda p:size(BASE) if p==BASE or p.name.startswith('process_') else 0
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',BASE/'control_observe01/levels.csv']
    return g.supervise('transfer_audit',[sys.executable,str(Path(__file__).resolve()),'task'],30,1536,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v246 transfer review closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':review()
    else:raise SystemExit(0 if main() else 1)
