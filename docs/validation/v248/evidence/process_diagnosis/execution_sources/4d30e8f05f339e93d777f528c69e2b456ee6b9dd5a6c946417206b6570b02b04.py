"""模块职责：仅用已保存末空间残差，量化恢复节点占比；不把相关性当因果，不启动CFD。"""
from pathlib import Path
import csv,sys
import numpy as np
import runFASRecoveryV248 as run
from fasRecoveryFieldsV248 import state,residual
from runSU2RocketPhaseV222 import ROOT,read,save,sha

BASE=run.BASE;CASE=run.CASE;CLOSED=False


def diagnose():
    out=BASE/'diagnosis01';out.mkdir(exist_ok=False)
    audit=read(BASE/'audit02/report.json');assert audit['passed'] and audit['developmentSucceeded']
    a=state(CASE/'state_000200.bin');a=a[np.argsort(a[:,0])]
    raw,r=residual(CASE/'terminal_spatial.bin');r=r[np.argsort(raw[:,0])]
    with (CASE/'rejections.csv').open(encoding='utf8') as stream:events=list(csv.DictReader(stream))
    ids=np.unique([int(e['globalID']) for e in events if int(e['level'])==0]);assert len(ids)==35
    total=np.sum(r*r,axis=0);share=np.sum(r[ids]*r[ids],axis=0)/total
    groups={}
    for e in events:
        key=e['phase'];groups[key]=groups.get(key,0)+1
    top=np.argsort(abs(r[:,4]))[-20:][::-1]
    details=[dict(globalID=int(i),xyzM=a[i,1:4].tolist(),wallDistanceM=float(a[i,24]),
      physicalEnergyRHS=float(r[i,4]),fractionOfEnergyRHSSquared=float(r[i,4]**2/total[4]),
      hadRecordedRecovery=bool(i in ids)) for i in top]
    # 距离分箱只描述位置，不替代原始区域责任或边界层外缘判断。
    boundaries=[0,1e-4,1e-3,1e-2,float('inf')];bins=[]
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        mask=(a[:,24]>=lo)&(a[:,24]<hi)
        bins.append(dict(lowerWallDistanceM=lo,upperWallDistanceM=hi if np.isfinite(hi) else None,
          points=int(sum(mask)),sixSquaredRHSFraction=(np.sum(r[mask]*r[mask],axis=0)/total).tolist()))
    assert sum(b['points'] for b in bins)==len(a)
    report=dict(passed=True,sourceAuditSHA256=sha(BASE/'audit02/report.json'),sourceResidualSHA256=sha(CASE/'terminal_spatial.bin'),
      completedUpdate=200,allPoints=len(a),recoveryPoints=len(ids),recoveryEventsByPhase=groups,
      sixEquations=['mass','momentumX','momentumY','momentumZ','physicalEnergy','SA'],
      recoveryPointsSquaredRHSFraction=share.tolist(),topEnergyRHSPoints=details,wallDistanceBins=bins,
      correlationNotCausation=True,recoveryDisturbanceNotBoundaryFlux=True,newCFDUpdates=0,physicalQualification=False)
    save(out/'report.json',report)
    print(dict(passed=True,recoveryPointsSquaredRHSFraction=share.tolist(),topEnergyRHSPoint=details[0]))


def main():
    if CLOSED:raise SystemExit('v248 diagnosis closed; follow NEXT_DECISION')
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/fasRecoveryFieldsV248.py',run.DOC/'PLAN.md',
      BASE/'audit02/report.json',CASE/'state_000200.bin',CASE/'terminal_spatial.bin',CASE/'rejections.csv']
    return run.setup_guard().supervise('diagnosis',[sys.executable,str(Path(__file__).resolve()),'task'],30,1536,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v248 diagnosis closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':diagnose()
    else:raise SystemExit(0 if main() else 1)
