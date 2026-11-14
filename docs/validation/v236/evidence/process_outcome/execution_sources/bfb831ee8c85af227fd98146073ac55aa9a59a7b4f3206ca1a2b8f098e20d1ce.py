"""模块职责：只读复用真实事件，连接旧admissibility记录、物理状态修正及原生提前停止保存差异。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runRestartBranchV236 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame
BASE=run.BASE;CASE=run.CASE;DOC=run.DOC


def rows(p):
    with p.open(encoding='utf8') as f:return list(csv.DictReader(f))


def main():
    out=BASE/'outcome01';out.mkdir(exist_ok=False)
    audit=read(BASE/'audit03/report.json');assert audit['passed']
    events=rows(CASE/'rejections.csv');old=rows(CASE/'admissibility.csv');assert len(events)==len(old)==466
    raw=read_frame(CASE/'state_000000.bin',335726,True)['values'];a=raw[np.argsort(raw[:,0])]
    counts={};windows={};max_energy_error=0;corrections=[]
    names=['rho','rhou','rhov','rhow','rhoE'];source_e_min=43997.3628006
    for r,s in zip(events,old):
        assert int(r['localID'])==int(s['local_point']) and int(r['globalID'])==int(raw[int(r['localID']),0])
        assert int(s['dimension'])==3 and int(s['recovered'])==int(r['recovered'])==1 and float(s['k'])==0 and s['reason']==r['reason']
        arrays=[]
        for key in ('trial','old'):
            x=np.array([float(r[key+str(j)]) for j in range(5)]);y=np.array([float(s[key+'_'+name]) for name in names])
            assert np.array_equal(x.view(np.uint64),y.view(np.uint64))
            e=x[4]/x[0]-.5*np.sum((x[1:4]/x[0])**2)
            max_energy_error=max(max_energy_error,abs(e-float(s[key+'_e'])))
            arrays.append(x)
        gid=int(r['globalID']);step=int(r['update']);delta=(arrays[1]-arrays[0])*a[gid,4];delta[4]-=300000*delta[0]
        counts[gid]=counts.get(gid,0)+1;win=(step-1)//50;windows[win]=windows.get(win,0)+1
        corrections.append([step,gid]+delta.tolist())
    assert max_energy_error<2e-9
    correction=np.array(corrections);unique=sorted(counts)
    details=[dict(globalID=i,events=counts[i],xyz=a[i,1:4].tolist(),volumeM3=float(a[i,4]),wallDistanceM=float(a[i,24])) for i in unique]
    np.savetxt(str(out/'observed_state_corrections.csv'),correction,delimiter=',',comments='',
        header='branchUpdate,globalID,deltaMassKg,deltaMomentumX,deltaMomentumY,deltaMomentumZ,deltaPhysicalEnergyJ',fmt='%.17g')
    save(out/'rejected_nodes.json',details)
    # 模块一：被强壁替换的残差行与补足闭合所需边界项符号相反；前审核的同向直接差不当作反力误差。
    ledger=rows(CASE/'wall_stages.csv');oriented=[]
    for r in ledger:
        if r['stage']!='before_strong':continue
        force=np.array([float(r['pNout'+d])-float(r['tauNout'+d]) for d in ('X','Y','Z')])
        eliminated=np.array([float(r['wallRow'+d]) for d in ('X','Y','Z')])
        oriented.append(dict(update=int(r['update']),eliminatedSpatialMomentumRow=eliminated.tolist(),outwardPressureMinusShear=force.tolist(),
            requiredBoundaryRow=(-eliminated).tolist(),oppositeSignClosureDifference=float(np.linalg.norm(eliminated+force)/max(np.linalg.norm(force),1)),
            independentGlobalMomentumBalanceProved=False))
    save(out/'oriented_wall_momentum.json',oriented)
    # 模块二：准确源码说明StopCalc未传给SetResultFiles，完整末帧已存但原生CSV计划未实现。
    source=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/iteration/CIteration.cpp'
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json');assert sha(source)==tree[str(source)]
    text=source.read_text(encoding='utf8');start=text.index('void CIteration::Output(');end=text.index('\n}',start)+2
    snippet=text[start:end];assert 'InnerIter);' in snippet and snippet.count('StopCalc')==1
    (out/'CIteration_Output.txt').write_text(snippet+'\n',encoding='utf8')
    absent=not list(CASE.glob('restart*.csv'));assert absent
    record=dict(passed=True,eventsMatchedToOriginalAdmissibility=len(events),uniqueRejectedNodes=len(unique),eventsBy50UpdateWindow=windows,
        eventInternalEnergyMaxAbsoluteDifferenceJkg=max_energy_error,allFiveStateOldTrialRecordsBitEqual=True,
        rejectedNodeCoordinateMin=a[unique,1:4].min(axis=0).tolist(),rejectedNodeCoordinateMax=a[unique,1:4].max(axis=0).tolist(),
        observedCumulativeStateCorrection=correction[:,2:].sum(axis=0).tolist(),
        correctionMeaning='volume times returned-minus-trial at each observed recovery; event sum, not time-integrated boundary flux or watts',
        wallMomentumSignSupplement=True,firstAuditSameDirectionDifferenceNotPhysicalReactionError=True,
        finalOppositeSignMomentumClosureDifference=oriented[-1]['oppositeSignClosureDifference'],
        nativeFinalRestartMissing=absent,plannedNativeRestartAcceptance=False,completeBinaryFinalFramePreserved=True,
        nativeSourceSHA256=sha(source),nativeOutputStopArgumentUnused=True,
        futureSaveFix='after native StartSolver returns, call existing COutput::SetResultFiles with force_writing=true; do not rerun this branch to fill a CSV',
        physicalQualification=False,fullConservationProved=False,uniqueRootCauseProved=False)
    save(out/'report.json',record);print(json.dumps(record))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.configure_guard();locks=[Path(__file__).resolve(),ROOT/'scripts/runRestartBranchV236.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            ROOT/'scripts/singleDomainFrame.py',BASE/'audit03/report.json',CASE/'rejections.csv',CASE/'admissibility.csv',CASE/'state_000000.bin',CASE/'wall_stages.csv',
            ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/iteration/CIteration.cpp',ROOT/'outputs/v167_validation/build01/source_tree.json']
        raise SystemExit(0 if run.guard.supervise('outcome',[sys.executable,str(Path(__file__).resolve()),'task'],30,1024,locks,BASE,'audit') else 1)
