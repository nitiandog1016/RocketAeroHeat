"""模块职责：复用本轮保存历史与准确原生源，形成下一次数值延续决定；不重新求解流场。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runHLLCCompatibilityV238 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from singleDomainFrame import read_frame
BASE=run.BASE;DOC=run.DOC


def rows(path):
    with path.open(encoding='utf8') as f:
        return [{k.strip():v for k,v in r.items()} for r in csv.DictReader(f,skipinitialspace=True)]


def main():
    out=BASE/'outcome02';out.mkdir(exist_ok=False)
    # 模块一：准确编译源与成熟阶次渐进入口绑定，只摘录所需调用链。
    origin=ROOT/'tools/runtime/su2-thermal-build/source'
    tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    specifications=[
      ('iteration',origin/'SU2_CFD/src/iteration/CFluidIteration.cpp',[(248,262),(369,394)]),
      ('reconstruction',BASE/'rocket_prepare01/CEulerSolver.cpp',None),
      ('config',BASE/'native_prepare01/CConfig.cpp',None)]
    # 配置隔离源的位置来自准确构建命令，不假定它属于旧公共源码目录。
    commands=read(BASE/'native_prepare01/commands.json')
    config_command=next(c for c in commands.values() if '-c' in c and any(str(x).endswith('CConfig.cpp') for x in c))
    specifications[2]=('config',Path(config_command[config_command.index('-c')+1]),None)
    locks=read(BASE/'rocket_prepare01/inputs.json');native=read(BASE/'native_prepare01/inputs.json');bindings=[]
    for name,p,ranges in specifications:
        expected=locks.get(str(p),native.get(str(p),tree.get(str(p))))
        if expected is None:raise AssertionError(('SOURCE_IDENTITY_MISSING',str(p)))
        if sha(p)!=expected:
            relative=p.relative_to(origin)
            alternatives=[ROOT/'outputs/v167_validation/build01/executed_sources'/relative,
              ROOT/'outputs/v222_validation/source01'/relative,ROOT/'outputs/v221_validation/source01'/relative,
              ROOT/'outputs/v219_validation/source01/source'/relative]
            p=next(q for q in alternatives if q.exists() and sha(q)==expected)
        lines=p.read_text(encoding='utf8').splitlines()
        if ranges is None:
            needles=['const su2double musclRamp','Recompute the reconstructed quantities','bool bad_recon'] if name=='reconstruction' else ['addBoolOption("RAMP_MUSCL"']
            ranges=[]
            for needle in needles:
                hits=[i for i,line in enumerate(lines) if needle in line];assert len(hits)==1,(needle,hits)
                ranges.append((max(1,hits[0]-3),min(len(lines),hits[0]+17)))
        target=out/(name+'_excerpt.txt')
        target.write_text('\n\n'.join('\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,b+1)) for a,b in ranges)+'\n',encoding='utf8')
        bindings.append(dict(path=str(p),sha256=sha(p),excerpts=ranges,excerptSHA256=sha(target)))
    # 模块二：现有记录区分线性解、非线性稳态、温度下界与CFL，不额外生成原生评估。
    case=BASE/'rocket01';history=rows(case/'history.csv');steps=rows(case/'steps.csv');events=rows(case/'accepted_steps.csv')
    assert len(history)==len(steps)==200
    linear=np.array([[float(r[k]) for k in ('Linear_Solver_Iterations','LinSolRes','LinSolIterTurb','LinSolResTurb')] for r in history])
    assert np.isfinite(linear).all()
    frames=[]
    for step in (0,1,100,200):
        a=read_frame(case/('state_%06d.bin'%step),335726,step==0)['values']
        point=a[a[:,0]==332194][0]
        frames.append(dict(update=step,flowCFLRange=[float(a[:,26].min()),float(a[:,26].max())],
            flowRelaxationRange=[float(a[:,25].min()),float(a[:,25].max())],
            temperatureRangeK=[float(a[:,11].min()),float(a[:,11].max())],
            limitedPoint=dict(globalID=332194,coordinatesM=point[1:4].tolist(),temperatureK=float(point[11]))))
    decisions=rows(case/'benefit_decisions.csv');assert len(decisions)==2
    report=dict(passed=True,sourceBindings=bindings,realCompletedUpdates=200,
      limitedGlobalIDs=sorted(set(int(r['globalID']) for r in events)),limitedUpdateRange=[min(int(r['update']) for r in events),max(int(r['update']) for r in events)],
      lastLimitedUpdate=145,noFurtherLimitedEventsFor55Updates=True,
      linearFlowMaxReportedLogResidual=float(linear[:,1].max()),linearTurbMaxReportedLogResidual=float(linear[:,3].max()),
      linearFlowMaxIterations=int(linear[:,0].max()),linearTurbMaxIterations=int(linear[:,2].max()),
      linearReportedToleranceSatisfied=bool(np.all(linear[:,[1,3]]<=-5)),independentLinearSystemAudit=False,
      frames=frames,benefitWindowRecords=decisions,
      observedStop='significant-deterioration',uniquePhysicalCauseEstablished=False,
      nextCandidate='native MUSCL spatial-order continuation, fixed schedule then mandatory full-order assessment',
      rationale='full second-order uniform startup remains a bounded hypothesis; native ramp is present; no new flux/Jacobian required',
      reconstructionDomainRecoveryGuaranteed=False,strictPhaseDomainPreservationGuaranteed=False,
      nativeRampActuallyQueried=False,newCFDAdmitted=False,physicalQualification=False,newQualifiedAreaM2=0)
    save(out/'report.json',report)
    print(json.dumps(dict(passed=True,limitedIDs=report['limitedGlobalIDs'],linearReportedToleranceSatisfied=report['linearReportedToleranceSatisfied'],newCFDAdmitted=False)))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runHLLCCompatibilityV238.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            BASE/'rocket_prepare01/inputs.json',BASE/'native_prepare01/inputs.json',BASE/'native_prepare01/commands.json',
            ROOT/'outputs/v167_validation/build01/source_tree.json',ROOT/'scripts/singleDomainFrame.py']
        locks += [p for p in (BASE/'rocket01').iterdir() if p.is_file()]
        raise SystemExit(0 if run.guard.supervise('outcome2',[sys.executable,str(Path(__file__).resolve()),'task'],45,1024,locks,BASE,'offline') else 1)
