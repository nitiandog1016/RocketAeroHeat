"""模块职责：核对已有整体控制、准确调用源和备用资源；不重新运行任何原生求解。"""
from pathlib import Path
import csv,json,sys
import numpy as np
from phaseAwareViscousAir import PhaseAwareViscousAir
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v234_validation';DOC=ROOT/'docs/validation/v234'
CLOSED=False


def rows(path):
    with path.open(encoding='utf8') as f:return list(csv.DictReader(f))


def audit():
    out=BASE/'audit01';out.mkdir(exist_ok=False);errors={};air=PhaseAwareViscousAir()
    def check(name,a,b,tol=2e-9):
        a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
        assert np.isfinite(a).all() and np.isfinite(b).all(),name
        e=float(np.max(np.abs(a-b)/np.maximum(1.,np.abs(b))))
        errors[name]=max(errors.get(name,0.),e);assert e<=tol,(name,e)
    def nums(r,keys):return np.array([float(r[k]) for k in keys])
    def key(r):return (int(r['scenario']),)+tuple(float(r[k]) for k in ('x','y','z'))
    # 模块一：单区已有四场的物性/能量/SA身份、分阶段热壁与接受；不重复旧物性全表。
    case=BASE/'single02';native=read(case/'native.json')
    assert native==dict(zones=1,scenarios=4,flowUpdates=0,currentStateBitChanges=0,actualTransferDataCalled=False)
    nr=rows(case/'nodes.csv');assert len(nr)==96
    nodes={key(r):r for r in nr};assert len(nodes)==96
    qr=air.evaluate(1200.,350.123);a=float(qr['soundSpeedMps'])
    nu=3*float(air.interpolated_transport(1200.,350.123)['dynamicViscosityPas'])/float(qr['densityKgm3'])
    for r in nr:
        s=int(r['scenario']);x=nums(r,['x','y','z']);t=350.123+(x@np.array([.25,1.5,-.75]) if s else 0)
        u=a*np.array([2.,.4,.3])+(np.array([4*x[1],3*x[2],2*x[0]]) if s else 0)
        q=air.evaluate(1200.,t);tr=air.interpolated_transport(1200.,t);rho=float(q['densityKgm3'])
        E=float(q['internalEnergyJkg'])+300000+.5*(u@u)
        check('fiveNativeState',nums(r,['U%d'%i for i in range(5)]),np.r_[rho,rho*u,rho*E])
        check('physicalTotalH',float(r['P6'])-300000,float(q['enthalpyJkg'])+.5*(u@u))
        check('primitive',nums(r,['P%d'%i for i in range(8)]),np.r_[t,u,1200.,rho,E+1200./rho,float(q['soundSpeedMps'])])
        nt=(-1 if s==2 else 1)*nu*(1+.05*x[0] if s else 1)
        check('SA',float(r['SA'])/nu,nt/nu)
        mu=float(tr['dynamicViscosityPas']);chi=nt/(mu/rho)
        mut=0 if nt<0 else rho*nt*chi**3/(chi**3+7.1**3)
        check('transportRatios',[float(r['P8'])/mu,float(r['P10'])/float(tr['conductivityWmK'])],[1.,1.])
        check('muT',float(r['P9']),mut)
        if s<3:
            G=np.zeros((7,3))
            if s:
                G[0]=[.25,1.5,-.75];G[1:4]=[[0,4,0],[0,0,3],[2,0,0]]
                G[5]=-rho/t*G[0];G[6]=float(q['cpJkgK'])*G[0]+u@G[1:4]
            check('prescribedGradients',nums(r,['G%d%d'%(j,d) for j in range(7) for d in range(3)]),G.ravel())
        assert float(r['wallDistance'])>=0 and float(r['volume'])>0
    stages=rows(case/'stages.csv');assert len(stages)==288
    stage_map={(key(r),r['stage']):r for r in stages};thermal=[];reactions=[]
    for r in rows(case/'walls.csv'):
        n=nodes[key(r)];A=float(r['area']);d=float(r['distance']);delta=float(r['Tnear'])-300
        k=float(r['k']);keff=k+float(n['P11'])*float(n['P9'])/.9
        check('wallEnergyIncrement',float(r['REwall']),keff*delta/d*A)
        check('nativeLaminarWallQ',float(r['q']),k*delta/d)
        check('strongOldWall',nums(r,['oldMx','oldMy','oldMz','oldSA','postRx','postRy','postRz']),np.zeros(7))
        before=stage_map[(key(r),'before_wall')];after=stage_map[(key(r),'after_wall')]
        check('sameStageWallEnergy',(float(after['R4'])-300000*float(after['R0']))-(float(before['R4'])-300000*float(before['R0'])),float(r['REwall']))
        check('sameStageMomentum',nums(before,['R1','R2','R3']),nums(r,['preRx','preRy','preRz']))
        check('wallDensityUnchanged',float(after['R0']),float(before['R0']))
        thermal.append(dict(scenario=int(r['scenario']),xyz=nums(r,['x','y','z']).tolist(),
          energyIncrement=float(r['REwall']),nativeLaminarQArea=float(r['q'])*A,
          turbulentHeatContribution=(keff-k)*delta/d*A))
        reactions.append(dict(scenario=int(r['scenario']),xyz=nums(r,['x','y','z']).tolist(),
          projectedMomentumRow=nums(r,['preRx','preRy','preRz']).tolist(),
          pressureArea=(float(n['P4'])*nums(r,['nx','ny','nz'])).tolist(),steadyBalanceClaim=False))
    for r in rows(case/'recovery.csv'):
        assert int(r['SAUnchangedBits'])==1
        if int(r['trial'])==0:assert [int(r[k]) for k in ('returned','threw','flowRestoredBits')]==[0,0,5]
        else:assert int(r['threw'])==1 and int(r['returned'])==-1
    assert len(rows(case/'recovery.csv'))==2
    save(out/'single_wall_energy.json',thermal);save(out/'same_stage_reaction.json',reactions)
    # 模块二：失败定位与准确源码绑定；缺失证据明确为null，不把partial CSV当完整双区结果。
    log=(BASE/'process_capture4/stdout.txt').read_text(encoding='utf8')
    assert 'SIGSEGV' in log and 'CGeneralAvgGrad_Flow::ComputeResidual' in log and 'BC_Fluid_Interface' in log
    assert '"interface"' in log and '$2 = 0' in log and '$3 = 0' in log
    assert read(BASE/'process_multi4/report.json')['exitCode']==3221225477
    old=ROOT/'tools/runtime/su2-thermal-build/source';tree=read(ROOT/'outputs/v167_validation/build01/source_tree.json')
    known=read(ROOT/'outputs/v220_validation/build02/dependencies.json');bindings=[]
    specs=[('SU2_CFD/include/solvers/CFVMFlowSolverBase.inl','void CFVMFlowSolverBase<V, FlowRegime>::BC_Fluid_Interface',162),
      ('SU2_CFD/src/numerics/flow/flow_diffusion.cpp','CNumerics::ResidualType<> CGeneralAvgGrad_Flow::ComputeResidual',155),
      ('SU2_CFD/src/solvers/CNSSolver.cpp','su2double thermal_conductivity = nodes->GetThermalConductivity(iPoint) + Cp',65),
      ('SU2_CFD/include/solvers/CFVMFlowSolverBase.inl','/*--- Compute non-dimensional velocity and y+ ---*/',44)]
    snippets=[]
    for i,(rel,anchor,count) in enumerate(specs):
        p=ROOT/'outputs/v220_validation/source'/rel
        if not p.exists():p=old/rel
        expected=known.get(str(p),tree.get(str(old/rel)));assert expected
        if sha(p)!=expected:p=ROOT/'outputs/v167_validation/build01/executed_sources'/rel
        assert sha(p)==expected,str(p)
        lines=p.read_text(encoding='utf8').splitlines();starts=[i for i,line in enumerate(lines) if anchor in line];assert len(starts)==1
        first=starts[0];last=min(len(lines),first+count);excerpt='\n'.join('%6d %s'%(j+1,lines[j]) for j in range(first,last))+'\n'
        target=out/('%02d_%s.txt'%(i,Path(rel).name));target.write_text(excerpt,encoding='utf8');snippets.append(excerpt)
        bindings.append(dict(path=str(p),sha256=expected,firstLine=first+1,lastLine=last,excerptSHA256=sha(target)))
    assert 'visc_numerics->SetSecondary' not in snippets[0] and 'conv_numerics->SetSecondary' in snippets[0]
    assert 'S_i[iVar+2]+S_j[iVar+2]' in snippets[1]
    assert 'nodes->GetGradient_Primitive(iPoint)' in snippets[0]
    assert 'RAH_PHASE_AIR_MZ' not in snippets[2] and 'RAH_PHASE_AIR)' in snippets[2]
    save(out/'source_binding.json',bindings)
    original=read(ROOT/'docs/validation/v222/evidence/final_report.json')
    r=read(ROOT/'outputs/v222_validation/process_engine/report.json');assert r['passed'] and original['interfacePassed'] and not original['physicalQualification']
    save(out/'route_comparison.json',dict(selected='existing-single-connected-phase-Air-SA-neg',
      existingSingleDomain=dict(cells=705490,points=335726,wallFaces=19752,peakMeasuredRSSBytes=r['peakObservedGroupRSSBytes'],
        elapsedSeconds=r['elapsedSeconds'],completedUpdates=200,converged=False,resumeAdmitted=False,referenceReportSHA256=sha(ROOT/'outputs/v222_validation/process_engine/report.json')),
      dualNative=dict(eligible=False,matchingCoordinatesObserved=True,interfaceFluxBalance=None,completeWallDistanceComparison=None,
        reason='general-EOS viscous interface lacks secondary-state assignment; native SIGSEGV at first interface',
        adapterThermalJacobianDispatchIncomplete=True,donorGradientTransferVerified=False),
      fullEilmer=dict(admitted=False,knownMinimumAllocationBytes=8743123908,existingWSLLimitBytes=8*1024**3),
      noClaimSingleDomainFaster=True,noFullRocketExecutionAdmitted=True))
    save(out/'report.json',dict(passed=True,singleControlVerified=True,singleNodeStates=len(nr),singleWallRows=len(thermal),
      errorMetrics=errors,nativeStateBitChanges=0,flowUpdates=0,dualCompatibilityPassed=False,
      nativeCrashPhase='interface/scenario0/zone0',dualInterfacePairFluxVerified=False,dualThermalWallVerified=False,
      fullSASourceChainIncluded=False,wallManufacturedEddyViscosityExplainsFluxDifference=True,
      uniqueFaultingPointerProven=False,sourceMissingViscousSecondaryConfirmed=True,missingMZThermalJacobianDispatchConfirmed=True,
      productionDualRouteClosed=True,singleFallbackSelectedNotAdmitted=True,newQualifiedAreaM2=0,physicalQualification=False))
    print(json.dumps(dict(singleControlVerified=True,dualCompatibilityPassed=False,maxScaledError=max(errors.values()))))


def main():
    if CLOSED:raise SystemExit('v234 offline audit closed; follow NEXT_DECISION')
    guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC;size=guard.size_bytes;guard.size_bytes=lambda _:size(BASE)
    locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',ROOT/'scripts/phaseAwareViscousAir.py',ROOT/'scripts/phaseAwareDryAir.py',ROOT/'data/thermo/phase_aware_dry_air_v1.json']
    locks+=list((BASE/'single02').glob('*.csv'))+[BASE/'single02/native.json',BASE/'process_capture4/stdout.txt',BASE/'process_multi4/report.json']
    return guard.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],120,256,locks,BASE)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v234 offline audit closed; follow NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':audit()
    else:raise SystemExit(0 if main() else 1)
