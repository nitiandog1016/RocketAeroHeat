"""模块职责：独立审核冻结场的内部消去、壁热和强行替换账；不把代数一致授为稳态资格。"""
from pathlib import Path
import csv,json,sys
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from su2PhaseAirFields import inverse_energy
from phaseAwareViscousAir import PhaseAwareViscousAir
BASE=ROOT/'outputs/v223_validation';DOC=ROOT/'docs/validation/v223';CASE=BASE/'frozen02'

def scaled(a,b):
    a=np.asarray(a);b=np.asarray(b)
    return float(np.max(abs(a-b)/np.maximum(1,np.maximum(abs(a),abs(b)))))

def main():
    out=BASE/'audit01';out.mkdir(exist_ok=False)
    native=read(CASE/'native.json');run=read(BASE/'process_query2/report.json')
    assert run['passed'] and run['inputsStable'] and not run['liveGroupAfterCleanup']
    assert native['spaceEvaluations']==2 and native['flowUpdates']==0
    assert native['stagedNativeRHSMismatches']==native['currentStateBitChanges']==0
    stages={}
    with (CASE/'stages.csv').open(encoding='utf8') as f:
        reader=csv.reader(f);next(reader)
        for row in reader:stages[row[0]]=np.array([float(v) for v in row[1:]])
    field=np.loadtxt(str(CASE/'field.csv'),delimiter=',',skiprows=1)
    wall=np.loadtxt(str(CASE/'wall.csv'),delimiter=',',skiprows=1)
    field=field[np.argsort(field[:,0])];wall=wall[np.argsort(wall[:,0])]
    assert field.shape==(335726,13) and wall.shape==(9878,21)
    assert np.isfinite(field).all() and np.isfinite(wall).all()
    assert np.array_equal(field[:,0],np.arange(len(field))) and np.array_equal(wall[:,0],np.arange(len(wall)))
    assert (field[:,1]>0).all() and (wall[:,1]>0).all() and (wall[:,7]>0).all()
    # 模块一：准确检查点与当前守恒状态。保留符号位，不以allclose掩盖状态修改。
    original=ROOT/'outputs/v222_validation/rocket01/restart_flow_000199.csv'
    expected=read(ROOT/'outputs/v222_validation/audit01/report.json')['outputSHA256'][original.name]
    assert sha(original)==expected
    saved=np.loadtxt(str(original),delimiter=',',skiprows=1,usecols=range(4,10))
    assert np.array_equal(field[:,2:8].copy().view(np.uint64),saved.copy().view(np.uint64))
    air=PhaseAwareViscousAir();rho=field[:,2]
    velocity=field[:,3:6]/rho[:,None]
    temp=inverse_energy(air,field[:,6]/rho-.5*np.sum(velocity**2,axis=1)-300000)
    pressure=rho*air.R*temp
    air.evaluate(pressure,temp)
    ids=wall[:,0].astype(int);neighbors=wall[:,6].astype(int)
    assert np.array_equal(field[ids,3:6],np.zeros((len(ids),3))) and np.all(field[ids,7]==0)
    thermo=max(scaled(temp[ids],wall[:,8]),scaled(pressure[ids],wall[:,9]),
               scaled(rho[ids],wall[:,10]),scaled(temp[neighbors],wall[:,12]))
    k=air.interpolated_transport(pressure[ids],temp[ids])['conductivityWmK']
    conductivity=scaled(k,wall[:,11]);assert thermo<2e-11 and conductivity<2e-11
    # 模块二：对流入口已融合边黏性；各内部阶段的有符号五量应成对消去到舍入量级。
    internal={name:(abs(stages[name][:5])/np.maximum(1,stages[name][5:])).tolist()
              for name in ('convective','viscous','source','fluid_interface')}
    assert max(max(v) for v in internal.values())<2e-12
    final=stages['native_dispatch'][:5];far=stages['farfield'][:5]-stages['fluid_interface'][:5]
    change=stages['isothermal'][:5]-stages['fixed_values'][:5]
    assert np.array_equal(stages['native_dispatch'],stages['isothermal'])
    assert change[0]==0
    # 模块三：壁热用原生法向邻点+独立线性k，质量不变时能量坐标平移不改变边界热量。
    q=k*(temp[neighbors]-300)/wall[:,7];power=float(np.sum(wall[:,1]*q))
    heat_error=scaled(q,wall[:,17]);power_error=scaled(power,change[4])
    removed=-np.sum(wall[:,13:16],axis=0);replacement_error=scaled(removed,change[1:4])
    assert max(heat_error,power_error,replacement_error)<2e-11
    rhs=field[:,8:13];sum_error=scaled(np.sum(rhs,axis=0),final)
    assert sum_error<2e-11
    physical=rhs.copy();physical[:,4]-=300000*physical[:,0]
    transform=lambda r:np.r_[r[:4],r[4]-300000*r[0]]
    # 原生完整应力牵引来自Pressure/Friction_Forces；Cf只含切向分量，不能伪称独立完整重建。
    force=np.array(native['bodyForceCoefficients'])*native['forceRefN']
    force_gap=removed-force
    volume=field[:,1]
    inventory=np.sum(field[:,2:7]*volume[:,None],axis=0);inventory[4]-=300000*inventory[0]
    report=dict(passed=True,model='rah-su2-phase-air-saneg-eplus300k-v1',points=len(field),wallNodes=len(wall),
        flowUpdates=0,spaceEvaluations=2,currentStateBitChanges=0,stagedNativeRHSMismatches=0,
        currentCheckpointSixVariablesBitwiseEqual=True,sourceCheckpointSHA256=expected,
        inputSHA256={p.name:sha(p) for p in (CASE/'field.csv',CASE/'wall.csv',CASE/'stages.csv',CASE/'native.json')},
        thermalScaledDifference=thermo,conductivityScaledDifference=conductivity,internalCancellationScaled=internal,
        internalConvectionIncludesViscosity=True,wallHeatFluxScaledDifference=heat_error,
        wallHeatPowerW=power,wallEnergyStageDifferenceW=float(change[4]),wallEnergyScaledDifference=power_error,
        wallMomentumReplacementN=removed.tolist(),wallMomentumReplacementScaledDifference=replacement_error,
        nativePhysicalBodyForceN=force.tolist(),reactionMinusNativePhysicalBodyForceN=force_gap.tolist(),
        reactionForceRelativeNorm=float(np.linalg.norm(force_gap)/np.linalg.norm(force)),
        nativeBodyForceIndependentlyReconstructed=False,
        farfieldPhysicalFlux=transform(far).tolist(),finalPhysicalResidualSum=transform(final).tolist(),
        physicalResidualRMS=np.sqrt(np.mean(physical**2,axis=0)).tolist(),
        domainPhysicalConservativeInventory=inventory.tolist(),nodeRhsSumScaledDifference=sum_error,
        wallFluidTemperatureRangeK=[float(min(temp[ids])),float(max(temp[ids]))],
        frozenHeatFluxRangeWm2=[float(min(q)),float(max(q))],
        frozenCurrentFlowOperator=True,fullFlowSAIterationReproduced=False,localTimeStepOrRejectionTraceObserved=False,
        thermalSteadyBalancePassed=False,fullPhysicalConservationQualified=False,
        physicalQualification=False,newQualifiedAreaM2=0,uniqueCauseEstablished=False)
    save(out/'report.json',report)
    print(json.dumps({k:report[k] for k in ('passed','wallHeatPowerW','reactionForceRelativeNorm','finalPhysicalResidualSum')}))

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
        actual=guard.size_bytes;guard.size_bytes=lambda unused:actual(BASE)
        locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md']
        locks += [CASE/n for n in ('native.json','field.csv','wall.csv','stages.csv')]
        raise SystemExit(0 if guard.supervise('audit',[sys.executable,str(Path(__file__).resolve()),'task'],80,256,locks,BASE,'offline') else 1)
