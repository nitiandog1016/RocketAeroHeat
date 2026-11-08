"""准备温变输运与独立Couette热壁参照；不启动原生求解器。"""
from pathlib import Path
import json,time,shutil
import numpy as np
from scipy.integrate import solve_bvp
from prepareOpenFoamControl import ROOT,header,save,sha
from phaseAwareViscousAir import PhaseAwareViscousAir,TRANSPORT

BASE=ROOT/'outputs/v207_validation'
def field(name,kind,dimensions,values,boundaries):
    vals=['('+' '.join('%.17g'%x for x in row)+')' for row in values] if kind=='vector' else ['%.17g'%x for x in values]
    return header(name,'volVectorField' if kind=='vector' else 'volScalarField')+'dimensions '+dimensions+';\ninternalField nonuniform List<'+kind+'> '+str(len(vals))+'\n(\n'+'\n'.join(vals)+'\n);\nboundaryField {\n'+boundaries+'\n}\n'

def main():
    start=time.monotonic();out=BASE/'prepare01';out.mkdir(parents=True,exist_ok=False);gas=PhaseAwareViscousAir()
    if not json.loads((BASE/'source01/report.json').read_text())['passed']:raise ValueError('来源未通过')
    names=['scripts/phaseAwareViscousAir.py','scripts/phaseAwareDryAir.py','scripts/prepareViscousThermalControl.py',
        'docs/validation/v207/THERMAL_CONTROL_PLAN.md','data/thermo/phase_aware_dry_air_v1.json']
    before={n:sha(ROOT/n) for n in names}
    for n in names:
        target=out/'executed_source'/n;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/n,target)
    # 模块一：验证新输运插值，不重跑v206真实EOS或旧量热查询。
    mid=(gas.T[1:]+gas.T[:-1])/2;reference=gas.transport(1,mid);interp=gas.interpolated_transport(1,mid)
    errors={key:float(np.max(abs(interp[key]/reference[key]-1))) for key in interp}
    if max(errors.values())>1e-4:raise ValueError('新输运插值误差超标')
    # 模块二：独立成熟边值法解变量输运Couette，u/T缩放避免绝对量级偏差。
    gap=.001;speed=100.;p=100000.;z=np.linspace(0,1,33)
    def ode(z,state,parameters):
        tau,C=parameters;T=300*state[1];props=gas.transport(p,T)
        return np.vstack((np.full_like(z,tau*gap/speed)/props['dynamicViscosityPas'],
            (C-tau*speed*state[0])*gap/(300*props['conductivityWmK'])))
    def bc(left,right,parameters):return np.array([left[0],right[0]-1,left[1]-1,right[1]-2])
    props=gas.transport(p,450);tau0=float(props['dynamicViscosityPas'])*speed/gap
    oracle=solve_bvp(ode,bc,z,np.vstack((z,1+z)),p=[tau0,float(props['conductivityWmK'])*300/gap+tau0*speed/2],tol=1e-9,bc_tol=1e-11,max_nodes=2000)
    if not oracle.success:raise ValueError('独立边值参照未收敛 '+oracle.message)
    y=(np.arange(32)+.5)*gap/32;state=oracle.sol(y/gap);U=np.zeros((32,3));U[:,0]=state[0]*speed;T=state[1]*300
    tau,C=map(float,oracle.p);walls=dict(cold=C,hot=-C+tau*speed)
    boundary_error=float(abs(bc(oracle.sol(0),oracle.sol(1),oracle.p)).max())
    power_error=abs(sum(walls.values())-tau*speed)/max(1,abs(tau*speed))
    if boundary_error>1e-8 or power_error>1e-8:raise ValueError('参照边界或黏性功账未过')
    case=BASE/'thermal01';case.mkdir(exist_ok=False)
    for name in ['0','constant','system']:(case/name).mkdir()
    (case/'constant/physicalProperties').write_text(header('physicalProperties')+'thermoType { type hePsiThermo; mixture pureMixture; transport icoTabulated; thermo hIcoTabulated; equationOfState perfectGas; specie specie; energy rahPhaseAwareInternalEnergy; }\nmixture {\n'+gas.dictionary()+'}\n')
    (case/'constant/momentumTransport').write_text(header('momentumTransport')+'simulationType laminar;\n')
    (case/'constant/thermophysicalTransport').write_text(header('thermophysicalTransport')+'simulationType laminar; laminar { model Fourier; }\n')
    cyclic='\n'.join(n+' { type cyclic; }' for n in ['xMinus','xPlus','zMinus','zPlus'])
    (case/'0/p').write_text(field('p','scalar','[1 -1 -2 0 0 0 0]',np.full(32,p),'cold { type zeroGradient; } hot { type zeroGradient; }\n'+cyclic))
    (case/'0/T').write_text(field('T','scalar','[0 0 0 1 0 0 0]',T,'cold { type fixedValue; value uniform 300; } hot { type fixedValue; value uniform 600; }\n'+cyclic))
    (case/'0/U').write_text(field('U','vector','[0 1 -1 0 0 0 0]',U,'cold { type fixedValue; value uniform (0 0 0); } hot { type fixedValue; value uniform (100 0 0); }\n'+cyclic))
    (case/'system/blockMeshDict').write_text(header('blockMeshDict')+'''convertToMeters 1;
vertices ((0 0 0) (.02 0 0) (.02 .001 0) (0 .001 0) (0 0 .02) (.02 0 .02) (.02 .001 .02) (0 .001 .02));
blocks (hex (0 1 2 3 4 5 6 7) (1 32 1) simpleGrading (1 1 1));
edges ();
boundary (
cold { type wall; faces ((0 1 5 4)); }
hot { type wall; faces ((3 7 6 2)); }
xMinus { type cyclic; neighbourPatch xPlus; faces ((0 4 7 3)); }
xPlus { type cyclic; neighbourPatch xMinus; faces ((1 2 6 5)); }
zMinus { type cyclic; neighbourPatch zPlus; faces ((0 3 2 1)); }
zPlus { type cyclic; neighbourPatch zMinus; faces ((4 5 6 7)); }
);
''')
    (case/'system/fvSchemes').write_text(header('fvSchemes')+'''fluxScheme Kurganov;
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; reconstruct(rho) upwind; reconstruct(U) upwind; reconstruct(T) upwind; }
snGradSchemes { default corrected; }
''')
    (case/'system/fvSolution').write_text(header('fvSolution')+'''solvers {
"rho.*" { solver diagonal; }
"(U|e).*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0; }
}
PIMPLE {}
''')
    # 原生库实际SHA在构建后/运行前绑定，这里只给固定目标路径。
    native='/mnt/e/RocketAeroHeat/outputs'
    (case/'system/controlDict').write_text(header('controlDict')+f'''solver shockFluid;
libs ("{native}/v206_validation/native01/library/librahPhaseAwareAir.so" "{native}/v207_validation/native01/library/librahPhaseAwareViscousAir.so");
startFrom startTime; startTime 0; stopAt endTime; endTime 4e-8; deltaT 5e-9;
writeControl timeStep; writeInterval 2; purgeWrite 0; writeFormat ascii; writePrecision 17;
timeFormat general; timePrecision 12; writeCompression off; runTimeModifiable false; adjustTimeStep no;
functions {{
phaseGuard {{ type rahPhaseStateGuard; model rah-phase-aware-dry-air-of13-v1; executeControl timeStep; executeInterval 1; executeAtStart true; writeControl none; }}
energy {{ type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e phi); writeControl writeTime; }}
heat {{ type wallHeatFlux; libs ("libfieldFunctionObjects.so"); patches (cold hot); executeControl timeStep; executeInterval 1; writeControl writeTime; }}
}}
''')
    save(case/'thermo_contract.json',gas.viscous_contract())
    save(out/'oracle.json',dict(method='SciPy solve_bvp; variable-mu/kappa exact planar steady balances',yM=y.tolist(),temperatureK=T.tolist(),velocityMps=U.tolist(),pressurePa=p,
        wallHeatIntoWallWm2=walls,wallShearPa=tau,viscousPowerPerAreaWm2=tau*speed,wallAreaM2=.0004,
        boundaryMaxError=boundary_error,powerMaxScaledError=power_error,solverResidualMax=float(oracle.rms_residuals.max()),
        constantEnergyFluxC=C,gridNodes=int(oracle.x.size),isExperiment=False,physicalQualification=False))
    stable={n:sha(ROOT/n)==h for n,h in before.items()}
    size=sum(x.stat().st_size for d in [out,case] for x in d.rglob('*') if x.is_file())
    if not all(stable.values()) or time.monotonic()-start>60 or size>16*1024**2:raise ValueError('准备预算/输入改变')
    report=dict(passed=True,transport=TRANSPORT,interpolationMidpoints=len(mid),interpolationRelativeErrors=errors,
        inputSHA256=before,inputStable=stable,caseInputSHA256={str(x.relative_to(ROOT)):sha(x) for x in case.rglob('*') if x.is_file()},
        referenceBoundaryError=boundary_error,referencePowerError=power_error,physicalQualification=False,
        CFDRuns=0,builds=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(dict(passed=True,interpolationErrors=errors,seconds=report['elapsedSeconds'],wallHeat=walls)))

if __name__=='__main__':main()
