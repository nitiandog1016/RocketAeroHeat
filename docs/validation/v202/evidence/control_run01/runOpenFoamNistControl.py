"""新具名NIST低温原生控制；只运行事前登记的一项，不自动续算。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha,header
from nistIntegratedAir import AIR,evaluate,contract,MODEL
from openFoamFields import audit_state
from openFoamSynchronousControl import PhysicalTimeGuard

BASE=ROOT/'outputs/v202_validation'


def prepare(case):
    """模块一：复用已过原网格/算子，仅在新目录设置新气体及低温工况。"""
    old=ROOT/'outputs/v199_validation/small01'
    for folder in ['constant','system','0']:shutil.copytree(old/folder,case/folder)
    for name in ['exact_mesh.npz','patches.json','wall_identity.json']:shutil.copyfile(old/name,case/name)
    for p in (old/'constant/polyMesh').iterdir():assert sha(p)==sha(case/'constant/polyMesh'/p.name)
    library=BASE/'native01/library/librahNistAir.so'
    native=json.loads((BASE/'native01/report.json').read_text(encoding='utf-8'))
    assert native['passed'] and sha(library)==native['registrationLibrarySHA256']
    (case/'constant/physicalProperties').write_text(header('physicalProperties')+
        'thermoType { type hePsiThermo; mixture pureMixture; transport const; thermo hIcoTabulated; equationOfState perfectGas; specie specie; energy sensibleInternalEnergy; }\nmixture {\n'+AIR.dictionary()+'}\n')
    save(case/'thermo_contract.json',contract())
    gas=evaluate(1000.,180.);velocity=-3*float(gas['soundSpeedMps'])
    values={'p':'1000','T':'180','U':'(%.17g 0 0)'%velocity}
    dims={'p':'1 -1 -2 0 0 0 0','T':'0 0 0 1 0 0 0','U':'0 1 -1 0 0 0 0'}
    patches=json.loads((case/'patches.json').read_text())
    for name,value in values.items():
        rows=[]
        for patch in patches:
            tag=patch['name']
            if tag=='inflow':bc='type fixedValue; value uniform '+value+';'
            elif tag=='outer':bc='type freestream; freestreamValue uniform '+value+'; value uniform '+value+';'
            elif tag in ['wall','span'] and name=='U':bc='type slip;'
            else:bc='type zeroGradient;'
            rows.append(tag+' { '+bc+' }')
        (case/'0'/name).write_text(header(name,'volVectorField' if name=='U' else 'volScalarField')+
            'dimensions ['+dims[name]+'];\ninternalField uniform '+value+';\nboundaryField {\n'+'\n'.join(rows)+'\n}\n')
    (case/'system/controlDict').write_text(header('controlDict')+f'''solver shockFluid;
libs ("{library.as_posix()}");
startFrom startTime; startTime 0; stopAt endTime; endTime 0.004; deltaT 1e-7;
writeControl adjustableRunTime; writeInterval 0.0005; purgeWrite 0; writeFormat ascii;
writePrecision 17; timeFormat general; timePrecision 12; writeCompression off;
runTimeModifiable false; adjustTimeStep yes; maxCo 0.2; maxDeltaT 2e-6;
functions {{ energy {{ type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e phi); writeControl writeTime; }} }}
''')
    return library,dict(pPa=1000.,TK=180.,Mach=3.,velocityMps=[velocity,0.,0.])


def main():
    """模块二：来源/回归/实际库绑定后执行唯一控制，结果资格分开记录。"""
    from runOpenFoamProcess import run
    start=time.monotonic();out=BASE/'control_run01';out.mkdir(exist_ok=False)
    inputs=[BASE/'native01/report.json',BASE/'callback_check02/report.json',BASE/'reference02/report.json']
    for p in inputs:assert json.loads(p.read_text(encoding='utf-8'))['passed']
    case=BASE/'control01';case.mkdir(exist_ok=False);library,condition=prepare(case)
    names=['runOpenFoamNistControl.py','nistIntegratedAir.py','openFoamFields.py','openFoamSynchronousControl.py',
        'runOpenFoamProcess.py','openFoamConservation.py']
    for name in names:shutil.copyfile(ROOT/'scripts'/name,out/name)
    locks=[*[ROOT/'scripts'/name for name in names],*inputs,library,
        ROOT/'data/thermo/nist_air_shomate_v1.json',ROOT/'docs/validation/v202/CONTROL_PLAN.md',
        *[p for p in case.rglob('*') if p.is_file()]]
    save(out/'admission.json',dict(model=MODEL,condition=condition,endTimeSeconds=.004,saveIntervalSeconds=.0005,
        maximumSteps=18000,maximumSeconds=180,memoryMiB=2048,outputMiB=96,automaticContinuation=0,
        inputSHA256={str(p):sha(p) for p in locks},physicalQualification=False))
    guard=PhysicalTimeGuard(case,out/'saved_guard.json',maximum_steps=18000,auditor=lambda f:audit_state(f,evaluate))
    process=run(['foamRun'],case,out/'native',180,locks,96,2048,guard)
    log=(out/'native/stdout.txt').read_text();dispatch=all(x in log for x in ['hIcoTabulated','sensibleInternalEnergy','perfectGas','hePsiThermo'])
    actual=sorted(v['timeSeconds'] for v in guard.checked.values());expected=[i*.0005 for i in range(1,9)]
    matches=len(actual)==8 and np.allclose(actual,expected,rtol=0,atol=1e-13)
    report=dict(passed=bool(process['passed'] and dispatch and matches and 'End\n' in log),
        model=MODEL,condition=condition,actualSolverDispatchVerified=dispatch,saved=guard.checked,process=process,
        elapsedSeconds=time.monotonic()-start,physicalQualification=False,newQualifiedHeatAreaM2=0)
    save(out/'report.json',report)
    print(json.dumps({k:report[k] for k in ['passed','actualSolverDispatchVerified','elapsedSeconds']}))
    return report['passed']


if __name__=='__main__':raise SystemExit(0 if main() else 1)
