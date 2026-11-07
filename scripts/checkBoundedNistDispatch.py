"""已有完整控制场的一次原生新映射更新；验证分派及守恒，非发展重跑。"""
from pathlib import Path
import json,re,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save,header
from openFoamSynchronousControl import PhysicalTimeGuard
from openFoamFields import audit_state,read_state
from openFoamConservation import geometry,flux_and_residual,conservative,native_mass_flux
from nistIntegratedAir import evaluate,R


def main():
    raise ValueError('v203唯一原生新映射一步分派已通过，入口关闭；不重复控制。')
    from runOpenFoamProcess import run
    start=time.monotonic();base=ROOT/'outputs/v203_validation';out=base/'dispatch01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');source=ROOT/'outputs/v202_validation/control01';case=base/'dispatch_case01';case.mkdir(exist_ok=False)
    library=base/'bounded01/library/librahBoundedNist.so';br=json.loads((base/'bounded01/report.json').read_text());assert br['passed'] and sha(library)==br['registrationLibrarySHA256']
    old=json.loads((ROOT/'outputs/v202_validation/control_run01/report.json').read_text());assert old['passed']
    for folder in ['constant','system']:shutil.copytree(source/folder,case/folder)
    for name in ['exact_mesh.npz','patches.json','wall_identity.json','thermo_contract.json']:shutil.copyfile(source/name,case/name)
    initial=case/'0';initial.mkdir()
    for name,h in old['saved']['0.004']['fieldsSHA256'].items():
        p=source/'0.004'/name;assert sha(p)==h
        text=p.read_text();text,n=re.subn(r'location\s+"0.004";', 'location "0";',text);assert n==1
        (initial/name).write_text(text)
    props=case/'constant/physicalProperties';text=props.read_text();assert text.count('energy sensibleInternalEnergy;')==1
    props.write_text(text.replace('energy sensibleInternalEnergy;','energy rahBoundedNistInternalEnergy;'))
    (case/'system/controlDict').write_text(header('controlDict')+f'''solver shockFluid;
libs ("{library.as_posix()}");
startFrom startTime; startTime 0; stopAt endTime; endTime 1e-9; deltaT 1e-9;
writeControl timeStep; writeInterval 1; purgeWrite 0; writeFormat ascii;
writePrecision 17; timeFormat general; timePrecision 12; writeCompression off;
runTimeModifiable false; adjustTimeStep no;
functions {{ energy {{ type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e phi); writeControl writeTime; }} }}
''')
    names=['checkBoundedNistDispatch.py','openFoamFields.py','openFoamSynchronousControl.py','openFoamConservation.py','nistIntegratedAir.py','runOpenFoamProcess.py']
    for name in names:shutil.copyfile(ROOT/'scripts'/name,out/name)
    locks=[*[ROOT/'scripts'/name for name in names],library,base/'bounded01/report.json',ROOT/'docs/validation/v203/DISPATCH_PLAN.md',
        *[source/'0.004'/n for n in old['saved']['0.004']['fieldsSHA256']],*[p for p in case.rglob('*') if p.is_file()]]
    save(out/'admission.json',dict(seconds=30,steps=1,dt=1e-9,sourceTime=.004,continuousRestart=False,automaticContinuation=0,inputSHA256={str(p):sha(p) for p in locks}))
    guard=PhysicalTimeGuard(case,out/'saved_guard.json',maximum_steps=1,auditor=lambda f:audit_state(f,evaluate))
    process=run(['foamRun'],case,out/'native',30,locks,64,2048,guard)
    save(out/'run_report.json',dict(passed=process['passed'],process=process,saved=guard.checked))
    assert process['passed'] and set(guard.checked)=={'1e-09'}
    text=(out/'native/stdout.txt').read_text();assert 'rahBoundedNistInternalEnergy' in text and len(re.findall(r'^Time = ',text,re.M))==1
    g=geometry(case);f,flux,res=flux_and_residual(case/'0',g,evaluate,R);following=read_state(case/'1e-09',len(g['cells']));audit_state(following,evaluate)
    # 同时检查全部输出边界与空间残差，不只看内部单元。
    flux_and_residual(case/'1e-09',g,evaluate,R)
    gas=evaluate(1000.,180.);rho=float(gas['densityKgm3']);speed=3*float(gas['soundSpeedMps'])
    scales=np.array([rho,*([rho*speed]*3),rho*(abs(float(gas['internalEnergyJkg']))+speed**2)])
    delta=conservative(following)-conservative(f);defect=delta+1e-9*res;errors=np.max(abs(defect)/scales,axis=0)
    phi=float(np.max(abs(native_mass_flux(case/'1e-09',g)-flux[:,0]))/(scales[1]*g['area'].max()))
    balance=abs(np.sum(delta*g['volume'][:,None],axis=0)+1e-9*flux[len(g['neighbor']):].sum(axis=0))/(scales*g['volume'].sum())
    assert max(errors)<2e-10 and phi<2e-10 and max(balance)<2e-10
    result=dict(passed=True,actualSolverDispatchVerified=True,internalAndBoundaryEOSPassed=True,conservativeUpdateMaxScaled=errors.tolist(),
        nativePhiMaxScaled=phi,globalBalanceScaled=balance.tolist(),process=process,saved=guard.checked,elapsedSeconds=time.monotonic()-start,
        physicalQualification=False,newQualifiedHeatAreaM2=0)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','conservativeUpdateMaxScaled','nativePhiMaxScaled','elapsedSeconds']}))


if __name__=='__main__':main()
