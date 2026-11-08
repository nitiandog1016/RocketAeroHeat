"""标准干空气的唯一原全壁发展：固定预算、原面和联合物性域。"""
from pathlib import Path
import json,re,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save,header
from standardDryAir import StandardDryAir,MODEL,DOMAIN,SOURCE as TABLE
from openFoamAirThermo import evaluate as old_evaluate
from openFoamFields import audit_state,read_state
from openFoamConservation import geometry,boundary_values,native_mass_flux
from openFoamSynchronousControl import PhysicalTimeGuard
from runOpenFoamRocketDevelopment import validate_preparation,local_path

BASE=ROOT/'outputs/v205_validation'
SOURCE=ROOT/'outputs/v200_validation/rocket03'
LIBRARY=ROOT/'outputs/v204_validation/native01/library/librahStandardDryAir.so'
PROFILE='of13-standard-dry-air-original-rocket-sync-v1'
ENERGY='rahStandardDryAirInternalEnergy'
TIMES=np.arange(1,7)*.000075
GAS=StandardDryAir()

def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8'))

def prerequisites():
    """模块一：只复用已过证据并检查准确身份，不重跑旧验证。"""
    evidence=[ROOT/'outputs/v204_validation/native_audit02/report.json',
        ROOT/'outputs/v204_validation/domain_tests04/report.json',
        ROOT/'outputs/v201_validation/budget01/report.json',
        ROOT/'outputs/v200_validation/rocket_relocation_audit01/report.json',
        ROOT/'outputs/v200_validation/rocket_check03/report.json']
    native,python,budget,geo,check=[read_json(p) for p in evidence]
    if not all(r.get('passed') for r in [native,python,budget,geo,check]):raise ValueError('前置证据未过')
    if native['model']!=MODEL or native['domain']!=DOMAIN or not native['nativeTemplateVerified']:
        raise ValueError('新模型原生模板身份错误')
    if sha(LIBRARY)!=native['registrationLibrarySHA256']:raise ValueError('实际库SHA改变')
    mixture=ROOT/'outputs/v204_validation/native01/mixture.dict'
    if sha(mixture)!=native['tableSHA256'] or mixture.read_text()!=GAS.dictionary():raise ValueError('原生Cp表改变')
    for p,h in python['sourceSHA256'].items():
        if sha(local_path(p))!=h:raise ValueError('最终Python源/表改变 '+p)
    validate_preparation(budget,geo)
    return evidence,budget

class Guard(PhysicalTimeGuard):
    """模块二：实际Co、完整内部场和全部边界均使用新联合域。"""
    def __init__(self,case,record,expected_co):
        super().__init__(case,record,count=64982,maximum_steps=40000,auditor=lambda f:audit_state(f,GAS.evaluate))
        self.expected_co=expected_co;self.initial_co=None;self.g=geometry(case);self.boundary_checked=set()
    def __call__(self,log,final):
        text=Path(log).read_text(errors='replace')
        m=re.search(r'Courant Number mean: [\deE.+-]+ max: ([\deE.+-]+)',text)
        if m:
            self.initial_co=float(m.group(1))
            if abs(self.initial_co/self.expected_co-1)>5e-6:return dict(reason='initial-Co-mismatch',observed=self.initial_co,expected=self.expected_co)
        elif final:return dict(reason='missing-initial-Co')
        problem=super().__call__(log,final)
        if problem:return problem
        for name,row in self.checked.items():
            if float(name)>.00045+1e-12:return dict(reason='physical-end-exceeded')
            if name in self.boundary_checked:continue
            try:
                folder=self.case/name;state=read_state(folder,self.count)
                boundary={n:boundary_values(folder,n,state[n],self.g) for n in state}
                row['boundaryAudit']=audit_state(boundary,GAS.evaluate)
                native_mass_flux(folder,self.g)
                self.boundary_checked.add(name)
            except (ValueError,OSError) as error:return dict(reason='saved-boundary-rejected',time=name,error=str(error))
        save(self.record,dict(frames=self.checked,physicalQualification=False))
        return None

def execute():
    """模块三：独立新模型任务、准确源快照与事前准入；目录防重复启动。"""
    raise ValueError('v205唯一真实候选因输入内能低于100K表界关闭；不得自动降低下界、补场或扫描重跑。')
    from runOpenFoamProcess import run
    start=time.monotonic();evidence,budget=prerequisites()
    BASE.mkdir(exist_ok=True)
    out=BASE/'rocket_run01';out.mkdir(exist_ok=False)
    case=BASE/'rocket01';shutil.copytree(SOURCE,case)
    condition=read_json(ROOT/'outputs/v194_validation/rocket_design03/manifest.json')['condition']
    p,T=condition['pressurePa'],condition['temperatureK']
    a=float(GAS.evaluate(p,T)['soundSpeedMps']);old_a=float(old_evaluate(p,T)['soundSpeedMps'])
    u=np.array(condition['velocityDirection'])*condition['Mach']*a
    expected_co=budget['initialMaxCoAtDt1e7']*a/old_a
    value='('+' '.join(format(x,'.17g') for x in u)+')'
    path=case/'0/U';text,count=re.subn(r'uniform\s+\([^;]+\);','uniform '+value+';',path.read_text())
    if count!=3:raise ValueError('原速度初场/边界格式改变')
    path.write_text(text)
    (case/'constant/physicalProperties').write_text(header('physicalProperties')+
        'thermoType { type hePsiThermo; mixture pureMixture; transport const; thermo hIcoTabulated; equationOfState perfectGas; specie specie; energy '+ENERGY+'; }\nmixture {\n'+GAS.dictionary()+'}\n')
    save(case/'thermo_contract.json',GAS.contract())
    (case/'system/controlDict').write_text(header('controlDict')+f'''solver shockFluid;
libs ("{LIBRARY.as_posix()}");
startFrom startTime; startTime 0; stopAt endTime; endTime 0.00045; deltaT 1e-7;
writeControl adjustableRunTime; writeInterval 0.000075; purgeWrite 0; writeFormat ascii;
writePrecision 17; timeFormat general; timePrecision 12; writeCompression off;
runTimeModifiable false; adjustTimeStep yes; maxCo 0.2; maxDeltaT 2e-6;
functions {{ energy {{ type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e phi); writeControl writeTime; }} }}
''')
    for f in (SOURCE/'constant/polyMesh').iterdir():
        if sha(f)!=sha(case/'constant/polyMesh'/f.name):raise ValueError('原网格改变')
    for n in ['exact_mesh.npz','wall_mapping.json','patches.json','system/fvSchemes','system/fvSolution','0/p','0/T']:
        if sha(SOURCE/n)!=sha(case/n):raise ValueError('原输入改变 '+n)
    names=['runStandardDryAirRocket.py','standardDryAir.py','runOpenFoamRocketDevelopment.py','prepareOpenFoamControl.py',
        'openFoamAirThermo.py','openFoamFields.py','openFoamConservation.py','openFoamSynchronousControl.py',
        'runOpenFoamProcess.py','runEilmerProcess.py']
    for n in names:shutil.copyfile(ROOT/'scripts'/n,out/n)
    locks=[*evidence,LIBRARY,TABLE,ROOT/'docs/validation/v205/RUN_PLAN.md',
        ROOT/'outputs/v194_validation/rocket_design03/manifest.json',
        *[ROOT/'scripts'/n for n in names],*[f for f in case.rglob('*') if f.is_file()]]
    admission=dict(profile=PROFILE,model=MODEL,domain=DOMAIN,energyType=ENERGY,contract=GAS.contract(),
        case=case.relative_to(ROOT).as_posix(),condition=condition,velocityMps=u.tolist(),initialMaxCoExpected=expected_co,
        initialDeltaTExpected=.2/expected_co*1e-7,endTimeSeconds=.00045,saveTimesSeconds=TIMES.tolist(),
        maximumSteps=40000,maximumSeconds=900,parentSeconds=960,memoryMiB=2048,outputMiB=160,automaticContinuation=0,
        registrationLibrarySHA256=sha(LIBRARY),inputSHA256={str(f):sha(f) for f in locks},physicalQualification=False)
    save(out/'admission.json',admission)
    guard=Guard(case,out/'saved_guard.json',expected_co)
    process=run(['foamRun'],case,out/'native',900,[*locks,out/'admission.json'],160,2048,guard)
    times=sorted(float(k) for k in guard.checked);log=(out/'native/stdout.txt').read_text(encoding='utf-8')
    complete=len(times)==6 and np.allclose(times,TIMES,rtol=0,atol=1e-13)
    dispatch=bool(re.search(r'energy\s+'+ENERGY+r'\s*;',log) and re.search(r'thermo\s+hIcoTabulated\s*;',log))
    result=dict(passed=bool(process['passed'] and complete and dispatch and 'End\n' in log),profile=PROFILE,model=MODEL,domain=DOMAIN,
        process=process,saved=guard.checked,actualSolverDispatchVerified=dispatch,initialMaxCo=guard.initial_co,
        admissionSHA256=sha(out/'admission.json'),expectedSavedTimesComplete=bool(complete),
        numericallyConverged=False,attachedEdgeQualified=False,physicalQualification=False,newQualifiedHeatAreaM2=0,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    print(json.dumps(dict(passed=result['passed'],seconds=result['elapsedSeconds'],stop=process['stopReason'],saved=times)))
    return result['passed']

if __name__=='__main__':raise SystemExit(0 if execute() else 1)
