"""低温物性闭合后的唯一原全壁Euler发展，不重开旧模型失败。"""
from pathlib import Path
import json,re,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save,header
from nistIntegratedAir import AIR,evaluate,contract,MODEL
from openFoamAirThermo import evaluate as old_evaluate
from openFoamFields import audit_state
from openFoamSynchronousControl import PhysicalTimeGuard
from runOpenFoamRocketDevelopment import validate_preparation,local_path

BASE=ROOT/'outputs/v202_validation';SOURCE=ROOT/'outputs/v200_validation/rocket03'
PROFILE='of13-nist-original-rocket-sync-v1';TIMES=np.arange(1,7)*.000075


def prerequisites():
    """模块一：已过模型/控制/原几何与实际库准确绑定，不用品牌或exit0替代。"""
    files=[BASE/'native01/report.json',BASE/'control_run01/report.json',BASE/'control_audit01/report.json',
        ROOT/'outputs/v201_validation/budget01/report.json',ROOT/'outputs/v200_validation/rocket_relocation_audit01/report.json']
    reports=[json.loads(p.read_text(encoding='utf-8')) for p in files]
    if not all(r.get('passed') for r in reports):raise ValueError('新气体/控制/几何前置不完整')
    native,control,audit,budget,geo=reports
    if not control.get('actualSolverDispatchVerified') or not audit.get('numericallyConverged'):raise ValueError('实际分派/真实残差未过')
    if control.get('model')!=MODEL or native.get('contract')!=contract():raise ValueError('模型合同改变')
    validate_preparation(budget,geo)
    lib=BASE/'native01/library/librahNistAir.so'
    if sha(lib)!=native['registrationLibrarySHA256']:raise ValueError('实际注册库改变')
    if (BASE/'native01/mixture.dict').read_text(encoding='utf-8')!=AIR.dictionary():raise ValueError('表数据与原生验证不一致')
    for row in native['query']['inputSHA256'].items():
        name,digest=row
        if name.endswith('/nistIntegratedAir.py') or name.endswith('/nist_air_shomate_v1.json'):
            if sha(local_path(name))!=digest:raise ValueError('物性执行源改变')
    return files,lib,budget


class Guard(PhysicalTimeGuard):
    """模块二：首Co核对齐次缩放；全场/压力域以新EOS明确审核。"""
    def __init__(self,case,record,expected_co):
        super().__init__(case,record,count=64982,maximum_steps=40000,auditor=lambda f:audit_state(f,evaluate))
        self.expected_co=expected_co;self.initial_co=None
    def __call__(self,log,final):
        text=Path(log).read_text(errors='replace')
        m=re.search(r'Courant Number mean: [\deE.+-]+ max: ([\deE.+-]+)',text)
        if m:
            self.initial_co=float(m.group(1))
            if abs(self.initial_co/self.expected_co-1)>5e-6:return dict(reason='initial-Co-mismatch',observed=self.initial_co,expected=self.expected_co)
        elif final:return dict(reason='missing-initial-Co')
        problem=super().__call__(log,final)
        if problem:return problem
        if any(float(t)>.00045+1e-12 for t in self.checked):return dict(reason='physical-end-exceeded')
        return None


def main():
    """模块三：新目录/有限预算/原面保持；不重编引擎、不自动追加。"""
    raise ValueError('v202真实NIST候选因原生100K表下界Fatal关闭；不续算/裁温/降低边界重跑。')
    from runOpenFoamProcess import run
    start=time.monotonic();evidence,library,oldbudget=prerequisites()
    out=BASE/'rocket_run01';out.mkdir(exist_ok=False);case=BASE/'rocket01';shutil.copytree(SOURCE,case)
    condition=json.loads((ROOT/'outputs/v194_validation/rocket_design03/manifest.json').read_text(encoding='utf-8'))['condition']
    p,T=condition['pressurePa'],condition['temperatureK'];newa=float(evaluate(p,T)['soundSpeedMps']);olda=float(old_evaluate(p,T)['soundSpeedMps'])
    u=np.array(condition['velocityDirection'])*condition['Mach']*newa
    expected_co=oldbudget['initialMaxCoAtDt1e7']*newa/olda
    # 同一物理边界及初温/初压；只替换新EOS的速度和热模型，旧源目录不变。
    path=case/'0/U';text=path.read_text();value='('+' '.join(format(x,'.17g') for x in u)+')'
    text,count=re.subn(r'uniform\s+\([^;]+\);','uniform '+value+';',text);assert count==3;path.write_text(text)
    (case/'constant/physicalProperties').write_text(header('physicalProperties')+
        'thermoType { type hePsiThermo; mixture pureMixture; transport const; thermo hIcoTabulated; equationOfState perfectGas; specie specie; energy sensibleInternalEnergy; }\nmixture {\n'+AIR.dictionary()+'}\n')
    save(case/'thermo_contract.json',contract())
    (case/'system/controlDict').write_text(header('controlDict')+f'''solver shockFluid;
libs ("{library.as_posix()}");
startFrom startTime; startTime 0; stopAt endTime; endTime 0.00045; deltaT 1e-7;
writeControl adjustableRunTime; writeInterval 0.000075; purgeWrite 0; writeFormat ascii;
writePrecision 17; timeFormat general; timePrecision 12; writeCompression off;
runTimeModifiable false; adjustTimeStep yes; maxCo 0.2; maxDeltaT 2e-6;
functions {{ energy {{ type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e phi); writeControl writeTime; }} }}
''')
    for pth in (SOURCE/'constant/polyMesh').iterdir():assert sha(pth)==sha(case/'constant/polyMesh'/pth.name)
    for rel in ['exact_mesh.npz','wall_mapping.json','patches.json','system/fvSchemes','system/fvSolution','0/p','0/T']:
        assert sha(SOURCE/rel)==sha(case/rel)
    names=['runOpenFoamNistRocket.py','runOpenFoamRocketDevelopment.py','nistIntegratedAir.py','openFoamAirThermo.py','openFoamFields.py',
        'openFoamSynchronousControl.py','runOpenFoamProcess.py','openFoamConservation.py']
    for n in names:shutil.copyfile(ROOT/'scripts'/n,out/n)
    locks=[*evidence,library,ROOT/'data/thermo/nist_air_shomate_v1.json',ROOT/'docs/validation/v202/ROCKET_PLAN.md',
        ROOT/'outputs/v194_validation/rocket_design03/manifest.json',*[ROOT/'scripts'/n for n in names],*[p for p in case.rglob('*') if p.is_file()]]
    admission=dict(profile=PROFILE,model=MODEL,case=case.relative_to(ROOT).as_posix(),condition=condition,velocityMps=u.tolist(),
        initialMaxCoExpected=expected_co,newOldSoundSpeedRatio=newa/olda,initialDeltaTExpected=.2/expected_co*1e-7,
        endTimeSeconds=.00045,saveTimesSeconds=TIMES.tolist(),maximumSteps=40000,maximumSeconds=900,parentSeconds=960,
        memoryMiB=2048,outputMiB=160,automaticContinuation=0,registrationLibrarySHA256=sha(library),
        inputSHA256={str(p):sha(p) for p in locks},physicalQualification=False)
    save(out/'admission.json',admission)
    guard=Guard(case,out/'saved_guard.json',expected_co)
    process=run(['foamRun'],case,out/'native',900,[*locks,out/'admission.json'],160,2048,guard)
    times=sorted(float(k) for k in guard.checked);log=(out/'native/stdout.txt').read_text(encoding='utf-8')
    complete=len(times)==6 and np.allclose(times,TIMES,rtol=0,atol=1e-13)
    dispatch='hIcoTabulated' in log and 'sensibleInternalEnergy' in log
    report=dict(passed=bool(process['passed'] and complete and dispatch and 'End\n' in log),profile=PROFILE,model=MODEL,
        process=process,saved=guard.checked,actualSolverDispatchVerified=dispatch,initialMaxCo=guard.initial_co,
        admissionSHA256=sha(out/'admission.json'),expectedSavedTimesComplete=bool(complete),
        numericallyConverged=False,attachedEdgeQualified=False,physicalQualification=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(dict(passed=report['passed'],seconds=report['elapsedSeconds'],stop=process['stopReason'],saved=times)))
    return report['passed']


if __name__=='__main__':raise SystemExit(0 if main() else 1)
