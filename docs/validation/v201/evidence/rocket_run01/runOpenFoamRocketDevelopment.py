"""真实全壁唯一同步Euler发展：预算、准确来源与原生保存监督。"""
from pathlib import Path
import json,re,shutil,time
from prepareOpenFoamControl import ROOT,sha,save,header
from openFoamSynchronousControl import PhysicalTimeGuard

BASE=ROOT/'outputs/v201_validation'
SOURCE=ROOT/'outputs/v200_validation/rocket03'
PROFILE='of13-original-rocket-sync-development-v1'
END=.00045
PERIOD=.000225
MAX_STEPS=40000


def local_path(value):
    """模块一：同项目Windows/WSL回执路径相容，不改变SHA所指文件。"""
    s=str(value).replace('\\','/')
    for prefix in ['E:/RocketAeroHeat/','/mnt/e/RocketAeroHeat/']:
        if s.startswith(prefix):return ROOT/s[len(prefix):]
    p=Path(s)
    if not p.is_absolute():p=ROOT/p
    return p


def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def validate_preparation(budget,geometry_review):
    """模块二：只有已过全部原壁几何与准确初态成本的具名来源才能准备。"""
    if not budget.get('passed') or not geometry_review.get('passed'):raise ValueError('前置审核未通过')
    if budget.get('cellCount')!=64982 or budget.get('originalWallFaces')!=19752:raise ValueError('网格范围改变')
    if abs(budget.get('originalAreaM2',0)-1.9646292565715553)>1e-13:raise ValueError('原面积改变')
    if budget.get('targetMaxCo')!=.2 or not 1.8e-8<budget.get('initialTargetDeltaTSeconds',0)<1.9e-8:raise ValueError('时间步合同改变')
    if not geometry_review.get('originalWallExact') or not geometry_review.get('connectivityExact'):raise ValueError('原壁/连接未保持')
    for r,key in [(budget,'sourceSHA256'),(geometry_review,'inputSHA256')]:
        if not r.get(key):raise ValueError('缺少准确来源')
        for p,h in r[key].items():
            if sha(local_path(p))!=h:raise ValueError('来源SHA改变 '+p)
    return True


def control_text():
    """模块三：常规两帧和一个早期对象快照，时间推进仍全部由原生引擎承担。"""
    return header('controlDict')+'''solver shockFluid;
startFrom startTime; startTime 0; stopAt endTime; endTime 0.00045; deltaT 1e-7;
writeControl adjustableRunTime; writeInterval 0.000225; purgeWrite 0; writeFormat ascii;
writePrecision 17; timeFormat general; timePrecision 12; writeCompression off;
runTimeModifiable false; adjustTimeStep yes; maxCo 0.2; maxDeltaT 2e-6;
functions
{
 energy { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e phi); writeControl writeTime; }
 early { type writeObjects; libs ("libutilityFunctionObjects.so"); objects (p T rho U e phi);
         executeAtStart false; startTime 1e-8; endTime 2e-8; writeControl timeStep; writeInterval 1; }
}
'''


class RocketGuard(PhysicalTimeGuard):
    """模块四：继承物性/完整场检查，并核对真实首Co、早期输出与固定终点。"""
    def __init__(self,case,record,budget):
        super().__init__(case,record,count=64982,maximum_steps=MAX_STEPS)
        self.budget=budget;self.initial_co=None
    def __call__(self,log,final):
        text=Path(log).read_text(errors='replace')
        m=re.search(r'Courant Number mean: [\deE.+-]+ max: ([\deE.+-]+)',text)
        if m:
            self.initial_co=float(m.group(1))
            if abs(self.initial_co/self.budget['initialMaxCoAtDt1e7']-1)>5e-6:
                return dict(reason='initial-Courant-mismatch',observed=self.initial_co)
        elif final:return dict(reason='missing-initial-Courant')
        problem=super().__call__(log,final)
        if problem:return problem
        early=[float(t) for t in self.checked if float(t)<1e-6]
        if len(early)>1:return dict(reason='unexpected-extra-early-snapshot',times=early)
        if any(float(t)>END+1e-12 for t in self.checked):return dict(reason='physical-end-exceeded')
        return None


def main():
    """模块五：唯一目录避免重复启动；准确执行源先保存，失败不自动追加。"""
    from runOpenFoamProcess import run
    start=time.monotonic()
    budget=read_json(BASE/'budget01/report.json')
    geometry_review=read_json(ROOT/'outputs/v200_validation/rocket_relocation_audit01/report.json')
    validate_preparation(budget,geometry_review)
    evidence=[BASE/'budget01/report.json',ROOT/'docs/validation/v201/RUN_PLAN.md',
        ROOT/'outputs/v200_validation/rocket_relocation_audit01/report.json',
        ROOT/'outputs/v200_validation/rocket_preparation02/report.json',
        ROOT/'outputs/v200_validation/rocket_check03/report.json',
        ROOT/'outputs/v200_validation/restart_audit02/report.json',BASE/'admission_tests01/report.json']
    for p in evidence:
        if p.suffix=='.json' and not read_json(p).get('passed'):raise ValueError('缺前置通过 '+str(p))
    out=BASE/'rocket_run01';out.mkdir(exist_ok=False)
    case=BASE/'rocket01';shutil.copytree(SOURCE,case)
    (case/'system/controlDict').write_text(control_text(),encoding='utf-8')
    names=['runOpenFoamRocketDevelopment.py','readOpenFoamRocketState.py','runOpenFoamProcess.py',
        'openFoamSynchronousControl.py','openFoamFields.py','openFoamAirThermo.py','openFoamConservation.py']
    for n in names:shutil.copyfile(ROOT/'scripts'/n,out/n)
    # 原生时间控制与对象保存源码/实际动态库同时绑定，避免猜测timeStart旧拼法。
    native=['src/OpenFOAM/db/functionObjects/timeControl/timeControlFunctionObject.C',
        'src/OpenFOAM/db/functionObjects/timeControl/timeControl.C',
        'src/functionObjects/utilities/writeObjects/writeObjects.C']
    paths=[]
    for name in native:
        p=Path('/opt/openfoam13')/name
        if not p.is_file():raise ValueError('缺准确保存源码 '+str(p))
        q=out/'native_source'/name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q);paths.append(p)
    paths.append(Path('/opt/openfoam13/platforms/linux64GccDPInt32Opt/lib/libutilityFunctionObjects.so'))
    locks=[*evidence,*paths,*[ROOT/'scripts'/n for n in names],*[p for p in case.rglob('*') if p.is_file()]]
    admission=dict(profile=PROFILE,source='outputs/v200_validation/rocket03',case=case.relative_to(ROOT).as_posix(),
        endTimeSeconds=END,saveTimesSeconds=[PERIOD,END],earlyWindowSeconds=[1e-8,2e-8],
        maximumSteps=MAX_STEPS,maximumSeconds=900,parentSeconds=960,memoryMiB=2048,outputMiB=96,
        automaticContinuation=0,budget=budget,condition=read_json(evidence[3])['condition'],
        velocityMps=read_json(evidence[3])['velocityMps'],
        inputSHA256={str(p):sha(p) for p in locks},physicalQualification=False)
    save(out/'admission.json',admission)
    guard=RocketGuard(case,out/'saved_guard.json',budget)
    process=run(['foamRun'],case,out/'native',900,[*locks,out/'admission.json'],96,2048,guard)
    times=sorted(float(k) for k in guard.checked);log=(out/'native/stdout.txt').read_text(encoding='utf-8')
    early=[t for t in times if 1e-8<=t<=2e-8]
    complete=len(early)==1 and times==[early[0],PERIOD,END]
    passed=bool(process['passed'] and complete and 'End\n' in log)
    result=dict(passed=passed,profile=PROFILE,process=process,saved=guard.checked,initialMaxCo=guard.initial_co,
        admissionSHA256=sha(out/'admission.json'),expectedSavedTimesComplete=complete,
        numericallyConverged=False,attachedEdgeQualified=False,physicalQualification=False,newQualifiedHeatAreaM2=0,
        elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result)
    print(json.dumps(dict(passed=passed,seconds=result['elapsedSeconds'],stop=process['stopReason'],saved=times)))
    return passed


if __name__=='__main__':raise SystemExit(0 if main() else 1)
