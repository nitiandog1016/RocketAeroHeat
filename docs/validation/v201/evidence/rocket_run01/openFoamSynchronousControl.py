"""原生同步Euler/一致面状态控制；所有推进仍由固定OpenFOAM完成。"""
from pathlib import Path
import argparse,json,shutil,time,re
from prepareOpenFoamControl import ROOT,sha,save,header
from openFoamFields import read_state,audit_state

BASE=ROOT/'outputs/v199_validation';OLD=ROOT/'outputs/v198_validation/control03'
PROFILE='of13-sync-euler-upwind-nasa7-v1'


class PhysicalTimeGuard:
    """模块一：完整保存场与原生警告同时监督；不把物理时间误作伪步编号。"""
    def __init__(self,case,record,count=14892,maximum_steps=16000):
        self.case=Path(case);self.record=Path(record);self.count=count;self.maximum_steps=maximum_steps;self.checked={}
    def __call__(self,log,final):
        text=Path(log).read_text(errors='replace')
        for token in ['out of temperature range','Negative initial temperature','Maximum number of iterations exceeded','FOAM FATAL']:
            if token in text:return dict(reason='native-warning-or-fatal',token=token)
        if len(re.findall(r'^Time = ',text,re.M))>self.maximum_steps:return dict(reason='step-budget')
        candidates=[]
        for d in self.case.iterdir():
            if not d.is_dir():continue
            try:t=float(d.name)
            except ValueError:continue
            if t>0:candidates.append((t,d))
        for t,d in sorted(candidates):
            if d.name in self.checked:continue
            paths=[d/n for n in ['p','T','rho','U','e','phi']]
            complete=all(p.is_file() and p.read_bytes().rstrip().endswith(b'// ************************************************************************* //') for p in paths)
            if not complete and not final:continue
            try:
                row=audit_state(read_state(d,self.count));row['timeSeconds']=t
                row['fieldsSHA256']={p.name:sha(p) for p in paths};self.checked[d.name]=row
            except (ValueError,OSError) as error:return dict(reason='saved-state-rejected',timeSeconds=t,error=str(error))
            save(self.record,dict(frames=self.checked,physicalQualification=False))
        return None


def prepare(case,small):
    """模块二：仅改原生算子字典，原准确网格/气体/边界逐字复用。"""
    source=BASE/'source_review01/report.json';assert json.loads(source.read_text())['passed']
    prep=json.loads((ROOT/'outputs/v198_validation/preparation_repaired.json').read_text())
    for p,h in prep['inputSHA256'].items():assert sha(Path(p))==h
    shutil.copytree(OLD,case)
    p=case/'system/fvSchemes';s=p.read_text()
    for a,b in [('default localEuler;','default Euler;'),('reconstruct(rho) vanLeer;','reconstruct(rho) upwind;'),
        ('reconstruct(U) vanLeerV;','reconstruct(U) upwind;'),('reconstruct(T) vanLeer;','reconstruct(T) upwind;')]:
        assert s.count(a)==1;s=s.replace(a,b)
    p.write_text(s)
    p=case/'system/fvSolution';s=p.read_text()
    s=s.replace(' U {',' "U.*" {').replace(' e {',' "e.*" {')
    s=s.replace('PIMPLE { maxCo 0.8; rDeltaTSmoothingCoeff 0.02; }','PIMPLE {}');p.write_text(s)
    dt='1e-7' if small else '1e-7';end='4e-7' if small else '0.003'
    control='timeStep' if small else 'adjustableRunTime';interval='1' if small else '0.0005'
    (case/'system/controlDict').write_text(header('controlDict')+f'''solver shockFluid;
startFrom startTime; startTime 0; stopAt endTime; endTime {end}; deltaT {dt};
writeControl {control}; writeInterval {interval}; purgeWrite 0; writeFormat ascii;
writePrecision 17; timeFormat general; timePrecision 12; writeCompression off;
runTimeModifiable false; adjustTimeStep {"no" if small else "yes"}; maxCo 0.2; maxDeltaT 2e-6;
functions {{ energy {{ type writeObjects; libs ("libutilityFunctionObjects.so"); objects (e phi); writeControl writeTime; }} }}
''')
    for p in (OLD/'constant/polyMesh').iterdir():assert sha(p)==sha(case/'constant/polyMesh'/p.name)
    assert sha(OLD/'constant/physicalProperties')==sha(case/'constant/physicalProperties')
    return [source,ROOT/'docs/validation/v199/PLAN.md']


def execute(small=True,recover=False):
    """模块三：4步接口通过后才能另登记发展预算；此阶段默认不允许长算。"""
    raise ValueError('v199原计划被监测IO错误中断；本轮四步/发展重跑入口关闭。新任务须读取NEXT_DECISION并另登记来源与剩余预算。')
    # 原准确执行源在small_run01/development_run01/development_run02归档。
    from runOpenFoamProcess import run
    start=time.monotonic();used=0
    if recover:
        assert not small
        prior=json.loads((BASE/'development_run01/native/report.json').read_text())
        assert prior['stopReason']=='supervisor-error' and prior['supervisorError']=="OSError(61, 'No data available')"
        assert not prior['liveGroupAfterCleanup'] and all(prior['inputStable'].values())
        used=prior['elapsedSeconds']
    case=BASE/('small01' if small else ('development02' if recover else 'development01'))
    out=BASE/('small_run01' if small else ('development_run02' if recover else 'development_run01'));out.mkdir(parents=True,exist_ok=False)
    additional=[]
    if not small:
        additional=[BASE/'small_run01/report.json',BASE/'small_audit01/report.json',ROOT/'docs/validation/v199/DEVELOPMENT_PLAN.md']
        for p in additional[:2]:assert json.loads(p.read_text())['passed']
    evidence=prepare(case,small)+additional
    names=['openFoamSynchronousControl.py','openFoamFields.py','openFoamAirThermo.py','runOpenFoamProcess.py']
    for name in names:shutil.copyfile(ROOT/'scripts'/name,out/name)
    locks=[*[ROOT/'scripts'/n for n in names],*evidence,*[p for p in case.rglob('*') if p.is_file()]]
    save(out/'admission.json',dict(profile=PROFILE,steps=4 if small else 16000,initialDt=1e-7,
        physicalEndTimeSeconds=4e-7 if small else .003,budgetSeconds=60 if small else 240-used,memoryMiB=2048,outputMiB=96,
        inputSHA256={str(p):sha(p) for p in locks},automaticContinuation=0,physicalQualification=False))
    guard=PhysicalTimeGuard(case,out/'saved_guard.json',maximum_steps=4 if small else 16000)
    p=run(['foamRun'],case,out/'native',60 if small else 240-used,locks,96,2048,guard)
    expected=[i*1e-7 for i in range(1,5)] if small else [i*.0005 for i in range(1,7)]
    actual=sorted(v['timeSeconds'] for v in guard.checked.values())
    matches=len(actual)==len(expected) and all(abs(a-b)<1e-13 for a,b in zip(actual,expected))
    log=(out/'native/stdout.txt').read_text();passed=p['passed'] and matches and 'End\n' in log
    r=dict(passed=passed,profile=PROFILE,process=p,saved=guard.checked,elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    save(out/'report.json',r);print(json.dumps(dict(passed=passed,elapsedSeconds=r['elapsedSeconds'],stop=p['stopReason'])))
    return passed

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--development',action='store_true')
    parser.add_argument('--recover-observation',action='store_true');args=parser.parse_args()
    raise SystemExit(0 if execute(not args.development,args.recover_observation) else 1)
