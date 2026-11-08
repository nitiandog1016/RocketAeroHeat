"""真实失败拒绝回归：错误身份、假完成和旧模型不得向工程热链传递。"""
from pathlib import Path
import copy,json,subprocess,sys,time,shutil
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepareOpenFoamControl import save,sha
from auditPhaseAwareAirRocket import BASE,read_json,verify_identity,parse_energy_failure,require_developed_frame
from readOpenFoamNistRocketState import verify_identity as old_identity

def main():
    start=time.monotonic();out=BASE/'failure_tests01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py')
    run=read_json(BASE/'rocket_run01/report.json');admission=read_json(BASE/'rocket_run01/admission.json');case=BASE/'rocket01'
    verify_identity(case,run,admission)
    text=(BASE/'rocket_run01/native/stdout.txt').read_text(encoding='utf-8')
    observed=parse_energy_failure(text);rejected=[]
    def rejects(label,action):
        try:action()
        except ValueError as error:rejected.append(dict(case=label,message=str(error)));return
        raise AssertionError('未拒绝 '+label)
    # 模块一：真实任务副本变更身份/来源，不触碰生产数据。
    for field,value in [('profile','old-nist'),('model','NASA7'),('domain','widened'),
                        ('actualSolverDispatchVerified',False),('admissionSHA256','wrong')]:
        changed=copy.deepcopy(run);changed[field]=value
        rejects('run-'+field,lambda r=changed:verify_identity(case,r,admission))
    for field,value in [('model','old-nist'),('domain','other'),('energyType','sensibleInternalEnergy'),
                        ('case','outputs/other'),('velocityMps',[-1418.2630866,0,516.2055479]),
                        ('registrationLibrarySHA256','wrong')]:
        changed=copy.deepcopy(admission);changed[field]=value
        rejects('admission-'+field,lambda a=changed:verify_identity(case,run,a))
    changed=copy.deepcopy(run);changed['process']['liveGroupAfterCleanup']={'1':{}}
    rejects('live-group',lambda:verify_identity(case,changed,admission))
    changed=copy.deepcopy(run);changed['process']['inputStable']={'x':False}
    rejects('unstable-input',lambda:verify_identity(case,changed,admission))
    changed=copy.deepcopy(run);changed['process']['logExportPassed']=False
    rejects('unverified-log',lambda:verify_identity(case,changed,admission))
    rejects('old-reader-new-gas',lambda:old_identity(case,run,admission))
    rejects('failed-task-as-edge',lambda:require_developed_frame(run,'0.000075'))
    fake=copy.deepcopy(run);fake['passed']=True
    rejects('fake-normal-with-no-frame',lambda:require_developed_frame(fake,'0.000075'))
    fake['saved']={'0':{}}
    rejects('initial-as-developed',lambda:require_developed_frame(fake,'0'))
    fake['saved']={'0.000075':{}}
    rejects('fake-unqualified-edge',lambda:require_developed_frame(fake,'0.000075'))
    rejects('missing-fatal',lambda:parse_energy_failure(''))
    rejects('duplicated-fatal',lambda:parse_energy_failure(text+text))
    rejects('wrong-energy-order',lambda:parse_energy_failure('RAH_PHASE_ENERGY_OUTSIDE_TABLE e=2 elo=1 ehi=3 pOld=70 seed=60'))
    # 模块二：执行关闭CLI仅检查启动前拒绝，不能产生第二个原生任务。
    before=sha(BASE/'rocket_run01/report.json')
    closed=subprocess.run([sys.executable,str(ROOT/'scripts/runPhaseAwareAirRocket.py')],capture_output=True,timeout=10)
    if closed.returncode==0 or 'v206' not in closed.stderr.decode('utf-8',errors='replace'):raise AssertionError('候选未在启动前关闭')
    (out/'closed_cli_stderr.txt').write_bytes(closed.stderr)
    if sha(BASE/'rocket_run01/report.json')!=before:raise AssertionError('关闭CLI改变旧结果')
    protected={k:sha(ROOT/k.replace('\\','/'))==h for k,h in read_json(ROOT/'outputs/v176_validation/protected_initial.json').items()}
    if not all(protected.values()):raise AssertionError('保护文件改变')
    result=dict(passed=True,identityPositiveChecks=1,rejections=rejected,closedCLIVerified=True,
        observedFailure=observed,positiveDevelopedReaderTested=False,protected=protected,nativeRuns=0,
        physicalQualification=False,newQualifiedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(dict(passed=True,rejections=len(rejected),seconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
