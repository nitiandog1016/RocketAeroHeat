"""真实新模型失败拒绝与原面身份检查；不启动原生任务。"""
from pathlib import Path
import copy,json,sys,subprocess,shutil,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from readOpenFoamNistRocketState import read_frame,verify_identity,BASE
from prepareOpenFoamControl import save


def main():
    start=time.monotonic();out=BASE/'rocket_tests01';out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/'executed_test.py')
    run=json.loads((BASE/'rocket_run01/report.json').read_text());admission=json.loads((BASE/'rocket_run01/admission.json').read_text())
    rows=[]
    def reject(name,action):
        try:action()
        except (ValueError,AssertionError):rows.append(name);return
        raise AssertionError('未拒绝 '+name)
    verify_identity(BASE/'rocket01',run,admission)
    reject('failed_production',lambda:read_frame('.000075',.000075))
    reject('no_development_frame',lambda:read_frame('.000075',.000075,diagnostic=True))
    for name,t in [('0',0),('../0',0),('0.000075',.00015)]:reject('time_'+name,lambda n=name,t=t:read_frame(n,t,True))
    for key,value in [('profile','old'),('model','NASA7'),('actualSolverDispatchVerified',False),('admissionSHA256','0'*64)]:
        r=copy.deepcopy(run);r[key]=value;reject('run_'+key,lambda r=r:verify_identity(BASE/'rocket01',r,admission))
    for key,value in [('profile','old'),('model','NASA7'),('case','wrong'),('velocityMps',[-1500,0,0]),('condition',{}),('registrationLibrarySHA256','0'*64)]:
        a=copy.deepcopy(admission);a[key]=value;reject('admission_'+key,lambda a=a:verify_identity(BASE/'rocket01',run,a))
    for key,value in [('liveGroupAfterCleanup',{'1':{}}),('inputStable',{'input':False}),('logExportPassed',False)]:
        r=copy.deepcopy(run);r['process'][key]=value;reject('process_'+key,lambda r=r:verify_identity(BASE/'rocket01',r,admission))
    for script in ['runOpenFoamNistControl.py','runOpenFoamNistRocket.py']:
        p=subprocess.run([sys.executable,str(ROOT/'scripts'/script)],capture_output=True,timeout=10)
        assert p.returncode!=0 and b'ValueError' in p.stderr and b'ModuleNotFoundError' not in p.stderr
        rows.append('closed_'+script)
    save(out/'report.json',dict(passed=True,rejections=rows,rejectionCount=len(rows),positiveIdentityChecks=1,
        positiveDevelopmentReadTested=False,CFDRuns=0,elapsedSeconds=time.monotonic()-start))
    print(json.dumps(dict(passed=True,rejections=len(rows),elapsedSeconds=time.monotonic()-start)))


if __name__=='__main__':main()
