"""唯一准确保存场恢复；保留原失败身份，只完成登记的剩余物理区间。"""
from pathlib import Path
import json,re,shutil,time
from prepareOpenFoamControl import ROOT,sha,save
from openFoamSynchronousControl import PhysicalTimeGuard

BASE=ROOT/'outputs/v200_validation';SOURCE=ROOT/'outputs/v199_validation/development02'
PROFILE='of13-sync-saved002-to003-v1'


def validate_source(review,source=SOURCE):
    """模块一：只接受具名准确时刻/物性/来源，不以任意有限场代替。"""
    if not review.get('passed') or review.get('continuousBitwiseEquivalent') is not False:raise ValueError('重启审核身份错误')
    if review.get('sourceFrame')!='outputs/v199_validation/development02/0.002':raise ValueError('错误来源时刻')
    if review.get('time',{}).get('index')!=7889 or review['time'].get('value')!=.002:raise ValueError('时间身份错误')
    if set(review.get('sourceFieldsSHA256',{}))!={'p','T','rho','U','e','phi'}:raise ValueError('缺完整保存量')
    for name,h in review['sourceFieldsSHA256'].items():
        if sha(source/'0.002'/name)!=h:raise ValueError('源字段SHA改变 '+name)
    if sha(source/'0.002/uniform/time')!=review['sourceTimeSHA256']:raise ValueError('时间SHA改变')
    if sha(source/'exact_mesh.npz')!=review['sourceGeometrySHA256']:raise ValueError('几何SHA改变')
    return True


def main():
    raise ValueError('v200唯一恢复已正常完成；入口关闭，直接进入真实原全壁目标，不重复控制计算。')
    # 准确实际执行源保存在restart_run01/runOpenFoamRestart.py。
    from runOpenFoamProcess import run
    start=time.monotonic();review_path=BASE/'restart_review02/report.json'
    review=json.loads(review_path.read_text());validate_source(review)
    for p,h in review['sourceSHA256'].items():assert sha(p)==h,p
    out=BASE/'restart_run01';out.mkdir(exist_ok=False);case=BASE/'restart01';case.mkdir(exist_ok=False)
    # 模块二：只复制准确源帧和必要输入，保持旧场，不复制历史发展序列。
    for name in ['constant','system','0.002']:shutil.copytree(SOURCE/name,case/name)
    for name in ['exact_mesh.npz','patches.json','thermo_contract.json','wall_identity.json','repair.json']:shutil.copyfile(SOURCE/name,case/name)
    control=case/'system/controlDict';text=control.read_text();assert text.count('startTime 0;')==1
    control.write_text(text.replace('startTime 0;','startTime 0.002;'))
    names=['runOpenFoamRestart.py','runOpenFoamProcess.py','openFoamSynchronousControl.py','openFoamFields.py','openFoamAirThermo.py']
    for name in names:shutil.copyfile(ROOT/'scripts'/name,out/name)
    locks=[review_path,ROOT/'docs/validation/v200/RESTART_PLAN.md',*[ROOT/'scripts'/n for n in names],
        *[p for p in case.rglob('*') if p.is_file()],*[SOURCE/'0.002'/n for n in review['sourceFieldsSHA256']],SOURCE/'0.002/uniform/time']
    save(out/'admission.json',dict(profile=PROFILE,source=review,maximumSteps=7000,maximumSeconds=120,
        parentSeconds=150,memoryMiB=2048,outputMiB=96,automaticContinuation=0,oldDevelopmentSeconds=61.204940375,
        inputSHA256={str(p):sha(p) for p in locks},physicalQualification=False))
    guard=PhysicalTimeGuard(case,out/'saved_guard.json',maximum_steps=7000)
    process=run(['foamRun'],case,out/'native',120,locks,96,2048,guard)
    times=sorted(float(k) for k in guard.checked)
    log=(out/'native/stdout.txt').read_text();passed=process['passed'] and times==[.002,.0025,.003] and 'End\n' in log
    result=dict(passed=passed,profile=PROFILE,process=process,saved=guard.checked,source=review,
        oldRunPlanAcceptance=False,continuousBitwiseEquivalent=False,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(dict(passed=passed,seconds=result['elapsedSeconds'],stop=process['stopReason'],saved=times)))
    return passed


if __name__=='__main__':raise SystemExit(0 if main() else 1)
