"""几何长度修正后的唯一真实Euler候选；保留旧关闭入口与独立数值身份。"""
from pathlib import Path
import argparse
import json
import re
import shutil
import time
import numpy as np
from hybridRocketEuler import GEOMETRIC_PROFILE as PROFILE,BUDGET,prepare,check_inputs,SavedFieldGuard
from prepareHybridRocketEulerNative import local,save
from readEilmerNative import sha,need

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v196_validation'
PACKAGE=ROOT/'tools/runtime/eilmer-geometric-length-v1'
EVIDENCE=[ROOT/'docs/validation/v196/ROCKET_PLAN.md',BASE/'small_validation02/report.json',
          BASE/'observation03/report.json',BASE/'build_result.json']

def check_admission(a):
    """模块一：固定候选、预算、真实小例和准确包，拒绝换包/改域/放宽温区。"""
    need(a['schema']=='RAH.GeometricRocketEulerAdmission.v1' and a['profile']==PROFILE,'修正候选身份错误')
    need(a['budget']==BUDGET and a['temperatureRangeK']==[200.,6000.],'修正候选预算/温域改变')
    for key,expected in [('case',BASE/'rocket01'),('preparation',BASE/'preparation01/report.json'),('record',BASE/'run01')]:
        need(local(a[key]).resolve()==expected.resolve(),'候选目录不符: '+key)
    need({local(p).resolve() for p in a['evidenceSHA256']}=={p.resolve() for p in EVIDENCE},'准入证据集合改变')
    for p,d in a['evidenceSHA256'].items():need(sha(local(p))==d,'准入证据改变')
    need(json.loads(EVIDENCE[1].read_text())['passed'] and json.loads(EVIDENCE[2].read_text())['passed']
         and json.loads(EVIDENCE[3].read_text())['passed'],'小例/原生观测/构建未通过')
    rp=PACKAGE/'runtime.json';need(sha(rp)==a['runtimeSHA256'],'运行包来源改变')
    pkg=json.loads(rp.read_text())
    need(pkg['schema']=='RAH.EilmerGeometryRuntime.v1' and pkg['profile']=='rah-geometric-length-v1'
         and pkg['sameNumericsAsBaseline'] is False,'运行包数值身份错误')
    for p,d in pkg['filesSHA256'].items():need(sha(local(p))==d,'运行包文件改变')
    need(sha(local(pkg['patch']))==pkg['patchSHA256'],'几何修正来源改变')
    need(sha(local(pkg['buildReceipt']))==pkg['buildReceiptSHA256'],'构建回执改变')
    return pkg

def prepare_candidate():
    """模块二：复用原完整几何和初场，不重新生成网格或原生prep。"""
    prepare(ROOT/'outputs/v194_validation/rocket01',ROOT/'outputs/v194_validation/preparation01/report.json',
            ROOT/'outputs/v194_validation/native_audit01/report.json',BASE/'rocket01',BASE/'preparation01',profile=PROFILE)
    admission=dict(schema='RAH.GeometricRocketEulerAdmission.v1',profile=PROFILE,budget=BUDGET,
        case=str(BASE/'rocket01'),preparation=str(BASE/'preparation01/report.json'),record=str(BASE/'run01'),
        runtimeSHA256=sha(PACKAGE/'runtime.json'),temperatureRangeK=[200.,6000.],
        evidenceSHA256={str(p):sha(p) for p in EVIDENCE},physicalQualification=False)
    check_admission(admission);save(BASE/'admission.json',admission)
    return dict(passed=True)

def execute():
    """模块三：唯一排他原生任务，同时观察完整场和实时Step日志。"""
    raise ValueError('full-wall-euler-geometric-400-v1已因原生更新失败及异常场负温/负压关闭；禁止重跑或接续，须另有成熟稳健方案与预算')
    # 原准确执行源已保存于run01/executed_entry.py；以下保留执行过程供追溯。
    from runEilmerProcess import run
    started=time.monotonic();ap=BASE/'admission.json';a=json.loads(ap.read_text());pkg=check_admission(a)
    case=local(a['case']);prep=local(a['preparation']);out=local(a['record']);r=check_inputs(case,prep)
    need(r['profile']==PROFILE and not out.exists(),'准备身份不符或任务已运行')
    marker=case/'geometric-euler-attempt.json'
    with marker.open('x') as f:json.dump(dict(admissionSHA256=sha(ap),automaticContinuation=0),f)
    out.mkdir();shutil.copyfile(__file__,out/'executed_entry.py')
    shutil.copyfile(ROOT/'scripts/hybridRocketEuler.py',out/'executed_shared.py')
    guard=SavedFieldGuard(case,64982,out/'saved_field_guard.json',maximum_dt=.001)
    def observe():
        bad=guard()
        if bad:return bad
        logfile=out/'native/stdout.txt'
        if logfile.exists():
            for step,dt in re.findall(r'Step=\s*(\d+)\s+t=\s*\S+\s+dt=\s*(\S+)',logfile.read_text()):
                x=float(dt)
                if not np.isfinite(x) or not 0<x<=.001:
                    bad=dict(kind='invalid-live-step-dt',step=int(step),dt=str(dt));save(out/'time_guard.json',bad);return bad
        return None
    locks=[*pkg['filesSHA256'],PACKAGE/'runtime.json',ap,prep,marker,Path(__file__),ROOT/'scripts/hybridRocketEuler.py']
    locks += [local(p) for key in ['sourceSHA256','nativeInputSHA256'] for p in r[key]]+EVIDENCE
    process=run([pkg['binary'],'--max-cpus=2','--max-wall-clock=600'],case,out/'native',
                630-(time.monotonic()-started),128,locks,4096,4096,checkpoint_guard=observe)
    result=dict(schema='RAH.GeometricRocketEulerRun.v1',profile=PROFILE,processPassed=process['passed'],
        case=str(case),preparation=str(prep),admission=str(ap),admissionSHA256=sha(ap),
        outputSHA256={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()},
        physicalQualification=False,newSolvedHeatAreaM2=0,elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',result)
    print(json.dumps(dict(processPassed=process['passed'],nativeExitCode=process['nativeExitCode'],
                         stopReason=process['stopReason'],elapsedSeconds=result['elapsedSeconds'])))
    return result

def validate_terminal(case,preparation,rp):
    """模块四：正常且来源完整才允许研究读取，热/外缘资格仍独立未通过。"""
    case=Path(case).resolve();rp=Path(rp);r=json.loads(rp.read_text());p=json.loads((rp.parent/'native/report.json').read_text())
    need(r['schema']=='RAH.GeometricRocketEulerRun.v1' and r['profile']==PROFILE and r['processPassed']
         and p['passed'] and p['nativeExitCode']==0 and p['inputsStable'] and p['stopReason'] is None
         and not p['liveGroupAfterCleanup'],'异常修正Euler结果不准生产导出')
    need(local(r['case']).resolve()==case and local(p['cwd']).resolve()==case
         and local(r['preparation']).resolve()==Path(preparation).resolve(),'运行来源目录错配')
    ap=local(r['admission']);need(sha(ap)==r['admissionSHA256'],'准入记录改变')
    a=json.loads(ap.read_text());pkg=check_admission(a)
    need(local(a['record']).resolve()==rp.parent.resolve() and local(a['case']).resolve()==case,'准入任务身份错配')
    need(p['argv']==[pkg['binary'],'--max-cpus=2','--max-wall-clock=600']
         and p['inputSHA256'].get(pkg['binary'])==pkg['filesSHA256'][pkg['binary']],'实际运行包错配')
    for path,d in r['outputSHA256'].items():need(sha(local(path))==d,'原生结果改变')
    need(check_inputs(case,preparation)['profile']==PROFILE,'准备数值身份错配')
    need('FINAL-STEP: 400' in (rp.parent/'native/stdout.txt').read_text(),'未完成具名400更新')
    return r,p

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','execute']);args=parser.parse_args()
    result=prepare_candidate() if args.stage=='prepare' else execute()
    if args.stage=='prepare':print(json.dumps(result))
    raise SystemExit(0 if result.get('passed',result.get('processPassed',False)) else 1)
