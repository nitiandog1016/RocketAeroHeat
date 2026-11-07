"""原全壁复数NK唯一有限候选：绑定真实小控制、预算和完整原面。"""
from pathlib import Path
import argparse
import json
import shutil
import time
import numpy as np
from readEilmerNative import sha, need, read_field
from prepareHybridRocketEulerNative import local, save, verify_design
from hybridNKTetraControl import NK
from hybridRocketEuler import SavedFieldGuard

ROOT=Path(__file__).resolve().parents[1]; BASE=ROOT/'outputs/v197_validation'
CASE=BASE/'rocket01'; OLD=ROOT/'outputs/v194_validation/rocket01'
PACKAGE=ROOT/'tools/runtime/eilmer-nk-geometric-v1'
PROFILE='full-wall-euler-nk-120-v1'
BUDGET=dict(maxSteps=120,nativeSeconds=600,supervisorSeconds=630,parentSeconds=650,
            cpuCount=2,addressSpaceMiB=4096,groupRSSMiB=4096,outputMiB=192,automaticContinuation=0)
EVIDENCE=[BASE/'audit01/report.json',BASE/'build_result.json',
          ROOT/'docs/validation/v197/ENTRY_DECISION.md',ROOT/'docs/validation/v197/ROCKET_PLAN.md']


def check_admission(a):
    """模块一：具名唯一候选；预算、温域、完整原面和准确小控制不得替换。"""
    need(a['schema']=='RAH.HybridNKRocketAdmission.v1' and a['profile']==PROFILE,'NK候选身份错误')
    need(a['budget']==BUDGET and a['temperatureRangeK']==[200.,6000.],'NK预算或温域改变')
    need(a['cells']==64982 and a['wallFaces']==19752 and a['physicalQualification'] is False,'NK原域或资格改变')
    need(a['savedSteps']==list(range(0,121,20)), 'NK保存计划改变')
    for key,expected in [('case',CASE),('preparation',BASE/'rocket_preparation01/report.json'),('record',BASE/'rocket_run01')]:
        need(local(a[key]).resolve()==expected.resolve(),'NK目录身份改变: '+key)
    need({local(p).resolve() for p in a['evidenceSHA256']}=={p.resolve() for p in EVIDENCE},'NK证据集合改变')
    for p,h in a['evidenceSHA256'].items():need(sha(local(p))==h,'NK证据改变')
    control=json.loads(EVIDENCE[0].read_text()); build=json.loads(EVIDENCE[1].read_text())
    need(control['passed'] and control['screeningPassed'] and build['passed'], 'NK控制/构建不合格')
    runtime=PACKAGE/'runtime.json';need(sha(runtime)==a['runtimeSHA256'],'NK运行包改变')
    pkg=json.loads(runtime.read_text())
    need(pkg['schema']=='RAH.EilmerNKGeometryRuntime.v1' and pkg['profile']=='rah-nk-geometric-v1'
         and pkg['numberType']=='complex' and pkg['officialFullNKStorage'] is True,'NK运行包身份错误')
    for p,h in pkg['filesSHA256'].items():need(sha(local(p))==h,'NK实际二进制或源改变')
    need(sha(local(pkg['buildReceipt']))==pkg['buildReceiptSHA256'],'NK构建回执改变')
    need(sha(ROOT/'patches/eilmer/geometric-length-v1.patch')==pkg['geometryPatchSHA256'],'NK几何补丁改变')
    return pkg


def prepare():
    """模块二：原网格/气体逐字复用，仅按新成熟算法重新原生准备初场。"""
    from runEilmerProcess import run
    start=time.monotonic();m=verify_design(OLD)
    oldprep=ROOT/'outputs/v194_validation/preparation01/report.json'
    oldaudit=ROOT/'outputs/v194_validation/native_audit01/report.json'
    p=json.loads(oldprep.read_text()); g=json.loads(oldaudit.read_text())
    need(p['passed'] and g['passed'] and g['originalFaceOwnerExact'],'原全域准备未通过')
    for path,h in p['nativeInputSHA256'].items():need(sha(local(path))==h,'原全域输入改变')
    need(not CASE.exists() and not (BASE/'rocket_preparation01').exists(),'禁止覆盖候选准备')
    admission=dict(schema='RAH.HybridNKRocketAdmission.v1',profile=PROFILE,budget=BUDGET,
        temperatureRangeK=[200.,6000.],cells=64982,wallFaces=19752,savedSteps=list(range(0,121,20)),
        case=str(CASE),preparation=str(BASE/'rocket_preparation01/report.json'),record=str(BASE/'rocket_run01'),
        runtimeSHA256=sha(PACKAGE/'runtime.json'),evidenceSHA256={str(p):sha(p) for p in EVIDENCE},physicalQualification=False)
    check_admission(admission)
    CASE.mkdir();out=BASE/'rocket_preparation01';out.mkdir();shutil.copyfile(__file__,out/'executed_entry.py')
    for name in ['mesh.su2','surface.json','rah-air.gas','grid.lua']:shutil.copyfile(OLD/name,CASE/name)
    (CASE/'lmrsim').mkdir();shutil.copytree(OLD/'lmrsim/grid',CASE/'lmrsim/grid')
    for name in ['blocks.list','mpimap']:shutil.copyfile(OLD/'lmrsim'/name,CASE/'lmrsim'/name)
    text=(OLD/'job.lua').read_text(encoding='utf-8').split('config.gasdynamic_update_scheme=')[0]
    need(text.count('config.solver_mode="time-marching"')==1,'原生模式替换不明确')
    text=text.replace('config.solver_mode="time-marching"','config.solver_mode="steady"')
    (CASE/'job.lua').write_text(text+NK,encoding='utf-8')
    locks=[Path(__file__),ROOT/'scripts/hybridNKTetraControl.py',*EVIDENCE,oldprep,oldaudit,
           *[local(p) for p in p['nativeInputSHA256']],CASE/'job.lua']
    process=run(['lmr','prep-sim','--job=job.lua'],CASE,out/'native',90-(time.monotonic()-start),192,locks,4096,4096)
    need(process['passed'],'NK原生准备失败')
    cfg=json.loads((CASE/'lmrsim/config').read_text());nk=json.loads((CASE/'lmrsim/nk-config').read_text())
    controlcfg=json.loads((BASE/'tetra01/lmrsim/config').read_text());controlnk=json.loads((BASE/'tetra01/lmrsim/nk-config').read_text())
    need(nk==controlnk,'NK完整参数与已审核控制不相同')
    for key in ['solver_mode','viscous','reacting','turbulence_model','flux_calculator','interpolation_order',
                'spatial_deriv_calc','spatial_deriv_locn','thermo_interpolator','unstructured_limiter',
                'apply_unstructured_limiter_stagnation_point_filter','adjust_invalid_cell_data','max_invalid_cells']:
        need(cfg[key]==controlcfg[key],'NK实际算子与控制不相同: '+key)
    oldcfg=json.loads((OLD/'lmrsim/config').read_text())
    need(cfg['block_0']==oldcfg['block_0'],'原完整域或边界条件改变')
    for name in ['mesh.su2','surface.json','rah-air.gas','lmrsim/grid/grid-0000']:
        need(sha(CASE/name)==sha(OLD/name),'原完整几何/气体改变')
    f=read_field(CASE,0,64982);before=read_field(OLD,0,64982)
    need(set(f)==set(before),'初场列改变')
    need(all(np.array_equal(f[k],before[k],equal_nan=True) for k in f),'初场与原准备不相同')
    result=dict(passed=True,profile=PROFILE,geometry=m['geometry'],condition=m['condition'],process=process,
        nativeInputSHA256={str(p):sha(p) for p in CASE.rglob('*') if p.is_file()},
        sourceSHA256={str(p):sha(p) for p in locks},runtimeSHA256=admission['runtimeSHA256'],
        originalInitialFieldsExactlyEqual=True,originalGeometryAndBoundaryExactlyPreserved=True,
        elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    save(out/'report.json',result);save(BASE/'rocket_admission.json',admission)
    print(json.dumps(dict(passed=True,elapsedSeconds=result['elapsedSeconds'])))
    return result


def execute():
    """模块三：排他且有限的原生复数任务，保存场独立检查，不自动重试或接续。"""
    raise ValueError('full-wall-euler-nk-120-v1已完成唯一120更新但质量/能量收益筛选失败；候选关闭，禁止重跑、接续或同参数扫描')
    # 原准确执行源已保存rocket_run01/executed_entry.py；以下保留完整运行逻辑供追溯。
    from runEilmerProcess import run
    start=time.monotonic();ap=BASE/'rocket_admission.json';a=json.loads(ap.read_text());pkg=check_admission(a)
    prep=local(a['preparation']);r=json.loads(prep.read_text());out=local(a['record'])
    need(r['passed'] and r['profile']==PROFILE and not out.exists(),'NK准备不合格或已执行')
    for key in ['sourceSHA256','nativeInputSHA256']:
        for p,h in r[key].items():need(sha(local(p))==h,'NK准备后输入改变')
    with (CASE/'nk-attempt.json').open('x') as f:json.dump(dict(profile=PROFILE,admissionSHA256=sha(ap),automaticContinuation=0),f)
    out.mkdir();shutil.copyfile(__file__,out/'executed_entry.py')
    guard=SavedFieldGuard(CASE,64982,out/'saved_field_guard.json')
    locks=[Path(__file__),ROOT/'scripts/hybridRocketEuler.py',PACKAGE/'runtime.json',ap,prep,
           *pkg['filesSHA256'],*[local(p) for key in ['sourceSHA256','nativeInputSHA256'] for p in r[key]]]
    process=run([pkg['binary'],'--max-cpus=2','--max-wall-clock=600'],CASE,out/'native',
        630-(time.monotonic()-start),192,locks,4096,4096,checkpoint_guard=guard)
    result=dict(schema='RAH.HybridNKRocketRun.v1',profile=PROFILE,case=str(CASE),preparation=str(prep),
        admission=str(ap),admissionSHA256=sha(ap),processPassed=process['passed'],process=process,
        outputSHA256={str(p):sha(p) for p in (CASE/'lmrsim').rglob('*') if p.is_file()},
        elapsedSeconds=time.monotonic()-start,physicalQualification=False,newSolvedHeatAreaM2=0)
    save(out/'report.json',result)
    print(json.dumps(dict(processPassed=process['passed'],exit=process['nativeExitCode'],
        stop=process['stopReason'],elapsedSeconds=result['elapsedSeconds'])))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','execute']);a=p.parse_args()
    r=prepare() if a.stage=='prepare' else execute()
    raise SystemExit(0 if r.get('passed',r.get('processPassed',False)) else 1)
