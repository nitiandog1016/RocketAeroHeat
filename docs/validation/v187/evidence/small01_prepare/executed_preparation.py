"""从准确准备包生成具名发展输入；复用原生网格，不重跑网格器、不迁移旧流场。"""
from pathlib import Path
import argparse,copy,json,os,shutil,time
from prepareEilmerRocketCase import validate_condition,make_job
from runEilmerProcess import run,sha,COMMIT
from readEilmerNative import need,audit_numerics,audit_run_scope
from eilmerExplicitProfile import audit_explicit_config,settings


def prepare(source,prepared,audit,case,record,condition,seconds,output_mib=1024):
    start=time.monotonic();source=Path(source).resolve();case=Path(case).resolve();record=Path(record).resolve()
    need(not case.exists() and not record.exists(),'拒绝覆盖')
    m=json.loads((source/'case_manifest.json').read_text());prior=json.loads(Path(prepared).read_text())
    need(m['sourceCommit']==COMMIT and prior['passed'] and prior['stage']=='prepare' and prior['case']==str(source),'源准备身份不符')
    need(prior['caseManifestSHA256']==sha(source/'case_manifest.json'),'源清单变化')
    geometry=json.loads(Path(audit).read_text())
    need(geometry['passed'] and geometry['originalCellOrderAndBoundaryBijection'] and
         geometry['sourceMeshSHA256']==m['mesh']['sourceSHA256'] and geometry['nativeGridSHA256']==sha(source/'lmrsim/grid/grid-0000'),
         '缺少相容的既有全域几何证据')
    # 只复用已绑定的配置/网格/初态；源对象可已运行，但其已审核初场必须逐字保持。
    for p,h in prior['nativeInputSHA256'].items():need(sha(p)==h,'源准备输入变化: '+p)
    c=validate_condition(condition);old=copy.deepcopy(m['condition'])
    changed=copy.deepcopy(c)
    if 'explicit' in c:
        changed.pop('explicit');old.pop('development',None);old.pop('maxNewtonSteps')
        need(changed==old,'显式入口不得修改原物理工况')
    else:
        changed.pop('development');changed['maxNewtonSteps']=old['maxNewtonSteps']
        need(changed==old and 'development' in c,'发展入口不得修改原物理工况')
    scoped=dict(condition=c,mesh=m['mesh']);audit_run_scope(scoped,case)
    extra_locks=[]
    if c.get('explicit',{}).get('profile')=='rah-explicit-workspace-v1':
        from eilmerWorkspaceRuntime import audit_workspace_runtime
        _,extra_locks=audit_workspace_runtime(c)
    if c.get('explicit',{}).get('profile') in ('rah-explicit-dense-v1','rah-explicit-dense-development-v1'):
        from eilmerDenseRuntime import audit_dense_runtime
        _,extra_locks=audit_dense_runtime(c)
    if c.get('explicit',{}).get('profile') in ('rah-explicit-local-development-v1','rah-explicit-dense-development-v1'):
        from eilmerExplicitAdmission import audit_admission
        admission,admission_locks=audit_admission(scoped,case);extra_locks+=admission_locks
        need(seconds==admission['budget']['prepareSeconds'] and output_mib==admission['budget']['outputMiB'],'准备预算不符')
        need(str(source)==admission['sourceCase'] and str(Path(prepared).resolve())==admission['evidence']['sourcePreparation']['path'] and
             str(Path(audit).resolve())==admission['evidence']['sourceGeometry']['path'],'准备入口不是登记来源')
        extra_locks.append(Path(__file__).with_name('eilmerExplicitAdmission.py'))
    need(sha(m['mesh']['nativeSourcePath'])==m['mesh']['sourceSHA256'],'原网格变化')
    case.mkdir(parents=True);record.mkdir(parents=True)
    shutil.copyfile(__file__,record/'executed_preparation.py')
    for name in ['gas.inp','grid.lua','rah-air.gas']:
        shutil.copyfile(source/name,case/name)
    # 仅硬链接不会被prep-sim写入的只读网格；新初场始终由原生单独写入。
    for p in (source/'lmrsim/grid').iterdir():
        dest=case/'lmrsim/grid'/p.name;dest.parent.mkdir(parents=True,exist_ok=True);os.link(p,dest)
    shutil.copyfile(source/'lmrsim/blocks.list',case/'lmrsim/blocks.list')
    (case/'job.lua').write_text(make_job(m['mesh']['nativeSourcePath'],c),encoding='utf8',newline='\n')
    m['condition']=c;m['status']='development-input';m['pending']=['native-phase-check','nonlinear-development','strict-heat-and-physical-validation']
    m['reusedPreparation']=dict(case=str(source),receipt=str(Path(prepared).resolve()),receiptSHA256=sha(prepared),
        geometryAudit=str(Path(audit).resolve()),geometryAuditSHA256=sha(audit))
    m['filesSHA256']={n:sha(case/n) for n in ['gas.inp','grid.lua','rah-air.gas','job.lua']}
    (case/'case_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    locks=[source/'case_manifest.json',Path(prepared),Path(audit),Path(__file__),case/'case_manifest.json',case/'job.lua',
           Path(__file__).with_name('prepareEilmerRocketCase.py'),Path(__file__).with_name('readEilmerNative.py')]
    locks += [Path(p) for p in prior['nativeInputSHA256']]
    locks += extra_locks
    if 'explicit' in c:locks.append(Path(__file__).with_name('eilmerExplicitProfile.py'))
    result=dict(stage='prepare',case=str(case),passed=False,caseManifestSHA256=sha(case/'case_manifest.json'),commands=[],physicsQualified=False,newQualifiedAreaM2=0)
    try:
        process=run(['lmr','prep-sim','--job=job.lua'],case,record/'process00',seconds-(time.monotonic()-start),output_mib,locks)
        result['commands'].append(process);need(process['passed'],'原生准备失败')
        if 'explicit' in c:
            cfg=json.loads((case/'lmrsim/config').read_text())
            audit_explicit_config(cfg,c,json.loads((case/'lmrsim/control').read_text()))
            previous_cfg=json.loads((source/'lmrsim/config').read_text())
            allowed=set(settings(c['explicit']))|{p+s for p in ['cfl','dt_plot'] for s in ['_schedule_length','_schedule_times','_schedule_values']}
            delta={key for key in set(cfg)|set(previous_cfg) if cfg.get(key)!=previous_cfg.get(key)}
            need(delta<=allowed,'显式配置修改了未登记的物理或空间离散字段')
            result['nativeConfigChangedKeys']=sorted(delta)
        else:
            audit_numerics(json.loads((case/'lmrsim/nk-config').read_text()),c)
            need((case/'lmrsim/config').read_bytes()==(source/'lmrsim/config').read_bytes(),'非NK物理配置改变')
        for name in ['native-input-state.json','lmrsim/snapshots/0000/fluid-0000','lmrsim/grid/grid-0000']:
            need(sha(case/name)==sha(source/name),'原工况/网格/初场不等价: '+name)
        result.update(passed=True,uniformInitialFieldBitwiseIdentical=True,nativeGridReused=True,
            nativeInputSHA256={str(p):sha(p) for p in case.rglob('*') if p.is_file()},
            reusedPreparationReceiptSHA256=sha(prepared))
    except Exception as e:result['error']=repr(e)
    result['elapsedSeconds']=time.monotonic()-start
    (record/'report.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--prepared',required=True)
    p.add_argument('--case',required=True);p.add_argument('--record',required=True);p.add_argument('--condition',required=True);p.add_argument('--seconds',type=float,required=True);p.add_argument('--audit',required=True)
    p.add_argument('--output-mib',type=int,default=1024)
    a=p.parse_args();r=prepare(a.source,a.prepared,a.audit,a.case,a.record,json.loads(Path(a.condition).read_text()),a.seconds,a.output_mib)
    print(json.dumps({k:r.get(k) for k in ['passed','elapsedSeconds','error']}));raise SystemExit(0 if r['passed'] else 1)
