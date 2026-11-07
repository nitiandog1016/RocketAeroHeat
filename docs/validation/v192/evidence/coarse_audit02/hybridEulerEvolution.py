# -*- coding: utf-8 -*-
"""有预算Euler发展合同；同场重启与细网格各具名，旧v191合同仍独立有效。"""
from pathlib import Path
import argparse
import json
import re
import shutil
import time
import numpy as np
from readEilmerNative import read_grid, read_field, compare_source, sha, need
from eilmerAirThermo import bind_model


def local_path(path):
    """跨Windows/WSL只转换已知E盘前缀，不推断其他路径。"""
    if Path('/mnt/e').exists(): return Path(str(path).replace('E:\\','/mnt/e/').replace('\\','/'))
    return Path(str(path).replace('/mnt/e/','E:/',1))


def audit_evolution_inputs(case, preparation):
    """模块一：新来源合同和实际算子逐项审核；不能由可调规模放行整箭旧失败。"""
    case=Path(case);r=json.loads(Path(preparation).read_text())
    need(r.get('schema')=='RAH.HybridEulerEvolutionPreparation.v1' and r.get('passed'),'发展准备未通过')
    need(local_path(r['case']).resolve()==case.resolve(),'发展目录身份错误')
    need(r['profile'] in ('coarse-restart-400-to-2000','refined-120-2-40'),'未登记的发展方案')
    for p,d in r['nativeInputSHA256'].items():need(sha(local_path(p))==d,'发展输入已变: '+p)
    for p,d in r['sourceSHA256'].items():need(sha(local_path(p))==d,'发展原始来源已变: '+p)
    cfg=json.loads((case/'lmrsim/config').read_text());ctl=json.loads((case/'lmrsim/control').read_text())
    baseline=json.loads(local_path(r['baselineConfig']).read_text())
    # 两种网格仅block计数/边界局部几何可不同；其余完整全局配置逐项一致。
    need({k:v for k,v in cfg.items() if not k.startswith('block_')}==
         {k:v for k,v in baseline.items() if not k.startswith('block_')},'全局实际算子/边界气体变化')
    need(ctl==r['expectedControl'],'实际控制与事前登记不符')
    geometry=r['geometry'];grid=read_grid(case/'lmrsim/grid/grid-0000');mapping,error=compare_source(grid,case/'mesh.su2')
    dims=geometry['gridIntervals'];nx,ny,nz=dims
    need(dims in ([60,2,20],[120,2,40]),'超出具名小控制网格范围')
    counts=(nx*ny*nz,(nx+1)*(ny+1)*(nz+1),nx*ny*(nz+1)+nx*(ny+1)*nz+(nx+1)*ny*nz)
    need((len(grid['cells']),len(grid['points']),len(grid['faces']))==counts,'原生网格规模错误')
    for name,b in grid['boundaries'].items():
        actual=cfg['block_0']['boundary_%d'%b['index']]
        expected=baseline['block_0']['boundary_%d'%b['index']]
        need(actual==expected,'实际边界变化: '+name)
    f=read_field(case,0,len(grid['cells']));volume=float(f['vol'].sum())
    need(np.all(f['vol']>0) and abs(volume-.04*(.6-.24*np.tan(np.radians(10))))<2e-13,'体积不符')
    model=bind_model(case/'rah-air.gas',sha(case/'rah-air.gas'))
    return dict(passed=True,cells=counts[0],vertices=counts[1],faces=counts[2],
        ghostCells=sum(len(b['ids']) for b in grid['boundaries'].values()),volumeM3=volume,
        sourceCoordinateError=error,actualGasModel=model,physicalQualification=False),grid,mapping


def prepare_restart(source, source_preparation, source_run, case, record):
    """模块二：保留真实保存2/400步与原历史，只在新目录改变终止/保存设置。"""
    from readHybridEulerNative import read_wall_frame
    source,case,record=map(Path,(source,case,record))
    need(not case.exists() and not record.exists(),'禁止覆盖发展目录')
    read_wall_frame(source,2,source_preparation,source_run)
    metadata=(source/'lmrsim/snapshots/snapshot-times-metadata').read_text()
    entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)",metadata,re.M)
    need([(int(i),int(s)) for i,s,t in entries]==[(0,0),(1,200),(2,400)],'旧帧不是400步终态')
    originals=[p for p in source.rglob('*') if p.is_file()]
    binding={str(p.resolve()):sha(p) for p in originals}
    shutil.copytree(source,case);record.mkdir(parents=True)
    # 原尝试标记逐字保留，只供来源审计；本轮用独立evolution-attempt防重复。
    ctl=json.loads((case/'lmrsim/control').read_text());ctl.update(max_step=2000,write_flow_solution_at_step=400)
    (case/'lmrsim/control').write_text(json.dumps(ctl,indent=2)+'\n')
    design=json.loads((case/'manifest.json').read_text());design['nativeConfig'].update(max_step=2000,write_flow_solution_at_step=400)
    design.update(candidate='coarse-restart-400-to-2000',restartSnapshot=2,restartStep=400,executionAuthorized=False)
    (case/'manifest.json').write_text(json.dumps(design,indent=2)+'\n')
    job=(case/'job.lua').read_text().replace('config.max_step=400','config.max_step=2000').replace('config.write_flow_solution_at_step=200','config.write_flow_solution_at_step=400')
    (case/'job.lua').write_text(job,encoding='utf8')
    design['filesSHA256']={n:sha(case/n) for n in design['filesSHA256']}
    (case/'manifest.json').write_text(json.dumps(design,indent=2)+'\n')
    # times/history会原生追加，不纳入本地不变输入；原目录逐字绑定仍必须稳定。
    mutable={'lmrsim/snapshots/snapshot-times-metadata','lmrsim/transient-residuals.dat','lmrsim/progress-time-marching'}
    inputs={str(p.resolve()):sha(p) for p in case.rglob('*') if p.is_file() and p.relative_to(case).as_posix() not in mutable}
    r=dict(schema='RAH.HybridEulerEvolutionPreparation.v1',profile='coarse-restart-400-to-2000',passed=True,
        case=str(case.resolve()),geometry=design['geometry'],expectedControl=ctl,
        baselineConfig=str((source/'lmrsim/config').resolve()),sourceSHA256=binding,nativeInputSHA256=inputs,
        sourcePreparation=str(Path(source_preparation).resolve()),sourceRun=str(Path(source_run).resolve()),
        sourceSnapshot=2,sourceStep=400,sourceTime=float(entries[-1][2]),
        historyPrefixSHA256=sha(case/'lmrsim/transient-residuals.dat'),physicalQualification=False)
    for p in [Path(source_preparation),Path(source_run)]:r['sourceSHA256'][str(p.resolve())]=sha(p)
    (record/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    shutil.copyfile(__file__,record/'executed_prepare.py')
    audit_evolution_inputs(case,record/'report.json')
    return r


def execute(admission, record):
    """模块三：精确准入、排他启动和固定运行包；不自动续算。"""
    from runEilmerProcess import run
    from eilmerWorkspaceRuntime import audit_workspace_runtime
    start=time.monotonic();admission,record=map(Path,(admission,record));a=json.loads(admission.read_text())
    need(a['schema']=='RAH.HybridEulerEvolutionAdmission.v1','发展准入版本不符')
    for p,d in a['evidenceSHA256'].items():need(sha(local_path(p))==d,'发展准入证据已变')
    case=Path(a['case']);prep=Path(a['preparation']);r=json.loads(prep.read_text());budget=a['budget']
    allowed={'coarse-restart-400-to-2000':(120,135,2000,2),'refined-120-2-40':(240,255,1600,0)}
    need(r['profile'] in allowed,'未登记方案');seconds,entry,steps,snapshot=allowed[r['profile']]
    need(budget==dict(nativeSeconds=seconds,entrySeconds=entry,maxStep=steps,memoryMiB=4096,
                     outputMiB=64,cpuCount=2,automaticContinuation=0),'发展预算不符')
    need(a['record']==str(record.resolve()) and not record.exists(),'发展回执目录错误')
    audit_evolution_inputs(case,prep)
    root=Path(__file__).resolve().parents[1];runtime=root/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json'
    package,locks=audit_workspace_runtime(dict(explicit=dict(profile='rah-explicit-workspace-v1',
        runtimePackage=str(runtime),runtimePackageSHA256=a['runtimeSHA256'])))
    with (case/'evolution-attempt.json').open('x') as f:json.dump(dict(admissionSHA256=sha(admission)),f)
    record.mkdir(parents=True);shutil.copyfile(__file__,record/'executed_entry.py')
    modules=[Path(__file__).with_name(n) for n in ['readHybridEulerNative.py','eilmerAirThermo.py','readEilmerNative.py','hybridExternalState.py']]
    for p in modules:shutil.copyfile(p,record/p.name)
    locks+=list(r['nativeInputSHA256'])+list(r['sourceSHA256'])+list(a['evidenceSHA256'])+modules+[Path(__file__),prep,admission,case/'evolution-attempt.json']
    process=run([package['binary'],'--max-cpus=2','--snapshot-start='+str(snapshot),'--max-wall-clock='+str(seconds)],
        case,record/'native',entry-(time.monotonic()-start),64,locks,4096,4096)
    result=dict(schema='RAH.HybridEulerControlRun.v1',profile=r['profile'],processPassed=process['passed'],
        elapsedSeconds=time.monotonic()-start,admissionSHA256=sha(admission),physicalQualification=False,
        outputSHA256={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()})
    (record/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(processPassed=process['passed'],elapsedSeconds=result['elapsedSeconds'])))
    return process['passed']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['restart','run'])
    for k in ['source','source-preparation','source-run','case','record','admission']:p.add_argument('--'+k)
    a=p.parse_args()
    if a.mode=='restart':prepare_restart(a.source,a.source_preparation,a.source_run,a.case,a.record)
    else:raise SystemExit(0 if execute(a.admission,a.record) else 1)
