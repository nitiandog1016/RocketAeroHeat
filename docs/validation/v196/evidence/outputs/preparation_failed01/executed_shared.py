# -*- coding: utf-8 -*-
"""真实全壁Euler有限发展、保存场守卫和逐原面读取；无热资格自动升级。"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import time
import numpy as np
from readEilmerNative import sha, need, read_grid, read_field, compare_source
from prepareHybridRocketEulerNative import local, verify_design, save
from eilmerAirThermo import bind_model, evaluate
from hybridEdgeSampling import map_original_wall, candidate
from hybridExternalState import exchange, fingerprint

ROOT=Path(__file__).resolve().parents[1]
PROFILE='full-wall-euler-400-v1'
GEOMETRIC_PROFILE='full-wall-euler-geometric-400-v1'
BUDGET=dict(nativeSeconds=600,entrySeconds=630,parentSeconds=660,maxStep=400,cpuCount=2,
            memoryMiB=4096,outputMiB=128,automaticContinuation=0)


def check_admission(a):
    """模块一：具名总预算和固定任务身份，不放宽旧小例/黏性合同。"""
    need(a['schema']=='RAH.HybridRocketEulerAdmission.v1' and a['profile']==PROFILE,'真实发展准入身份错误')
    need(a['budget']==BUDGET,'真实发展总预算错误')
    need(a['temperatureRangeK']==[200.,6000.],'物性接口温域改变')
    for p,d in a['evidenceSHA256'].items():need(sha(local(p))==d,'准入来源变化: '+p)


def prepare(source, original_receipt, original_audit, case, record, profile=PROFILE):
    """模块二：复制已审核原生输入，只改已登记三处运行/保存参数，不重新prep。"""
    need(profile in (PROFILE,GEOMETRIC_PROFILE),'未知真实发展准备身份')
    source,original_receipt,original_audit,case,record=map(lambda p:Path(p).resolve(),
        (source,original_receipt,original_audit,case,record))
    need(not case.exists() and not record.exists(),'禁止覆盖真实发展目录')
    m=verify_design(source);r=json.loads(original_receipt.read_text());audit=json.loads(original_audit.read_text())
    need(r['passed'] and r['schema']=='RAH.HybridRocketEulerNativePreparation.v1' and audit['passed'],'原生准备未审核')
    need(local(r['case']).resolve()==source and audit['sourceSHA256'][str(original_receipt)]==sha(original_receipt),
         '原生准备/审核来源不符')
    for p,d in r['nativeInputSHA256'].items():need(sha(local(p))==d,'原准备输入已变')
    sourcefiles=[p for p in source.rglob('*') if p.is_file()]+[original_receipt,original_audit]
    hashes={str(p):sha(p) for p in sourcefiles}
    shutil.copytree(source,case);record.mkdir(parents=True)
    cfg=json.loads((source/'lmrsim/config').read_text());cfg['dt_plot_schedule_values']=[.0002]
    ctl=json.loads((source/'lmrsim/control').read_text());ctl.update(max_step=400,write_flow_solution_at_step=200)
    save(case/'lmrsim/config',cfg);save(case/'lmrsim/control',ctl)
    files={str(p):sha(p) for p in case.rglob('*') if p.is_file()}
    result=dict(schema='RAH.HybridRocketEulerDevelopmentPreparation.v1',profile=profile,passed=True,
        case=str(case),originalCase=str(source),originalPreparation=str(original_receipt),originalAudit=str(original_audit),
        sourceSHA256=hashes,nativeInputSHA256=files,geometry=m['geometry'],condition=m['condition'],
        actualOverrides={'control':dict(max_step=400,write_flow_solution_at_step=200),
                         'config':dict(dt_plot_schedule_values=[.0002])},
        luaIdentity='unchanged original preparation only; development uses exactly named JSON overrides',
        executionAuthorized=False,physicalQualification=False)
    shutil.copyfile(__file__,record/'executed_prepare.py');save(record/'report.json',result)
    check_inputs(case,record/'report.json')
    return result


def check_inputs(case, preparation):
    """原基线与复制品SHA逐一核对，并验证完整配置只发生已准入变化。"""
    case=Path(case).resolve();r=json.loads(Path(preparation).read_text())
    need(r['schema']=='RAH.HybridRocketEulerDevelopmentPreparation.v1' and r['profile'] in (PROFILE,GEOMETRIC_PROFILE) and r['passed'],
         '真实发展准备未通过')
    need(local(r['case']).resolve()==case,'真实发展目录不符')
    for key in ['sourceSHA256','nativeInputSHA256']:
        for p,d in r[key].items():need(sha(local(p))==d,'真实发展来源改变: '+p)
    original=local(r['originalCase']);m=verify_design(original)
    need(m['geometry']==r['geometry'] and m['condition']==r['condition'],'几何/工况合同改变')
    cfg=json.loads((original/'lmrsim/config').read_text());cfg['dt_plot_schedule_values']=[.0002]
    ctl=json.loads((original/'lmrsim/control').read_text());ctl.update(max_step=400,write_flow_solution_at_step=200)
    need(json.loads((case/'lmrsim/config').read_text())==cfg,'真实发展原生配置不符')
    need(json.loads((case/'lmrsim/control').read_text())==ctl,'真实发展原生控制不符')
    return r


class SavedFieldGuard:
    """模块三：只观察完整稳定文件；记录首次完整观察，不声称瞬时越界时刻。"""
    def __init__(self,case,cells,report,maximum_dt=None):
        self.case=Path(case);self.cells=cells;self.report=Path(report);self.seen={};self.rows=[]
        self.maximum_dt=maximum_dt
        self.metadata=self.case/'lmrsim/snapshots/fluid.metadata'
        text=self.metadata.read_text();need('revision-id: f53f460' in text,'守卫场版本不符')
        self.metadataSHA=sha(self.metadata)
        self.names=re.findall(r"^\s*- '([^']+)'\s*$",text,re.M)
        need(len(self.names)==len(set(self.names)) and all(k in self.names for k in ['p','rho','T','pos.x','pos.y','pos.z']),
             '守卫字段缺失')
        self.expected=len(self.names)*cells*8

    def __call__(self):
        need(sha(self.metadata)==self.metadataSHA,'守卫元数据改变')
        for path in sorted((self.case/'lmrsim/snapshots').glob('[0-9][0-9][0-9][0-9]/fluid-0000')):
            stat=path.stat();identity=(stat.st_size,stat.st_mtime_ns)
            if stat.st_size!=self.expected or self.seen.get(str(path))==identity:continue
            raw=path.read_bytes();after=path.stat()
            if len(raw)!=self.expected or (after.st_size,after.st_mtime_ns)!=identity:continue
            fields=dict(zip(self.names,np.frombuffer(raw,dtype='<f8').reshape(len(self.names),self.cells)))
            T=fields['T'];bad=~np.isfinite(T)|(T<200)|(T>6000)
            for name in self.names:
                if name=='dt_local':continue
                bad|=~np.isfinite(fields[name])
            for name in ['p','rho','vol']:bad|=fields[name]<=0
            ids=np.flatnonzero(bad);finite=T[np.isfinite(T)]
            row=dict(snapshot=int(path.parent.name),path=str(path),sha256=hashlib.sha256(raw).hexdigest(),
                temperatureRangeK=[float(finite.min()),float(finite.max())] if len(finite) else None,
                outsideTemperatureCount=int(np.sum(~np.isfinite(T)|(T<200)|(T>6000))),invalidCells=int(len(ids)),
                firstInvalidCellZeroBased=int(ids[0]) if len(ids) else None,
                kind='first-observation-of-complete-saved-field; not first-instantaneous-crossing')
            if len(ids):
                i=ids[0];row['firstInvalidPositionM']=[float(fields['pos.'+k][i]) if np.isfinite(fields['pos.'+k][i]) else None for k in 'xyz']
            self.rows.append(row);self.seen[str(path)]=identity
            save(self.report,dict(rows=self.rows,metadataSHA256=self.metadataSHA,physicalQualification=False))
            if len(ids):return row
        # 首次候选暴露了有限场之外的伪时间异常；只检查完整已观察帧的元数据。
        # 不改变推进，不将修正后的守卫追称为首次运行已启用。
        times=self.case/'lmrsim/snapshots/snapshot-times-metadata'
        if self.maximum_dt is not None and times.exists():
            entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)\s*\n\s+dt:\s*(\S+)",times.read_text(),re.M)
            completed={row['snapshot'] for row in self.rows}
            for index,step,t,dt in entries:
                if int(index) not in completed:continue
                values=np.array([float(t),float(dt)])
                if not np.isfinite(values).all() or values[0]<0 or not 0<values[1]<=self.maximum_dt:
                    rejection=dict(kind='invalid-saved-pseudo-time',snapshot=int(index),step=int(step),
                        pseudoTime=float(t) if np.isfinite(values[0]) else None,
                        dt=float(dt) if np.isfinite(values[1]) else None,maximumDt=self.maximum_dt)
                    save(self.report,dict(rows=self.rows,timeRejection=rejection,metadataSHA256=self.metadataSHA,physicalQualification=False))
                    return rejection
        return None


def execute(admission, record):
    """模块四：唯一排他任务，原生进程组/RSS/输出及完整场守卫共用预算。"""
    raise ValueError('full-wall-euler-400-v1已因首步伪时间3e98关闭；必须先核验实际四面体时间尺度并建立新具名合同，禁止重跑本候选')
    # 以下保留原已登记流程作来源；当前在上方明确拒绝，没有可调用的绕过入口。
    from runEilmerProcess import run
    from eilmerWorkspaceRuntime import audit_workspace_runtime
    start=time.monotonic();admission,record=map(lambda p:Path(p).resolve(),(admission,record))
    a=json.loads(admission.read_text());check_admission(a);case=local(a['case']);prep=local(a['preparation'])
    need(local(a['record']).resolve()==record and not record.exists(),'真实发展回执路径错误')
    r=check_inputs(case,prep)
    runtime=ROOT/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json'
    package,locks=audit_workspace_runtime(dict(explicit=dict(profile='rah-explicit-workspace-v1',
        runtimePackage=str(runtime),runtimePackageSHA256=a['runtimeSHA256'])))
    with (case/'rocket-euler-attempt.json').open('x') as f:json.dump(dict(admissionSHA256=sha(admission),automaticContinuation=0),f)
    record.mkdir(parents=True);shutil.copyfile(__file__,record/'executed_entry.py')
    guard=SavedFieldGuard(case,r['geometry']['cells'],record/'saved_field_guard.json',maximum_dt=.001)
    locks += [local(p) for key in ['sourceSHA256','nativeInputSHA256'] for p in r[key]]
    locks += [local(p) for p in a['evidenceSHA256']]+[admission,prep,Path(__file__),case/'rocket-euler-attempt.json']
    process=run([package['binary'],'--max-cpus=2','--max-wall-clock=600'],case,record/'native',
        630-(time.monotonic()-start),128,locks,4096,4096,checkpoint_guard=guard)
    output={str(p):sha(p) for p in (case/'lmrsim').rglob('*') if p.is_file()}
    result=dict(schema='RAH.HybridRocketEulerRun.v1',profile=PROFILE,processPassed=process['passed'],
        case=str(case),preparation=str(prep),admissionSHA256=sha(admission),admission=str(admission),
        outputSHA256=output,physicalQualification=False,newSolvedHeatAreaM2=0,
        elapsedSeconds=time.monotonic()-start)
    save(record/'report.json',result)
    print(json.dumps(dict(processPassed=process['passed'],stopReason=process['stopReason'],
                         nativeExitCode=process['nativeExitCode'],elapsedSeconds=result['elapsedSeconds'])))
    return result


def validate_terminal(case,preparation,run_receipt):
    """模块五：生产读取拒绝异常/守卫停止/错包和旧帧；单独诊断不得绕过此入口。"""
    case=Path(case).resolve();rp=Path(run_receipt);record=json.loads(rp.read_text())
    if record.get('profile')==GEOMETRIC_PROFILE:
        from geometricRocketEuler import validate_terminal as validate_geometric
        return validate_geometric(case,preparation,rp)
    process=json.loads((rp.parent/'native/report.json').read_text())
    need(record['schema']=='RAH.HybridRocketEulerRun.v1' and record['profile']==PROFILE and record['processPassed']
         and process['passed'] and process['nativeExitCode']==0 and process['inputsStable']
         and process['stopReason'] is None and not process['liveGroupAfterCleanup'],'异常真实Euler结果不准生产导出')
    need(local(record['case']).resolve()==case and local(process['cwd']).resolve()==case,'运行结果目录错配')
    need(local(record['preparation']).resolve()==Path(preparation).resolve(),'运行准备身份不符')
    ap=local(record['admission']);need(sha(ap)==record['admissionSHA256'],'准入记录改变')
    a=json.loads(ap.read_text());check_admission(a)
    runtime=ROOT/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json';pkg=json.loads(runtime.read_text())
    need(sha(runtime)==a['runtimeSHA256'] and process['argv']==[pkg['binary'],'--max-cpus=2','--max-wall-clock=600']
         and process['inputSHA256'].get(pkg['binary'])==pkg['filesSHA256'][pkg['binary']], '实际求解包不符')
    for p,d in record['outputSHA256'].items():need(sha(local(p))==d,'保存结果已变: '+p)
    check_inputs(case,preparation)
    return record,process


def read_wall(case,snapshot,preparation,run_receipt=None):
    """模块六：原三角面精确读取同EOS壁邻状态；不给初场或有限值热资格。"""
    need(type(snapshot) is int and snapshot>=0,'帧编号不符')
    case=Path(case).resolve();r=check_inputs(case,preparation)
    extra=[]
    if snapshot>0:
        need(run_receipt is not None,'发展帧缺运行来源')
        validate_terminal(case,preparation,run_receipt)
        extra=[Path(run_receipt),Path(run_receipt).parent/'native/report.json']
    grid=read_grid(case/'lmrsim/grid/grid-0000');mapping,_=compare_source(grid,case/'mesh.su2')
    f=read_field(case,snapshot,r['geometry']['cells']);need(np.all(f['vol']>0),'原生非正体积')
    gas=evaluate(f['p'],f['T']);errors={}
    for name,key in [('rho','densityKgm3'),('e','internalEnergyJkg'),('a','soundSpeedMps')]:
        errors[name]=float(np.max(abs(f[name]-gas[key])/np.maximum(1,abs(gas[key]))))
        need(errors[name]<2e-10,'全场EOS不相容: '+name)
    need(np.all(f['massf-air']==1),'纯空气身份改变')
    surface=json.loads((case/'surface.json').read_text());position=np.column_stack([f['pos.'+k] for k in 'xyz'])
    location=map_original_wall(surface,grid,mapping['wall'],position);owners=np.array(location['nativeOwnerCellsZeroBased'])
    model=bind_model(case/'rah-air.gas',sha(case/'rah-air.gas'));condition=r['condition']
    fields={key:gas[key][owners].tolist() for key in ['pressurePa','temperatureK','densityKgm3','internalEnergyJkg','enthalpyJkg','entropyJkgK']}
    U=np.column_stack([f['vel.'+k][owners] for k in 'xyz']);fields['velocityMps']=U.tolist()
    fields['totalEnthalpyJkg']=(gas['enthalpyJkg'][owners]+np.sum(U*U,axis=1)/2).tolist()
    V=np.array(surface['vertices']);F=np.array(surface['facesOneBased'])-1;tri=V[F]
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    sources=[case/'lmrsim/snapshots'/('%04d'%snapshot)/'fluid-0000',case/'lmrsim/snapshots/fluid.metadata',
        case/'lmrsim/grid/grid-0000',case/'lmrsim/config',case/'lmrsim/control',case/'surface.json',case/'rah-air.gas',
        Path(preparation),Path(__file__),Path(__file__).with_name('eilmerAirThermo.py'),
        Path(__file__).with_name('hybridEdgeSampling.py'),Path(__file__).with_name('hybridExternalState.py')]+extra
    hashes={str(p.resolve()):sha(p) for p in sources}
    frame=dict(schema='RAH.HybridEulerFrame.v1',solver='Euler',turbulenceModel='none',
        surfaceSHA256=fingerprint(surface),conditionSHA256=fingerprint(condition),modelSHA256=fingerprint(model),
        sourceFilesSHA256=hashes,snapshot=snapshot,faceIDs=surface['faceIDs'],originalAreaM2=float(area.sum()),
        faceCentersM=tri.mean(1).tolist(),fields=fields)
    exchanged=exchange(frame,surface,condition,model,hashes,evaluate)
    result=candidate(exchanged,surface,grid,mapping,f)
    result.update(fullVolumeEOSErrors=errors,stateIdentity='uniform-initial-not-solved' if snapshot==0 else 'unqualified-Euler-development',
                  faceAreasM2=area.tolist(),fullVolumeTemperatureRangeK=[float(f['T'].min()),float(f['T'].max())],
                  numericsProfile=r['profile'],physicalQualification=False)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='mode',required=True)
    for name,keys in [('prepare',['source','original-receipt','original-audit','case','record']),('execute',['admission','record'])]:
        s=sub.add_parser(name)
        for key in keys:s.add_argument('--'+key,required=True)
    a=vars(p.parse_args());mode=a.pop('mode');result=globals()[mode](**a)
    if mode=='prepare':print(json.dumps(dict(passed=result['passed'])))
    raise SystemExit(0 if result.get('passed',result.get('processPassed',False)) else 1)
