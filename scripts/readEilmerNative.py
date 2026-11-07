"""独立读取固定Eilmer原生几何/初场/壁面载荷；保留原SU2面身份和入壁热流正号。"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re
import struct
import time
import numpy as np
from auditEilmerCompletion import audit_completion,completion_log_for_receipt

TYPES={4:10,5:13,6:12,7:14}
WIDTH={10:4,13:6,12:8,14:5}


def need(ok,message):
    if not ok: raise ValueError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def read_grid(path):
    """模块一：按固定源码rawWrite合同读取，而不是借用上游网格解释器。"""
    data=Path(path).read_bytes();magic=b'unstructured_grid 1.0'
    need(data.startswith(magic),'非Eilmer原生三维网格');pos=len(magic)
    def integer():
        nonlocal pos
        need(pos+4<=len(data),'原生网格截断')
        n=struct.unpack_from('<i',data,pos)[0];pos+=4;return n
    def array():
        nonlocal pos
        n=integer();need(0<=n<=(len(data)-pos)//4,'原生数组长度无效')
        value=np.frombuffer(data,dtype='<i4',count=n,offset=pos);pos+=4*n;return value
    def string():
        nonlocal pos
        n=integer();need(0<=n<=len(data)-pos,'原生字符串截断')
        value=data[pos:pos+n].decode('utf8');pos+=n;return value
    label=string();dim,nv,nf,nc,nb=[integer() for _ in range(5)]
    need(dim==3 and min(nv,nf,nc,nb)>0,'当前只支持非空三维混合网格')
    need(24*nv+8*nf+32*nc<=len(data),'原生网格计数与大小不符')
    points=np.frombuffer(data,dtype='<f8',count=3*nv,offset=pos).reshape(nv,3);pos+=24*nv
    need(np.isfinite(points).all(),'原生坐标非有限')
    faces=np.full((nf,4),-1,dtype=np.int32)
    for i in range(nf):
        ids=array();need(len(ids) in (3,4) and min(ids)>=0 and max(ids)<nv,'原生面节点无效')
        need(len(set(ids))==len(ids),'原生面重复节点');faces[i,:len(ids)]=ids
    cells=np.full((nc,8),-1,dtype=np.int32);kinds=np.zeros(nc,dtype=np.int8)
    owners=np.full((nf,2),-1,dtype=np.int32)
    for i in range(nc):
        kind=TYPES.get(integer());need(kind is not None,'原生单元类型不支持')
        ids,face_ids,signs=array(),array(),array()
        need(len(ids)==WIDTH[kind] and min(ids)>=0 and max(ids)<nv,'原生单元节点无效')
        need(len(face_ids)==len(signs) and len(set(face_ids))==len(face_ids),'单元面数/重复错误')
        need(min(face_ids)>=0 and max(face_ids)<nf and np.all(np.abs(signs)==1),'面ID/符号无效')
        side=(signs<0).astype(np.int8)
        need(np.all(owners[face_ids,side]==-1),'同向重复单元面')
        owners[face_ids,side]=i;cells[i,:len(ids)]=ids;kinds[i]=kind
    boundaries={};is_boundary=np.zeros(nf,dtype=bool)
    for ib in range(nb):
        tag=string();ids,signs=array(),array()
        need(tag not in boundaries and len(ids)==len(signs)>0,'边界重复或长度不符')
        need(min(ids)>=0 and max(ids)<nf and np.all(np.abs(signs)==1),'边界ID/符号无效')
        need(len(set(ids))==len(ids) and not np.any(is_boundary[ids]),'边界重复面')
        side=(signs<0).astype(np.int8)
        need(np.all(owners[ids,side]>=0) and np.all(owners[ids,1-side]==-1),'边界外向符号/归属错误')
        is_boundary[ids]=True
        boundaries[tag]=dict(index=ib,ids=ids,signs=signs,owners=owners[ids,side])
    need(pos==len(data),'原生网格多余/未解析内容')
    need(np.array_equal(is_boundary,np.sum(owners>=0,axis=1)==1),'边界集合不等于全外边界')
    need(np.all(np.sum(owners>=0,axis=1)>=1),'存在无主面')
    return dict(points=points,faces=faces,cells=cells,kinds=kinds,boundaries=boundaries,label=label,
                internalFaces=int(np.sum(~is_boundary)))


def compare_source(grid,path):
    """模块二：逐原节点/单元核对；壁面按连接关系双射，不依赖空间最近邻。"""
    with Path(path).open(encoding='utf-8-sig') as f:
        def line():
            for s in f:
                s=s.split('%')[0].strip()
                if s:return s
            raise ValueError('原SU2截断')
        def head(key):
            k,v=line().split('=',1);need(k.strip()==key,'原SU2头错误');return v.strip()
        need(head('NDIME')=='3','原SU2不是三维')
        nc=int(head('NELEM'));need(nc==len(grid['cells']),'原单元数改变')
        for i in range(nc):
            row=list(map(int,line().split()));kind=row[0];need(kind in WIDTH,'源单元类型不支持')
            ids=row[1:1+WIDTH[kind]]
            need(len(row) in (WIDTH[kind]+1,WIDTH[kind]+2),'源单元宽度不符')
            if len(row)==WIDTH[kind]+2:need(row[-1]==i,'源单元ID不符')
            need(kind==grid['kinds'][i] and np.array_equal(ids,grid['cells'][i,:len(ids)]),'原单元顺序/节点解释改变')
        nv=int(head('NPOIN'));need(nv==len(grid['points']),'原节点数改变')
        error=0.
        for i in range(nv):
            row=line().split();xyz=np.array([float(v) for v in row[:3]])
            if len(row)==4:need(int(row[3])==i,'原节点ID不符')
            error=max(error,float(np.max(np.abs(xyz-grid['points'][i]))))
        need(error<=2e-14,'原节点坐标改变')
        nb=int(head('NMARK'));need(nb==len(grid['boundaries']),'边界集合改变');mapping={}
        for _ in range(nb):
            tag=head('MARKER_TAG');n=int(head('MARKER_ELEMS'))
            need(tag in grid['boundaries'],'原边界标签丢失')
            native=grid['boundaries'][tag];need(len(native['ids'])==n,'原边界面数改变')
            keys={tuple(sorted(v for v in grid['faces'][fid] if v>=0)):j for j,fid in enumerate(native['ids'])}
            need(len(keys)==n,'原生边界非双射')
            original_ids=[];rowmap=[];seen=set()
            for _ in range(n):
                row=list(map(int,line().split()));ids=row[1:];key=tuple(sorted(ids))
                need((row[0],len(ids)) in ((5,3),(9,4)),'原面类型错误')
                need(key in keys and key not in seen,'原面丢失或重复')
                seen.add(key);rowmap.append(keys[key]);original_ids.append(ids)
            mapping[tag]=dict(nativeRows=np.array(rowmap,dtype=int),originalNodes=original_ids)
    return mapping,error


def read_field(case,snapshot,ncells):
    """原生二进制为逐变量double数组，版本和列名由准确metadata提供。"""
    base=Path(case)/'lmrsim/snapshots'
    text=(base/'fluid.metadata').read_text()
    need('revision-id: f53f460' in text and 'version: "1.0"' in text,'场元数据版本不符')
    names=re.findall(r"^\s*- '([^']+)'\s*$",text,re.M)
    need(names and len(names)==len(set(names)),'场变量缺失或重复')
    values=np.fromfile(base/('%04d'%snapshot)/'fluid-0000',dtype='<f8')
    need(len(values)==len(names)*ncells,'场截断')
    fields=dict(zip(names,values.reshape(len(names),ncells)))
    for key,value in fields.items():
        # 尚未计算伪时间步的初场允许dt_local为NaN，物理场仍逐列严格有限。
        if key=='dt_local' and snapshot==0 and np.isnan(value).all():
            continue
        need(np.isfinite(value).all(),'场非有限值: '+key)
    return fields


def audit_run_scope(manifest,case=None):
    """未通过真实小例相容门槛的实数v1仅供有限诊断，禁止进入完整域或延长发展。"""
    c=manifest['condition']
    if 'explicit' in c:
        from eilmerExplicitProfile import validate_explicit,DEVELOPMENT_PROFILE,DENSE_PROFILE,DENSE_DEVELOPMENT_PROFILE
        validate_explicit(c['explicit'])
        if c['explicit']['profile']==DEVELOPMENT_PROFILE:
            # v184唯一原域在初始化耗尽内存；保留登记证据检查，但不再放行同一方案。
            raise ValueError('rah-explicit-local-development-v1原完整域初始化OOM，当前方案关闭；需独立内存改动及新证据')
        elif c['explicit']['profile']==DENSE_DEVELOPMENT_PROFILE:
            # v186确认缩减容量已生效，但原域仍在初始化OOM；同包候选关闭。
            raise ValueError('rah-explicit-dense-development-v1原域初始化仍OOM，当前候选关闭；需独立内存方案和证据')
        else:need(manifest['mesh']['cells']<=1000,'显式完整域尚未登记候选预算；当前仅最小接口')
        if c['explicit']['profile'] in (DENSE_PROFILE,DENSE_DEVELOPMENT_PROFILE):
            from eilmerDenseRuntime import audit_dense_runtime
            audit_dense_runtime(c)
    if c.get('development',{}).get('profile')=='rah-nk-real-v1':
        need(manifest['mesh']['cells']<=1000 and c['maxNewtonSteps']<=8,
             'rah-nk-real-v1小例比较未通过；完整域/延长计算关闭，需独立新方案证据')


def audit_numerics(nk,c):
    """按具名合同审核原生已生成的阶段与输出；旧小检查配置保持原门槛。"""
    need(nk['max_newton_steps']==c['maxNewtonSteps'] and nk['write_loads'],'步数或载荷输出不符')
    need(nk['use_physicality_check'] and nk['use_line_search'],'成熟状态检查未启用')
    if 'development' not in c:
        need(nk['min_cfl']==.1 and nk['max_cfl']==.1 and all(v[1]==.1 for v in nk['cfl_schedule']),
             '原生CFL没有按准备合同实际设置')
        return
    d=c['development']
    need(d['profile'] in ('rah-nk-development-v1','rah-nk-real-v1'),'发展方案错误')
    if d['profile']=='rah-nk-real-v1':
        # 明确实数Jv和预条件矩阵两处差分尺度，禁止保留复数默认1e-30。
        need(nk['use_scaling'] and nk['use_real_valued_frechet_derivative'],'实数缩放/导数合同错误')
        need(nk['frechet_derivative_perturbation']==2.0**-26 and
             nk['preconditioner_perturbation']==2.0**-26,'实数差分尺度错误')
        need(nk['use_preconditioner'] and nk['preconditioner']=='ilu' and nk['ilu_fill']==0,'实数预条件合同错误')
    need(nk['number_of_phases']==2 and nk['max_steps_in_initial_phases']==[d['initialOrderSteps']],'阶段切换错误')
    need(nk['min_cfl']==.01 and nk['max_cfl']==1000,'发展CFL合同错误')
    need(nk['steps_between_loads_update']==d['loadsEvery'] and nk['steps_between_snapshots']==d['snapshotsEvery'] and
         nk['total_snapshots']==d['snapshotCount'] and nk['write_snapshot_on_last_step'] and nk['write_loads_on_last_step'],'保存合同错误')
    for i,start in enumerate([.1,1.]):
        p=nk['NewtonKrylovPhase_%d'%i]
        need(p['use_auto_cfl'] and p['start_cfl']==start and p['residual_interpolation_order']==i+1 and
             p['jacobian_interpolation_order']==1 and p['linear_solve_tolerance']==.01 and
             p['threshold_relative_residual_for_cfl_growth']==.9 and p['auto_cfl_exponent']==.75,'原生阶段分派错误')


def audit_prepared(case):
    """核对真正落地的原生配置、M/方向/热边界、初场和原面全体身份。"""
    case=Path(case);m=json.loads((case/'case_manifest.json').read_text());c=m['condition']
    need(m['sourceCommit']=='f53f4609a0331d48efee69a4e4f3c3598378cc03','来源提交不符')
    need(m['mesh']['expectGmshOrderForWedges'] is False,'未采用已审核原生楔体顺序')
    need(sha(m['mesh']['nativeSourcePath'])==m['mesh']['sourceSHA256'],'原网格身份改变')
    for name,digest in m['filesSHA256'].items():
        need(sha(case/name)==digest,'准备输入变化: '+name)
    gridpath=case/'lmrsim/grid/grid-0000';grid=read_grid(gridpath)
    mapping,coorderr=compare_source(grid,m['mesh']['nativeSourcePath'])
    cfg=json.loads((case/'lmrsim/config').read_text())
    need(cfg['dimensions']==3 and cfg['viscous'] and cfg['turbulence_model']=='k_omega','原生物理模型不符')
    if 'explicit' in c:
        from eilmerExplicitProfile import audit_explicit_config
        audit_explicit_config(cfg,c,json.loads((case/'lmrsim/control').read_text()))
    else:
        audit_numerics(json.loads((case/'lmrsim/nk-config').read_text()),c)
    walltag=c['boundaryMarkers']['wall'];wallid=grid['boundaries'][walltag]['index']
    wall=cfg['block_0']['boundary_'+str(wallid)]
    actions=wall['pre_spatial_deriv_action_at_bndry_faces']
    need(wall['type']=='wall_no_slip_fixed_t' and any(a['type']=='zero_slip_wall_velocity' for a in actions),'热壁速度边界不符')
    tw=[a['Twall'] for a in actions if a['type']=='fixed_temperature']
    need(tw==[c['wallTemperatureK']],'给定壁温不符')
    native=json.loads((case/'native-input-state.json').read_text());speed=math.sqrt(sum(native[k]**2 for k in ('velx','vely','velz')))
    mach=speed/native['a']
    if c['flowMode']=='nativeMach':
        need(abs(mach-c['Mach'])<2e-13,'原生Mach不符')
        alpha,beta=map(math.radians,(c['alphaDeg'],c['betaDeg']))
        expected=np.array([-math.cos(alpha)*math.cos(beta),-math.sin(beta),math.sin(alpha)*math.cos(beta)])
        actual=np.array([native[k] for k in ('velx','vely','velz')])/speed
        need(np.max(np.abs(expected-actual))<2e-14,'来流方向不符')
        need(abs(native['tke']/(1.5*(c['turbulenceIntensity']*speed)**2)-1)<2e-13,'原生k不符')
        need(abs(native['rho']*native['tke']/native['mu']/native['omega']/c['eddyViscosityRatio']-1)<2e-13,'原生黏度比不符')
    field=read_field(case,0,len(grid['cells']))
    for variable,key in [('rho','rho'),('p','p'),('T','T'),('e','u'),('vel.x','velx'),('vel.y','vely'),('vel.z','velz'),('tq-tke','tke'),('tq-omega','omega')]:
        need(np.max(np.abs(field[variable]-native[key]))<=2e-13*max(1,abs(native[key])),'原生初场不符: '+variable)
    need(np.all(field['vol']>0),'原生单元体积非正')
    need(abs(field['vol'].sum()/m['mesh']['discreteVolumeM3']-1)<2e-10,'原生体积与独立离散体积不符')
    return dict(passed=True,nodes=len(grid['points']),cells=len(grid['cells']),internalFaces=grid['internalFaces'],
                typeCounts={str(k):int(np.sum(grid['kinds']==k)) for k in WIDTH},
                boundaryCounts={k:len(v['ids']) for k,v in grid['boundaries'].items()},
                maximumOriginalCoordinateDifference=coorderr,nativeDiscreteVolumeM3=float(field['vol'].sum()),
                nativeMinimumCellVolumeM3=float(field['vol'].min()),nativeMach=mach,
                initialLocalTimeStepDefined=bool(np.isfinite(field['dt_local']).all()),
                nativeFlowState=native,originalCellOrderAndBoundaryBijection=True,
                wallTemperatureGiven=True,wallTemperatureSolved=False,newQualifiedAreaM2=0,
                sourceMeshSHA256=sha(m['mesh']['nativeSourcePath']),nativeGridSHA256=sha(gridpath)),grid,mapping


def export_loads(case,loads,output,run_receipt):
    """模块三：入壁q=-outsign*(q_cond+q_diff)，全原面输出，不按有限值授予资格。"""
    case=Path(case).resolve();loads=Path(loads).resolve()
    need(run_receipt is not None,'缺少准确原生运行来源')
    receipt=json.loads(Path(run_receipt).read_text())
    need(receipt['stage']=='run' and receipt['passed'] and receipt['case']==str(case),'拒绝失败/错源运行')
    manifest=json.loads((case/'case_manifest.json').read_text())
    if manifest['condition'].get('explicit',{}).get('profile') in ('rah-explicit-dense-v1','rah-explicit-dense-development-v1'):
        from eilmerDenseRuntime import audit_dense_execution
        audit_dense_execution(manifest['condition'],receipt)
    # 旧回执可能只有启动器exit0；必须只读核对真正原生完成，失败历史不改写。
    completion=audit_completion(case,completion_log_for_receipt(case,run_receipt,receipt))
    need(receipt['caseManifestSHA256']==sha(case/'case_manifest.json'),'运行工况来源不符')
    for cmd in receipt['commands']:
        need(cmd['passed'] and cmd['inputsStable'] and cmd['nativeExitCode']==0 and cmd['stopReason'] is None,'原生异常运行')
        for path,digest in cmd['inputSHA256'].items():
            # 开发源码可继续演化；原任务准确源码另存回执，计算输入不得改变。
            p=Path(path)
            if p.is_relative_to(case):need(sha(p)==digest,'原运行输入改变')
    outputs=receipt.get('outputSHA256',{})
    need(str(loads) in outputs,'载荷未绑定本次原生输出')
    for path,digest in outputs.items():need(sha(path)==digest,'原运行输出改变: '+path)
    prepared,grid,mapping=audit_prepared(case)
    m=json.loads((case/'case_manifest.json').read_text())
    tag=m['condition']['boundaryMarkers']['wall'];boundary=grid['boundaries'][tag]
    need(loads.parent.parent==case/'lmrsim/loads' and
         loads.name=='blk-0000-bndry-%d-wall.dat'%boundary['index'],'载荷位置/边界错误')
    metadata=np.atleast_2d(np.loadtxt(case/'lmrsim/loads/loads-metadata',skiprows=1))
    row=metadata[metadata[:,0]==int(loads.parent.name)]
    limit=m['condition']['explicit']['maxSteps'] if 'explicit' in m['condition'] else m['condition']['maxNewtonSteps']
    need(len(row)==1 and row[0,2]>0 and row[0,2]<=limit,'载荷步骤来源不符')
    a=np.genfromtxt(loads,names=True,deletechars='',encoding='utf8');a=np.atleast_1d(a)
    required=['pos.x','pos.y','pos.z','n.x','n.y','n.z','area','outsign','p','T','q_total','q_cond','q_diff']
    need(all(k in a.dtype.names for k in required),'载荷列不完整')
    need(len(a)==len(boundary['ids']),'载荷面数不符')
    need(all(np.isfinite(a[k]).all() for k in a.dtype.names),'载荷含非有限值')
    positions=np.column_stack([a[k] for k in ('pos.x','pos.y','pos.z')]);normals=np.column_stack([a[k] for k in ('n.x','n.y','n.z')])
    centers=[];areas=[];geomnorm=[]
    for fid in boundary['ids']:
        ids=grid['faces'][fid];p=grid['points'][ids[ids>=0]];center=p.mean(axis=0)
        av=np.sum(np.cross(p-center,np.roll(p,-1,axis=0)-center),axis=0)/2
        area=float(np.linalg.norm(av));need(area>0,'原壁面面积非正')
        centers.append(center);areas.append(area);geomnorm.append(av/area)
    areas=np.array(areas)
    need(np.max(np.abs(positions-centers))<1e-10,'载荷位置不匹配原面')
    need(np.max(np.abs(a['area']/areas-1))<2e-10,'载荷面积不匹配原面')
    need(np.max(np.abs(normals-geomnorm))<2e-10,'原生面法向不匹配')
    need(np.array_equal(a['outsign'],boundary['signs']),'载荷外向符号不符')
    need(np.max(np.abs(a['T']-m['condition']['wallTemperatureK']))<1e-9,'实际面壁温不符')
    need(np.all(a['p']>0) and np.all(a['T']>0),'壁面热力状态非正')
    need(max(float(np.max(np.abs(a[k]))) for k in ('vel.x','vel.y','vel.z'))<1e-12,'无滑移壁速度不符')
    heat_scale=max(1.,float(np.max(np.abs(a['q_total']))))
    need(np.max(np.abs(a['q_total']-a['q_cond']-a['q_diff']))<2e-12*heat_scale,'热流分项不闭合')
    order=mapping[tag]['nativeRows'];heat=-a['outsign']*a['q_total']
    output=Path(output);need(not output.exists(),'不覆盖原有原面CSV')
    method='Eilmer5 thermally-perfect air / Wilcox k-omega'
    numerics={}
    if 'explicit' in m['condition']:
        profile=m['condition']['explicit']['profile']
        method += ' / '+profile+' TVD RK3 local pseudo-time (unconverged research)'
        numerics=dict(numericsProfile=profile,nativeNumberType='real',
                      timeMeaning='local-pseudo-time-not-flight-time',compatibilityStatus='interface-only')
    if m['condition'].get('development',{}).get('profile')=='rah-nk-real-v1':
        method += ' / rah-nk-real-v1 real finite-difference Jv+ILU (comparison failed)'
        numerics=dict(numericsProfile='rah-nk-real-v1',nativeNumberType='real',
                      compatibilityStatus='small-case-comparison-failed')
    with output.open('x',newline='',encoding='utf8') as f:
        writer=csv.writer(f);writer.writerow(['originalFaceId','nativeBoundaryRow','nativeFaceId','area_m2',
            'pressure_Pa','heatFluxIntoWall_W_m2','givenWallTemperature_K','wallFaceFluidTemperature_K',
            'native_q_total','outsign','method','status','physicsQualified'])
        for original,j in enumerate(order):
            writer.writerow([original+1,int(j),int(boundary['ids'][j]),float(a['area'][j]),float(a['p'][j]),
                float(heat[j]),m['condition']['wallTemperatureK'],float(a['T'][j]),float(a['q_total'][j]),
                int(a['outsign'][j]),method,'research',False])
    return dict(passed=True,wallFaces=len(a),areaM2=float(areas.sum()),heatPowerIntoWallW=float(np.dot(heat,a['area'])),
                negativeHeatFaces=int(np.sum(heat<0)),positiveHeatFaces=int(np.sum(heat>0)),
                zeroHeatFaces=int(np.sum(heat==0)),newQualifiedAreaM2=0,originalFaceOrder=True,
                sourceLoadsSHA256=sha(loads),csvSHA256=sha(output),
                sourceRunReceiptSHA256=sha(run_receipt),nativeStep=int(row[0,2]),
                outputBinding=receipt.get('outputBinding'),nativeLocalYPlusQualified=False,
                limitation='符号/身份/原生热分项接口审核；不是独立温度梯度重放、收敛或实验验证',**numerics)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',required=True)
    p.add_argument('--report',required=True);p.add_argument('--loads');p.add_argument('--csv');p.add_argument('--run-receipt')
    a=p.parse_args();dest=Path(a.report);need(not dest.exists(),'拒绝覆盖审核回执');started=time.monotonic()
    try:
        r=audit_prepared(a.case)[0] if a.loads is None else export_loads(a.case,a.loads,a.csv,a.run_receipt)
    except Exception as e:
        r=dict(passed=False,error=repr(e),newQualifiedAreaM2=0)
    r['elapsedSeconds']=time.monotonic()-started
    dest.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8');print(json.dumps(r))
    raise SystemExit(0 if r['passed'] else 1)
