"""Run14原生二维轴对称审核、壁面回读与公开测站对照；每弧度面积单独处理。"""
from pathlib import Path
import argparse
import csv
import gzip
import json
import math
import re
import time
import numpy as np
from runEilmerProcess import sha

R=8.31451/.0280134
CP=3.5*R


def need(ok,msg):
    if not ok:raise ValueError(msg)


def grid(case,block):
    p=Path(case)/('lmrsim/grid/grid-%04d.gz'%block)
    with gzip.open(p,'rt') as f:
        head=[f.readline().strip() for _ in range(6)]
        need(head[0]=='structured_grid 1.1' and head[2]=='dimensions: 2','网格类型不符')
        ni,nj,nk=[int(x.split(':')[1]) for x in head[3:]];need(nk==1,'不是二维域')
        a=np.loadtxt(f)
    need(a.shape==(ni*nj,3) and np.isfinite(a).all(),'网格截断/非有限')
    return a.reshape(nj,ni,3)


def field(case,block,snapshot):
    base=Path(case)/'lmrsim/snapshots';meta=(base/'fluid.metadata').read_text()
    need('revision-id: f53f460' in meta,'原生版本不符')
    names=re.findall(r"^\s*- '([^']+)'\s*$",meta,re.M)
    a=np.fromfile(base/('%04d'%snapshot)/('fluid-%04d'%block),dtype='<f8')
    need(len(a)==len(names)*160*96,'场大小不符')
    result=dict(zip(names,a.reshape(len(names),96,160)))
    for key,value in result.items():
        if key=='dt_local' and snapshot==0 and np.isnan(value).all():continue
        need(np.isfinite(value).all(),'场非有限: '+key)
    for key in ['p','rho','T','vol']:need(np.all(result[key]>0),'物理量非正: '+key)
    need(np.max(np.abs(result['p']/(result['rho']*R*result['T'])-1))<3e-13,'EOS不闭合')
    need(np.max(np.abs(result['e']/(2.5*R*result['T'])-1))<3e-13,'原生内能参考不符')
    return result


def audit(case,reference):
    start=time.monotonic();case=Path(case).resolve();m=json.loads((case/'case_manifest.json').read_text())
    need(m['schema']=='RAH.EilmerRun14Control.v1','工况身份错误')
    for name,digest in m['filesSHA256'].items():need(sha(case/name)==digest,'输入变化')
    cfg=json.loads((case/'lmrsim/config').read_text());nk=json.loads((case/'lmrsim/nk-config').read_text())
    need(cfg['axisymmetric'] and cfg['dimensions']==2 and cfg['turbulence_model']=='none' and cfg['viscous'],'模型不符')
    need(nk['max_newton_steps']==2000 and nk['max_steps_in_initial_phases']==[50],'步数合同不符')
    for i,startcfl in enumerate([.1,1.]):
        phase=nk['NewtonKrylovPhase_%d'%i]
        need(phase['use_auto_cfl'] and phase['start_cfl']==startcfl and phase['residual_interpolation_order']==i+1,'原生阶段合同不符')
    grids=[grid(case,b) for b in (0,1)];minjac=math.inf;minheight=math.inf;maxheight=0.;area=0.;volume=0.;nativevol=0.
    need(np.max(np.abs(grids[0][:,-1]-grids[1][:,0]))<1e-14,'两块接口不匹配')
    for b,g in enumerate(grids):
        need(g.shape==(97,161,3),'原生网格维数不符')
        p=np.stack([g[:-1,:-1,:2],g[:-1,1:,:2],g[1:,1:,:2],g[1:,:-1,:2]],axis=2)
        # 双线性四边形Jacobian在四参考角的正性；不以面积正代替无折叠。
        for k in range(4):
            e=p[:,:,(k+1)%4]-p[:,:,k];v=p[:,:,(k-1)%4]-p[:,:,k]
            j=e[:,:,0]*v[:,:,1]-e[:,:,1]*v[:,:,0]
            need(np.all(j>0),'四边形参考角翻转');minjac=min(minjac,float(j.min()))
        nxt=np.roll(p,-1,axis=2);cross=p[:,:,:,0]*nxt[:,:,:,1]-nxt[:,:,:,0]*p[:,:,:,1]
        moment=np.sum((p[:,:,:,1]+nxt[:,:,:,1])*cross)/6
        volume+=2*math.pi*float(moment)
        wall=g[0,:,:2];slope=0 if b==0 else math.tan(math.pi/6)
        y=.0325+slope*(wall[:,0]-.1017 if b else wall[:,0])
        need(np.max(np.abs(wall[:,1]-y))<2e-14,'原物理壁形状改变')
        edge=np.diff(wall,axis=0);ds=np.linalg.norm(edge,axis=1);area+=float(np.sum(2*math.pi*wall[:-1,1]*ds+math.pi*np.diff(wall[:,1])*ds))
        normal=np.column_stack([-edge[:,1],edge[:,0]])/ds[:,None]
        centers=p.mean(axis=2);wc=(wall[:-1]+wall[1:])/2
        distances=np.einsum('ij,ij->i',centers[0]-wc,normal)
        need(np.all(distances>0),'首层中心不在流体侧');minheight=min(minheight,float(distances.min()));maxheight=max(maxheight,float(distances.max()))
        f=field(case,b,0);nativevol+=2*math.pi*float(f['vol'].sum())
        expectedp=31.88 if b==0 else 31.88/5;expectedu=2304. if b==0 else 0.
        need(np.max(np.abs(f['p']-expectedp))<1e-13 and np.max(np.abs(f['vel.x']-expectedu))<1e-13,'初态与合同不符')
        bc=cfg['block_%d'%b]['boundary_south']
        actions=bc['pre_spatial_deriv_action_at_bndry_faces']
        need(bc['type']=='wall_no_slip_fixed_t' and [a['Twall'] for a in actions if a['type']=='fixed_temperature']==[295.2],'热壁错误')
    expectedarea=2*math.pi*.0325*.1017+math.pi*(.0325+.0325+math.tan(math.pi/6)*(.22-.1017))*(.22-.1017)/math.cos(math.pi/6)
    need(abs(area/expectedarea-1)<2e-13,'全周壁面积错误')
    gas=np.genfromtxt(case/'native-gas.csv',delimiter=',',names=True);comparisons=[]
    ref=json.loads(Path(reference).read_text())
    for a,known in zip(gas,ref['rows']):
        T=float(a['T']);need(T==known['T'],'温度锚点不符')
        mu=1.663e-5*(T/273)**1.5*380/(T+107);k=.0242*(T/273)**1.5*423/(T+150)
        need(abs(a['mu']/mu-1)<3e-15 and abs(a['k']/k-1)<3e-15,'原生Sutherland没有实际生效')
        need(abs(a['cp']/CP-1)<3e-15 and abs(a['e']/(2.5*R*T)-1)<3e-15,'原生定cp/能量错误')
        comparisons.append(dict(T=T,muRelativeToIndependent=float(a['mu']/known['referenceMu']-1),kRelativeToIndependent=float(a['k']/known['referenceK']-1)))
    need(all(abs(a['muRelativeToIndependent'])<.05 and abs(a['kRelativeToIndependent'])<.05 for a in comparisons if a['T']<=200),'低温输运5%筛选未过')
    return dict(passed=True,caseManifestSHA256=sha(case/'case_manifest.json'),cells=30720,wallFaces=320,
        fullCircumferenceWallAreaM2=area,independentAxisymmetricVolumeM3=volume,nativeAxisymmetricVolumeM3=nativevol,
        nativeVsExactMomentVolumeRelativeDifference=nativevol/volume-1,minimumCornerJacobian=minjac,
        firstCellNormalDistanceMinM=minheight,firstCellNormalDistanceMaxM=maxheight,
        nativeMach=2304/float(gas['a'][1]),nativeDensity=float(gas['rho'][1]),
        transportComparisons=comparisons,modelBiasAndMeshSensitivityPassed=False,newQualifiedAreaM2=0,
        elapsedSeconds=time.monotonic()-start)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--reference',required=True);p.add_argument('--report',required=True)
    a=p.parse_args();out=Path(a.report);need(not out.exists(),'回执已存在')
    try:r=audit(a.case,a.reference)
    except Exception as e:r=dict(passed=False,error=repr(e))
    out.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8');print(json.dumps(r));raise SystemExit(0 if r['passed'] else 1)
