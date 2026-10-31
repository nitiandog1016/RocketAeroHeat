"""从本程序已收敛平板解生成细网格守恒变量初场；不输入外部CFD，不推定物理资格。"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def config(p):
    pairs={}
    for line in p.read_text(encoding='utf-8-sig').splitlines():
        line=line.split('%')[0].strip()
        if not line or '=' not in line:continue
        key,value=line.split('=',1);assert key.strip() not in pairs
        value=value.strip()
        # 配置里的0.2与0.20000000000000001为同一double；不把打印格式当物性改变。
        if re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?',value):value=float(value)
        pairs[key.strip()]=value
    return pairs


def mesh_points(p):
    lines=p.read_text(encoding='ascii').splitlines();assert lines[0].strip()=='NDIME= 2'
    i=next(i for i,s in enumerate(lines) if s.startswith('NPOIN='));n=int(lines[i].split('=')[1])
    values=np.array([[float(x) for x in line.split()] for line in lines[i+1:i+1+n]])
    assert np.array_equal(values[:,2],np.arange(n));return values[:,:2]


def tensor_grid(points,values):
    x=np.unique(points[:,0]);y=np.unique(points[:,1]);assert len(x)*len(y)==len(points)
    ix=np.searchsorted(x,points[:,0]);iy=np.searchsorted(y,points[:,1])
    assert len(np.unique(ix*len(y)+iy))==len(points)
    g=np.empty((len(x),len(y),values.shape[1]));g[ix,iy]=values
    return x,y,g


def interpolate(x,y,g,points):
    tolerance=64*np.finfo(float).eps*max(1,np.max(np.abs(points)))
    lo=np.array([x[0],y[0]]);hi=np.array([x[-1],y[-1]])
    assert np.all(points>=lo-tolerance) and np.all(points<=hi+tolerance),'不能外推初始场'
    query=np.minimum(np.maximum(points,lo),hi)
    ix=np.clip(np.searchsorted(x,query[:,0],side='right')-1,0,len(x)-2)
    iy=np.clip(np.searchsorted(y,query[:,1],side='right')-1,0,len(y)-2)
    tx=(query[:,0]-x[ix])/(x[ix+1]-x[ix]);ty=(query[:,1]-y[iy])/(y[iy+1]-y[iy])
    weights=np.column_stack(((1-tx)*(1-ty),tx*(1-ty),(1-tx)*ty,tx*ty))
    assert np.all(weights>=0) and np.max(np.abs(weights.sum(axis=1)-1))<1e-14
    corners=np.stack((g[ix,iy],g[ix+1,iy],g[ix,iy+1],g[ix+1,iy+1]),axis=1)
    return np.einsum('ni,nij->nj',weights,corners),int(np.count_nonzero(query!=points))


def integral(x,y,g):
    # 节点中位双网格积分；初始化仅报告插值前后误差，不伪称严格保守映射。
    weights=lambda t: np.r_[np.diff(t)[0]/2,(t[2:]-t[:-2])/2,np.diff(t)[-1]/2]
    return np.einsum('i,j,ijk->k',weights(x),weights(y),g)


def prepare(source,target):
    source=Path(source).resolve();target=Path(target).resolve();target.relative_to(ROOT/'outputs')
    assert not (target/'engine_request.json').exists() and not (target/'initial_flow.csv').exists()
    origin=json.loads((source/'engine_process.json').read_text(encoding='utf8'))
    assert origin['status']=='terminal' and origin['exitCode']==0 and not origin['timedOut']
    assert all(origin['inputsStable'].values())
    for name in ('restart_flow.csv','history.csv','surface_flow.vtk'):
        assert sha(source/name)==origin['outputSHA256'][name],name
    with (source/'history.csv').open() as stream:rows=list(csv.reader(stream,skipinitialspace=True))
    last={key.strip():float(value) for key,value in zip(rows[0],rows[-1])}
    assert last['rms[Rho]']<-8 and last['rms[RhoE]']<-8
    assert 'All convergence criteria satisfied.' in (source/'engine_stdout.log').read_text(errors='replace')
    a=config(source/'case.cfg');b=config(target/'case.cfg')
    ignored={'ITER','OUTPUT_WRT_FREQ','WRT_RESTART_OVERWRITE','WRT_SURFACE_OVERWRITE'}
    assert {k:v for k,v in a.items() if k not in ignored}=={k:v for k,v in b.items() if k not in ignored},'物理/离散配置须一致'
    source_manifest=json.loads((source/'case_manifest.json').read_text(encoding='utf8'))
    target_manifest=json.loads((target/'case_manifest.json').read_text(encoding='utf8'))
    for key in ('length','height','upstreamLength','Mach','temperature','pressure','wallTemperature'):
        assert source_manifest['options'][key]==target_manifest['options'][key]
    assert target_manifest['options']['outerExtensionHeight']==0
    header=(source/'restart_flow.csv').open().readline().strip()
    assert header=='"PointID","x","y","Density","Momentum_x","Momentum_y","Energy"'
    raw=np.genfromtxt(str(source/'restart_flow.csv'),delimiter=',',skip_header=1)
    points=mesh_points(source/'mesh.su2');new_points=mesh_points(target/'mesh.su2')
    assert np.array_equal(raw[:,0],np.arange(len(points))) and raw.shape==(len(points),7)
    assert np.max(np.abs(points-raw[:,1:3]))<1e-14
    x,y,g=tensor_grid(points,raw[:,3:]);mapped,snaps=interpolate(x,y,g,new_points)
    identical,_=interpolate(x,y,g,points);assert np.array_equal(identical,raw[:,3:])
    rho=mapped[:,0];internal=mapped[:,3]-(mapped[:,1]**2+mapped[:,2]**2)/(2*rho)
    assert np.all(np.isfinite(mapped)) and np.all(rho>0) and np.all(internal>0),'初始场不能裁温度/密度补合法'
    T=internal/rho*(1.4-1)/287.05287
    wall=(new_points[:,1]==0)&(new_points[:,0]<=0)
    assert np.max(np.abs(mapped[wall,1:3]))<1e-10,'壁面必须保持无滑移'
    xn,yn,gn=tensor_grid(new_points,mapped);before=integral(x,y,g);after=integral(xn,yn,gn)
    reference=np.array([abs(before[0]),np.linalg.norm(before[1:3]),np.linalg.norm(before[1:3]),abs(before[3])])
    relative=np.abs(after-before)/reference;assert np.max(relative)<.005,'初场积分差超过运行前0.5%筛选'
    data=np.column_stack((np.arange(len(new_points)),new_points,mapped))
    np.savetxt(str(target/'initial_flow.csv'),data,delimiter=',',fmt=['%d']+['%.17e']*6,header=header,comments='')
    old=(target/'case.cfg').read_text(encoding='utf8')
    assert old.count('RESTART_SOL= NO')==1
    (target/'case_before_warm_start.cfg').write_bytes(old.encode('utf8'))
    text=old.replace('RESTART_SOL= NO','RESTART_SOL= YES')+'SOLUTION_FILENAME= initial_flow\nREAD_BINARY_RESTART= NO\n'
    (target/'case.cfg').write_bytes(text.encode('utf8'))
    inputs=[source/name for name in ('case.cfg','case_manifest.json','engine_process.json','restart_flow.csv','history.csv','mesh.su2')]
    inputs += [target/name for name in ('mesh.su2','case_before_warm_start.cfg','initial_flow.csv','case.cfg')]
    inputs.append(Path(__file__).resolve())
    report=dict(method='bilinear convex interpolation of dimensional conservative variables on rectilinear nodes',
                source=str(source),target=str(target),sourceNodes=len(points),targetNodes=len(new_points),
                sourceResiduals=last,integralBefore=before.tolist(),integralAfter=after.tolist(),
                integralRelativeChanges=relative.tolist(),roundoffOnlyCoordinateSnaps=snaps,
                minDensity=float(rho.min()),minTemperature=float(T.min()),maxTemperature=float(T.max()),
                wallNoSlipPreserved=True,wallNodeTemperatureIsNotPrescribedTw=True,
                pointwiseSourceIdentityPreserved=True,strictlyConservativeRemap=False,
                physicsQualified=False,inputSHA256={str(p):sha(p) for p in inputs})
    (target/'warm_start_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('source');parser.add_argument('target')
    args=parser.parse_args();report=prepare(args.source,args.target)
    print(json.dumps({k:report[k] for k in ('sourceNodes','targetNodes','integralRelativeChanges','minTemperature','maxTemperature')}))
