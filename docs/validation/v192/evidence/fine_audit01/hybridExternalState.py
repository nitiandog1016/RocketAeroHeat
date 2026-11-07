# -*- coding: utf-8 -*-
"""逐原面Euler状态交换；数值身份审核与工程热资格分开，禁止把q补零。"""
from pathlib import Path
import hashlib
import json
import numpy as np


def fingerprint(value):
    """模块一：工况/几何/物性使用规范JSON绑定，拒绝NaN和隐式重排。"""
    raw=json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode('utf-8')
    return hashlib.sha256(raw).hexdigest()


def require(ok,message):
    if not ok: raise ValueError(message)


def exchange(frame, surface, condition, model, source_files, thermo_evaluator):
    """审核一个已映射帧；thermo_evaluator必须由调用方绑定同一准确物性。

    这里只接收显式逐原面映射，不做空间最近邻、不解释未收敛场为资格。
    生产Eilmer读场/边界采样尚须独立实现；合成检查不能转授原生入口。
    """
    require(frame.get('schema')=='RAH.HybridEulerFrame.v1','外缘帧版本不符')
    require(frame.get('solver')=='Euler' and frame.get('turbulenceModel')=='none','不是具名Euler帧')
    for key,expected in [('surfaceSHA256',fingerprint(surface)),('conditionSHA256',fingerprint(condition)),
                         ('modelSHA256',fingerprint(model))]:
        require(frame.get(key)==expected,'过期或错误身份: '+key)
    require(frame.get('sourceFilesSHA256')==source_files,'源帧文件集合或SHA不符')
    for name,value in source_files.items():
        require(Path(name).is_file() and hashlib.sha256(Path(name).read_bytes()).hexdigest()==value,'源帧已改变')
    require(bool(source_files),'缺少可回读的源帧')
    require(type(frame.get('snapshot')) is int and frame['snapshot']>=0,'缺少原生帧序号')
    vertices=np.asarray(surface['vertices'],dtype=float)
    faces=np.asarray(surface['facesOneBased'])
    require(vertices.ndim==2 and vertices.shape[1]==3 and np.isfinite(vertices).all(),'顶点错误')
    require(faces.ndim==2 and faces.shape[1]==3 and np.issubdtype(faces.dtype,np.integer),'原三角面编号错误')
    require(faces.size>0 and faces.min()>=1 and faces.max()<=len(vertices),'原三角面越界')
    xyz=vertices[faces-1];areas=np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)/2
    require(np.all(areas>0),'退化面')
    n=len(faces)
    require(surface['faceIDs']==list(range(1,n+1)),'原面ID须完整唯一')
    require(frame.get('faceIDs')==surface['faceIDs'],'外缘缺面/重排/重复')
    require(np.isclose(areas.sum(),frame.get('originalAreaM2',float('nan')),rtol=0,atol=1e-12),'原面面积分母变化')
    centers=np.asarray(frame.get('faceCentersM'),dtype=float)
    require(centers.shape==(n,3) and np.allclose(centers,xyz.mean(axis=1),rtol=0,atol=1e-12),'映射中心不符')
    fields=frame.get('fields',{})
    keys=['pressurePa','temperatureK','densityKgm3','internalEnergyJkg','enthalpyJkg','totalEnthalpyJkg','entropyJkgK']
    values={key:np.asarray(fields.get(key),dtype=float) for key in keys}
    require(all(a.shape==(n,) and np.isfinite(a).all() for a in values.values()),'外缘状态不完整/非有限')
    velocity=np.asarray(fields.get('velocityMps'),dtype=float)
    require(velocity.shape==(n,3) and np.isfinite(velocity).all(),'三维速度不完整')
    p,T,rho=[values[key] for key in keys[:3]]
    require(np.all(p>0)&np.all(T>0)&np.all(rho>0),'非正原始状态')
    tmin,tmax=model['temperatureRangeK']
    require(np.all((T>=tmin)&(T<=tmax)),'超出本接口已声明温域')
    # 模块二：同一实际EOS独立回算；恒等式不能替代热物性审核。
    expected=thermo_evaluator(p,T)
    errors={}
    for name in ('densityKgm3','internalEnergyJkg','enthalpyJkg','entropyJkgK'):
        target=np.asarray(expected[name],dtype=float)
        require(target.shape==(n,) and np.isfinite(target).all(),'物性回调返回错误')
        error=float(np.max(abs(values[name]-target)/np.maximum(1,abs(target))))
        errors[name]=error
        require(error<=2e-10,'同EOS状态不符: '+name)
    energy=values['internalEnergyJkg']+p/rho
    h0=values['enthalpyJkg']+np.sum(velocity**2,axis=1)/2
    require(np.all(abs(energy-values['enthalpyJkg'])<=2e-10*np.maximum(1,abs(energy))),'h与e/p/rho不相容')
    require(np.all(abs(h0-values['totalEnthalpyJkg'])<=2e-10*np.maximum(1,abs(h0))),'三维动能/总焓不相容')
    # 模块三：完整交换不等于附着外缘已可用；Euler没有黏性热流。
    return dict(schema='RAH.HybridExternalState.v1',faceIDs=surface['faceIDs'],originalAreaM2=float(areas.sum()),
                surfaceSHA256=frame['surfaceSHA256'],conditionSHA256=frame['conditionSHA256'],
                model=model,modelSHA256=frame['modelSHA256'],sourceFilesSHA256=source_files,
                snapshot=frame['snapshot'],fields=fields,algebraicInterfacePassed=True,thermodynamicErrors=errors,
                method='mapped Euler state, same-EOS checked; not yet an attached boundary-layer edge',
                attachedEdgeQualified=False,engineeringHeatReady=False,physicalQualification=False,
                wallHeatFluxWm2=[None]*n,wallTemperatureSolvedK=[None]*n)
