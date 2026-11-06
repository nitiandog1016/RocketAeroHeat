"""模块用途：单wall空间壁温的原生弱热流及有限段守恒分解，分开投影面积与几何面积。"""
import csv
from pathlib import Path
import numpy as np

from nacaWallFlux import geometry as marker_geometry

METHOD='single-wall-spatial-Tw-native-normal-distance-v1'
SEGMENT_METHOD='nodal-flux-vector-projected-median-subtriangle-clipping-v1'


def read_temperature_table(path):
    """模块一：沿用原生严格SI表，拒绝外插、乱序或多余列。"""
    with Path(path).open(encoding='utf8') as stream:
        assert stream.readline().strip()=='RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K'
        rows=list(csv.reader(stream))
    assert 2<=len(rows)<=10000 and all(len(row)==2 for row in rows)
    values=np.asarray(rows,float)
    assert np.isfinite(values).all() and np.all(np.diff(values[:,0])>0) and np.all(values[:,1]>0)
    return values


def temperature_at(x,table):
    """原生表在X上插值，绝不把s插值或最近测站静默替换为实际边界。"""
    x=np.asarray(x,float);table=np.asarray(table,float)
    assert np.isfinite(x).all() and np.all(x>=table[0,0]) and np.all(x<=table[-1,0])
    return np.interp(x,table[:,0],table[:,1])


def geometry(points,cells,markers,node_s=None):
    """模块二：重构节点原生法向/邻点，并保存真实中位子面和线性s坐标。"""
    assert 'wall' in markers and not any(k.startswith('wall_') for k in markers)
    points=np.asarray(points,float)
    g=marker_geometry(points,cells,{'wall_spatial':markers['wall']})['wall_spatial']
    ids=g['ids'];inverse=np.full(len(points),-1,int);inverse[ids]=np.arange(len(ids))
    delta=points[g['normalNeighbor']]-points[ids]
    euclidean=np.linalg.norm(delta,axis=1);projected=g['distance'].copy()
    # SU2 GeometryToolbox::NormalDistance的真实路径，不删5%保护项。
    g.update(projectedDistance=projected,euclideanDistance=euclidean,
             distance=np.maximum(projected,.05*euclidean),
             distanceFloorTriggered=projected<.05*euclidean)
    vertices=[];owners=[];scalars=[]
    if node_s is not None:
        node_s=np.asarray(node_s,float)
        assert node_s.shape==(len(points),) and np.isfinite(node_s[ids]).all()
    for row in markers['wall']:
        face=np.asarray(row[1:],int);q=points[face];center=q.mean(axis=0)
        sf=node_s[face] if node_s is not None else None
        for k,node in enumerate(face):
            following=(k+1)%len(face);previous=(k-1)%len(face)
            a=q[k];b=(a+q[following])/2;c=(a+q[previous])/2
            # 两个三角片均朝向流体，回加等于该节点的原中位对偶法向。
            vertices.extend(([a,center,b],[a,c,center]));owners.extend([inverse[node]]*2)
            if sf is not None:
                scalars.extend(([sf[k],sf.mean(),(sf[k]+sf[following])/2],
                                [sf[k],(sf[k]+sf[previous])/2,sf.mean()]))
    vertices=np.asarray(vertices);owners=np.asarray(owners,int)
    vectors=np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0])/2
    geometric=np.linalg.norm(vectors,axis=1)
    projected_area=np.einsum('ij,ij->i',vectors,g['normal'][owners])
    assert geometric.min()>0 and np.all(projected_area>0)
    rebuilt=np.zeros((len(ids),3));np.add.at(rebuilt,owners,vectors)
    target=g['normal']*g['area'][:,None]
    relative=np.max(np.linalg.norm(rebuilt-target,axis=1)/g['area'])
    assert relative<2e-12
    native_areas=np.bincount(owners,weights=projected_area,minlength=len(ids))
    assert np.max(abs(native_areas/g['area']-1))<2e-12
    g.update(patchVertices=vertices,patchOwnerLocal=owners,patchVectors=vectors,
             patchGeometricArea=geometric,patchProjectedArea=projected_area,
             patchS=np.asarray(scalars) if node_s is not None else None,
             discreteGeometricArea=float(geometric.sum()),nativeProjectedArea=float(g['area'].sum()),
             normalReconstructionRelativeError=float(relative),method=METHOD)
    return g


def evaluate(g,temperature,conductivity,given_tw):
    """模块三：正入壁导热通量，givenTw和实际流体壁温分别保留。"""
    t=np.asarray(temperature,float);k=np.asarray(conductivity,float);tw=np.asarray(given_tw,float)
    assert t.ndim==k.ndim==1 and t.shape==k.shape and tw.shape==g['ids'].shape
    assert np.isfinite(t).all() and np.isfinite(k).all() and np.isfinite(tw).all()
    assert np.all(k>0) and np.all(tw>0)
    q=k[g['ids']]*(t[g['normalNeighbor']]-tw)/g['distance']
    power=q*g['area']
    patch_power=q[g['patchOwnerLocal']]*g['patchProjectedArea']
    discrepancy=abs(float(power.sum()-patch_power.sum()))/max(1.,float(np.sum(abs(power))))
    assert discrepancy<2e-12
    return dict(ids=g['ids'],givenTw=tw,wallFluidTemperature=t[g['ids']],qIntoWall=q,
                nodalPowerIntoWall=power,powerIntoWall=float(power.sum()),patchPowerIntoWall=patch_power,
                patchPowerPartitionScaledError=float(discrepancy),method=METHOD,
                temperatureIsGiven=True,solidWallTemperatureSolved=False,physicalQualified=False)


def _clip(polygon,bound,keep_above):
    """在子三角片的线性s参数上裁切；保留真实三维交点，不投影回光滑球面。"""
    result=[]
    for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
        ia=a[3]>=bound if keep_above else a[3]<=bound
        ib=b[3]>=bound if keep_above else b[3]<=bound
        if ia:result.append(a)
        if ia!=ib:
            fraction=(bound-a[3])/(b[3]-a[3])
            x=a+fraction*(b-a);x[3]=bound;result.append(x)
    return np.asarray(result,float).reshape((-1,4))


def segment_integrals(g,flux,bounds,sector_radians=np.pi/2):
    """模块四：裁切实际中位子面，对原生节点热流向量作守恒积分。

    原生只给节点q_i和|ΣN_i|，没有输出子面q。采用q_i*n_i·N片为明确的
    几何分解；全域严格回加原生功率。几何面积Σ|N片|另列，不冒充投影面积。
    """
    assert g['patchS'] is not None
    bounds=np.asarray(bounds,float)
    assert bounds.ndim==2 and bounds.shape[1]==2 and np.isfinite(bounds).all()
    assert np.all(bounds[:,1]>bounds[:,0]) and 0<sector_radians<=2*np.pi
    assert bounds.min()>=g['patchS'].min()-1e-14 and bounds.max()<=g['patchS'].max()+1e-14
    vertices=np.concatenate((g['patchVertices'],g['patchS'][:,:,None]),axis=2)
    pmin=g['patchS'].min(axis=1);pmax=g['patchS'].max(axis=1);out=[]
    q=flux['qIntoWall'];scale=2*np.pi/sector_radians
    for lower,upper in bounds:
        inside=(pmin>=lower)&(pmax<=upper)
        candidates=np.flatnonzero((pmax>lower)&(pmin<upper)&~inside)
        pa=float(g['patchProjectedArea'][inside].sum());ga=float(g['patchGeometricArea'][inside].sum())
        power=float(flux['patchPowerIntoWall'][inside].sum());partial=0
        for index in candidates:
            poly=_clip(vertices[index],lower,True)
            if len(poly)<3:continue
            poly=_clip(poly,upper,False)
            if len(poly)<3:continue
            xyz=poly[:,:3]-poly[0,:3]
            normal=np.sum(np.cross(xyz,np.roll(xyz,-1,axis=0)),axis=0)/2
            geometric=float(np.linalg.norm(normal));owner=g['patchOwnerLocal'][index]
            projected=float(np.dot(normal,g['normal'][owner]))
            assert projected>=-1e-12*geometric, '反向局部投影不属于此分解范围'
            pa+=projected;ga+=geometric;power+=float(q[owner]*projected);partial+=1
        assert ga>0 and pa>0
        out.append(dict(boundsM=[float(lower),float(upper)],powerIntoWallSectorW=power,
            symmetryScaledPowerIntoWallFullW=power*scale,discreteGeometricAreaSectorM2=ga,
            nativeProjectedAreaSectorM2=pa,meanQByDiscreteGeometricAreaWm2=power/ga,
            meanQByNativeProjectedAreaWm2=power/pa,fullyIncludedSubtriangles=int(inside.sum()),
            intersectedSubtriangles=partial,sectorRadians=float(sector_radians),fullRevolutionScale=float(scale),
            method=SEGMENT_METHOD,symmetryScalingIsExperimentalEvidence=False,physicalQualified=False))
    return out


def verify_applied(path,points,g,expected):
    """模块五：核对原生实际写入表的全局编号/坐标/Tw，不能只核对输入文件。"""
    with Path(path).open(encoding='utf8') as stream:rows=list(csv.DictReader(stream))
    ids=np.array([int(row['globalNodeID']) for row in rows])
    order=np.argsort(ids);assert len(ids)==len(g['ids']) and np.array_equal(ids[order],g['ids'])
    xyz=np.array([[float(row[k]) for k in ('X_m','Y_m','Z_m')] for row in rows])[order]
    values=np.array([float(row['givenTemperature_K']) for row in rows])[order]
    coord=float(np.max(abs(xyz-points[g['ids']])))
    difference=float(np.max(abs(values-expected)))
    assert coord<2e-12 and difference<2e-12
    return dict(nodes=len(ids),maximumCoordinateDifferenceM=coord,maximumGivenTwDifferenceK=difference)
