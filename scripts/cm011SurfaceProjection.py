"""将不同三角剖分的 CM-011 低速试件投影到由实际 STL 测得的求解壁面。"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from cm011Tunnel import D, REGIONS, meridian

FACET_TOL_M = 5e-6


def profile_positions(points, line, tree):
    """在加密的名义母线邻段中求最近点，得到顶点偏差和连续母线位置。"""
    xr=np.column_stack((points[:,0],np.hypot(points[:,1],points[:,2])))
    _,nearest=tree.query(xr)
    starts=np.column_stack((np.maximum(0,nearest-1),np.minimum(len(line)-2,nearest)))
    ab=line[starts+1]-line[starts]
    av=xr[:,None,:]-line[starts]
    fraction=np.clip(np.sum(av*ab,axis=2)/np.sum(ab*ab,axis=2),0,1)
    distance=np.linalg.norm(av-fraction[:,:,None]*ab,axis=2)
    take=np.argmin(distance,axis=1)
    i=np.arange(len(points))
    return distance[i,take],(starts[i,take]+fraction[i,take])*48/(len(line)-1)


def check_topology(vertices,faces):
    edges=np.vstack((faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]))
    pairs,index,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    winding=np.bincount(index,weights=np.sign(edges[:,1]-edges[:,0]))
    if (np.any(counts!=2) or np.any(winding!=0) or len(vertices)-len(pairs)+len(faces)!=2 or
            len(np.unique(vertices,axis=0))!=len(vertices)):
        raise ValueError("CM-011 表面有开口、非流形边、重复顶点或面方向不一致")
    graph=coo_matrix((np.ones(len(pairs)),pairs.T),shape=(len(vertices),len(vertices)))
    if connected_components(graph,directed=False,return_labels=False)!=1:
        raise ValueError("CM-011 表面须为单一连通体")


def measured_meridian(vertices,faces,nominal):
    """从内部发射射线，与原三角面求交，按方位平均得到半域体网所需的母线。"""
    center=.5*(nominal[0,0]+nominal[-1,0])
    local=vertices/D-np.array([center,0.,0.])
    radius=np.linalg.norm(local,axis=1)
    if np.any(radius==0): raise ValueError("导入表面经过试件内部中心")
    u=(local/radius[:,None])[faces]
    abc=local[faces]
    normal=np.cross(abc[:,1]-abc[:,0],abc[:,2]-abc[:,0])
    offsets=np.einsum("ij,ij->i",normal,abc[:,0])
    if np.any(offsets<=0): raise ValueError("CM-011 有退化、内向或径向折返的三角面")
    side=np.cross(u,np.roll(u,-1,axis=1))
    det=np.einsum("ij,ij->i",side[:,0],u[:,2])
    denom=1+np.sum(u*np.roll(u,-1,axis=1),axis=(1,2))
    angle=float(np.sum(2*np.arctan2(det,denom)))
    if not np.isclose(angle,4*np.pi,rtol=0,atol=1e-8):
        raise ValueError("CM-011 三角面没有完整单层径向覆盖")
    side/=np.linalg.norm(side,axis=2)[:,:,None]
    caps=u.sum(axis=1);caps/=np.linalg.norm(caps,axis=1)[:,None]
    cap_size=np.linalg.norm(u-caps[:,None,:],axis=2).max()
    if cap_size>=np.sqrt(2): raise ValueError("面元的径向投影跨度过大")
    tree=cKDTree(caps)
    direction=nominal-[center,0.]
    direction/=np.linalg.norm(direction,axis=1)[:,None]
    phi=np.arange(48)*2*np.pi/48
    dx,dr=direction.T
    rays=np.stack(np.broadcast_arrays(dx[:,None],dr[:,None]*np.cos(phi),dr[:,None]*np.sin(phi)),axis=-1)
    flat=rays.reshape(-1,3)
    nearby=tree.query_ball_point(flat,cap_size+1e-10)
    hits=np.empty(len(flat))
    for i,(ray,candidates) in enumerate(zip(flat,nearby)):
        ids=np.asarray(candidates,dtype=int)
        inside=np.all(np.einsum("ijk,k->ij",side[ids],ray)>=-2e-12,axis=1)
        ids=ids[inside]
        if not len(ids): raise ValueError("实际三角面未覆盖某一求解射线")
        t=offsets[ids]/np.einsum("ij,j->i",normal[ids],ray)
        if min(t)<=0 or np.ptp(t)*D>1e-11: raise ValueError("求解射线与多个不同深度的表面相交")
        hits[i]=np.mean(t)
    distances=hits.reshape(len(nominal),len(phi))
    spread=D*np.ptp(distances,axis=1).max()
    if spread>FACET_TOL_M: raise ValueError("实际 STL 的周向半径变化超出 CM-011 试件范围")
    wall=[center,0.]+direction*np.mean(distances,axis=1)[:,None]
    wall[[0,-1],1]=0
    info=dict(solidAngleSr=angle,rayCount=len(flat),maxAzimuthalSpreadM=float(spread),
              maxMeanMeridianShiftM=float(D*np.linalg.norm(wall-nominal,axis=1).max()))
    return wall,info


def project_vertices(vertices,coarse,center):
    """沿实际极向与求解母线线段相交，返回线性插值的连续站号和距离。"""
    rel=coarse-np.array([center,0.])
    xr=np.column_stack((vertices[:,0]/D-center,np.hypot(vertices[:,1],vertices[:,2])/D))
    angles=np.arctan2(rel[:,1],-rel[:,0])
    if np.any(np.diff(angles)<=0): raise ValueError("求解壁面站的径向次序无效")
    dirs=xr/np.linalg.norm(xr,axis=1)[:,None]
    left=np.clip(np.searchsorted(angles,np.arctan2(dirs[:,1],-dirs[:,0]))-1,0,len(rel)-2)
    delta=rel[left+1]-rel[left]
    t=-np.cross(rel[left],dirs)/np.cross(delta,dirs)
    if not np.isfinite(t).all() or np.any((t < -1e-12)|(t > 1+1e-12)):
        raise ValueError("导入顶点落在求解母线段之外")
    t=np.clip(t,0,1)
    gap=D*np.linalg.norm(rel[left]+t[:,None]*delta-xr,axis=1)
    return left+t,gap


def register_surface(vertices,faces,vertex_tol):
    idx=faces.astype(int)-1
    check_topology(vertices,idx)
    dense=meridian(256)[0]*D
    tree=cKDTree(dense)
    dist,station=profile_positions(vertices,dense,tree)
    if max(dist)>vertex_tol:
        raise ValueError("CM-011 节点超出名义轮廓，请核对尺度、朝向和局部几何")
    tri=vertices[idx]
    samples=np.vstack((tri.mean(axis=1),(.5*(tri+np.roll(tri,-1,axis=1))).reshape(-1,3)))
    facet_gap=profile_positions(samples,dense,tree)[0].max()
    if facet_gap>FACET_TOL_M:
        raise ValueError("CM-011 面心或边中点偏离外形过大，需要更细的 STL 剖分")
    nominal,_=meridian(2)
    wall,info=measured_meridian(vertices,idx,nominal)
    coordinate,gap=project_vertices(vertices,wall[::2],.5*(nominal[0,0]+nominal[-1,0]))
    _,bands=meridian(1)
    faces_s=station[idx]
    regions=bands[np.minimum(faces_s.mean(axis=1).astype(int),len(bands)-1)]
    boundaries=np.flatnonzero(np.diff(bands))+1
    crossing=(faces_s.min(axis=1)[:,None]<boundaries-1e-5)&(faces_s.max(axis=1)[:,None]>boundaries+1e-5)
    body=dict(V=vertices.tolist(),F=(idx+1).tolist(),wallStation=coordinate.tolist(),
              regionCode=regions.tolist(),regionNames=REGIONS,meridian=wall.tolist(),diameterM=D,grid="fine",extentD=12.)
    info.update(registrationMethod="CM-011 radial triangle intersections",vertexToleranceM=vertex_tol,
                maxVertexDistanceM=float(max(dist)),maxFacetSampleDistanceM=float(facet_gap),
                facetSampleToleranceM=FACET_TOL_M,maxVertexToSolverWallM=float(max(gap)),
                regionCrossingFaces=int(np.any(crossing,axis=1).sum()))
    return body,info
