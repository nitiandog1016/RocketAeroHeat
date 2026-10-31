"""SST固定版本的状态和近壁几何接口；不推定物理资格，不修改历史输入。"""
from pathlib import Path
import csv
import numpy as np
from scipy.spatial import cKDTree


PHASE_OUTPUT='native-noncompact-flow-primitive-before-turbulence-v1'


def phase_columns(nd):
    """固定官方输出注册顺序，二维/三维分开；未知扩展列仍拒绝。"""
    axes=['x','y','z'][:nd]
    return ['Pressure','Temperature','Mach','Pressure_Coefficient']+['Velocity_'+s for s in axes]+[
        'Laminar_Viscosity']+['Skin_Friction_Coefficient_'+s for s in axes]+['Heat_Flux','Y_Plus','Eddy_Viscosity']


def state(path, turbulent, phase_output=False):
    """按具名列绑定ρE、k/ω；只有显式合同才读取同次保存的原生primitive。"""
    with Path(path).open(encoding='utf-8-sig') as f:
        keys=next(csv.reader(f))
    nd=3 if 'z' in keys else 2
    expected=['PointID','x','y']+(['z'] if nd==3 else [])+['Density','Momentum_x','Momentum_y']
    if nd==3:expected+=['Momentum_z']
    expected+=['Energy']+(['Turb_Kin_Energy','Omega'] if turbulent else [])
    if phase_output:
        assert turbulent
        expected+=phase_columns(nd)
    assert keys==expected, '重启动列名/模型不匹配：'+str(keys)
    a=np.loadtxt(str(path),delimiter=',',skiprows=1)
    assert a.ndim==2 and a.shape[1]==len(keys) and np.isfinite(a).all()
    assert np.array_equal(a[:,0],np.arange(len(a)))
    rho=a[:,keys.index('Density')];assert min(rho)>0
    u=a[:,nd+2:2*nd+2]/rho[:,None]
    k=a[:,keys.index('Turb_Kin_Energy')] if turbulent else np.zeros(len(a))
    e=a[:,keys.index('Energy')]/rho-.5*np.sum(u*u,axis=1)-k
    assert min(e)>0 and min(k)>=0
    if turbulent:assert min(a[:,keys.index('Omega')])>0
    return a,keys,nd,rho,u,e


def mesh(path):
    """只读SU2原网格，保留编号；棱柱/四面体/六面体和二维四边形支持。"""
    with Path(path).open() as f:
        assert f.readline().startswith('NDIME=');f.seek(0)
        nd=int(f.readline().split('=')[1]);nc=int(f.readline().split('=')[1])
        cells=[list(map(int,f.readline().split())) for _ in range(nc)]
        n=int(f.readline().split('=')[1]);p=np.array([list(map(float,f.readline().split())) for _ in range(n)])
        assert np.array_equal(p[:,-1],np.arange(n));p=p[:,:nd]
        nm=int(f.readline().split('=')[1]);wall=[]
        for _ in range(nm):
            tag=f.readline().split('=')[1].strip();count=int(f.readline().split('=')[1])
            faces=[list(map(int,f.readline().split())) for _ in range(count)]
            if tag=='wall':wall=faces
    assert wall and nd in (2,3)
    return p,cells,wall


def triangle_distances(p,a,b,c):
    """点到三角形的精确距离：面内垂足或三条闭线段，不以最近顶点代替。"""
    ab=b-a;ac=c-a;ap=p-a
    aa=np.sum(ab*ab,axis=1);bb=np.sum(ab*ac,axis=1);cc=np.sum(ac*ac,axis=1)
    u=np.sum(ap*ab,axis=1);v=np.sum(ap*ac,axis=1);den=aa*cc-bb*bb
    assert np.all(den>0)
    s=(u*cc-v*bb)/den;t=(v*aa-u*bb)/den
    normal=np.cross(ab,ac);plane=abs(np.sum(ap*normal,axis=1))/np.linalg.norm(normal,axis=1)
    best=np.full(len(a),np.inf);inside=(s>=0)&(t>=0)&(s+t<=1);best[inside]=plane[inside]
    for lo,hi in ((a,b),(b,c),(c,a)):
        edge=hi-lo;w=np.clip(np.sum((p-lo)*edge,axis=1)/np.sum(edge*edge,axis=1),0,1)
        best=np.minimum(best,np.linalg.norm(p-lo-w[:,None]*edge,axis=1))
    return best


def wall_geometry(path):
    """分开核对正常邻点和SST最近邻壁距；搜索半径用三角包围球作完整性界。"""
    p,cells,wall=mesh(path);nd=p.shape[1];wallids=np.unique([v for row in wall for v in row[1:]])
    wallset=set(wallids);neighbors={int(i):set() for i in wallids}
    edges={5:((0,1),(1,2),(2,0)),9:((0,1),(1,2),(2,3),(3,0)),
           10:((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)),
           13:((0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)),
           12:((0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7))}
    for row in cells:
        typ=row[0];v=row[1:];assert typ in edges
        for a,b in edges[typ]:
            i,j=v[a],v[b]
            if i in wallset:neighbors[i].add(j)
            if j in wallset:neighbors[j].add(i)
    normal=np.zeros_like(p)
    for row in wall:
        v=np.array(row[1:])
        if nd==2:
            assert row[0]==3
            tangent=p[v[1]]-p[v[0]];n=np.array([-tangent[1],tangent[0]])/2
        else:
            assert row[0]==5, '首次SST三维适配要求原三角壁面'
            n=np.cross(p[v[1]]-p[v[0]],p[v[2]]-p[v[0]])/6
        normal[v]+=n
    unit=normal[wallids]/np.linalg.norm(normal[wallids],axis=1)[:,None]
    near=[];off=np.unique([j for i in wallids for j in neighbors[i] if j not in wallset])
    if nd==3:
        tri=np.array([row[1:] for row in wall]);a,b,c=(p[tri[:,j]] for j in range(3))
        center=(a+b+c)/3;radius=max(np.linalg.norm(p[tri]-center[:,None],axis=2).ravel())
        tree=cKDTree(center);vertices=cKDTree(p[wallids]);upper=vertices.query(p[off])[0]
        distances={}
        for j,bound in zip(off,upper):
            ix=tree.query_ball_point(p[j],radius+bound+1e-12)
            distances[j]=float(min(triangle_distances(p[j],a[ix],b[ix],c[ix])))
    else:
        seg=np.array([r[1:] for r in wall]);a=p[seg[:,0]];b=p[seg[:,1]];v=b-a
        distances={}
        for j in off:
            t=np.clip(np.sum((p[j]-a)*v,axis=1)/np.sum(v*v,axis=1),0,1)
            distances[j]=float(min(np.linalg.norm(p[j]-a-t[:,None]*v,axis=1)))
    nearest=[]
    for ii,i in enumerate(wallids):
        js=np.array(sorted(neighbors[i]));delta=p[js]-p[i]
        near.append(js[np.argmax(delta@unit[ii]/np.linalg.norm(delta,axis=1))])
        candidates=[distances[j] for j in js if j not in wallset and distances[j]>1e-16]
        assert candidates, '仅壁邻点的体积/面积回退尚未适配，不能用任意壁距代替'
        nearest.append(min(candidates))
    near=np.array(near);projection=abs(np.sum((p[near]-p[wallids])*unit,axis=1))
    assert min(projection)>0 and min(nearest)>0
    return dict(points=p,wallIDs=wallids,normalNeighbor=near,thermalDistance=projection,
                nearestWallDistance=np.array(nearest),nodeArea=np.linalg.norm(normal[wallids],axis=1))
