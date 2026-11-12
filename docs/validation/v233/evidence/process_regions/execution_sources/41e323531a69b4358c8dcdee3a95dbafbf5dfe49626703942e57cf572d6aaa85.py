"""模块职责：从既有四面体提取封闭区域及特征边界；不推进流动或授予外缘资格。"""
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy


def read_tetra_mesh(path):
    """模块一：读取明确的SU2四面体/三角形格式，保持母节点、单元和marker编号。"""
    rows=Path(path).read_text(encoding='utf8').splitlines()
    def section(key):
        found=[i for i,s in enumerate(rows) if s.split('=')[0].strip()==key]
        if len(found)!=1: raise ValueError('REGION_MESH_SECTION_'+key)
        j=found[0];return j,int(rows[j].split('=')[1].split()[0])
    j,n=section('NELEM');cells=[]
    for s in rows[j+1:j+1+n]:
        w=list(map(int,s.split()))
        if w[0]!=10 or len(w) not in (5,6):raise ValueError('REGION_REQUIRES_TETRA')
        if len(w)==6 and w[-1]!=len(cells):raise ValueError('REGION_CELL_ID')
        cells.append(w[1:5])
    j,n=section('NPOIN');xyz=np.empty((n,3));seen=set()
    for s in rows[j+1:j+1+n]:
        w=s.split();i=int(w[3])
        if i in seen or not 0<=i<n:raise ValueError('REGION_NODE_ID')
        seen.add(i);xyz[i]=list(map(float,w[:3]))
    if len(seen)!=n or not np.isfinite(xyz).all():raise ValueError('REGION_COORDINATES')
    markers={}
    for j,s in enumerate(rows):
        if s.split('=')[0].strip()!='MARKER_TAG':continue
        tag=s.split('=')[1].strip();n=int(rows[j+1].split('=')[1]);faces=[]
        if tag in markers:raise ValueError('REGION_DUPLICATE_MARKER')
        for line in rows[j+2:j+2+n]:
            w=list(map(int,line.split()))
            if len(w)!=4 or w[0]!=5:raise ValueError('REGION_REQUIRES_TRIANGLE')
            faces.append(w[1:])
        markers[tag]=np.asarray(faces,dtype=np.int64)
    cells=np.asarray(cells,dtype=np.int64)
    if np.any(cells<0) or np.any(cells>=len(xyz)):raise ValueError('REGION_CELL_NODE_RANGE')
    return xyz,cells,markers


def topology(xyz,cells,markers):
    """模块二：面规范键建立原单元邻接，法向统一从owner指向区域外部。"""
    n=len(cells);half=np.concatenate([cells[:,ids] for ids in ((0,1,2),(0,3,1),(1,3,2),(2,3,0))])
    half_owner=np.tile(np.arange(n),4);keys,first,inv,count=np.unique(np.sort(half,axis=1),axis=0,return_index=True,return_inverse=True,return_counts=True)
    if np.any(count>2):raise ValueError('REGION_NONMANIFOLD')
    faces=half[first];owner=half_owner[first];other=np.full(len(faces),-1,dtype=np.int64)
    duplicate=np.arange(len(inv))!=first[inv];other[inv[duplicate]]=half_owner[duplicate]
    center=xyz[cells].mean(axis=1);fc=xyz[faces].mean(axis=1)
    normal=np.cross(xyz[faces[:,1]]-xyz[faces[:,0]],xyz[faces[:,2]]-xyz[faces[:,0]])/2
    sign=np.einsum('ij,ij->i',normal,fc-center[owner])
    if np.any(sign==0):raise ValueError('REGION_ZERO_FACE_OR_CELL')
    normal[sign<0]*=-1
    delta=xyz[cells[:,1:]]-xyz[cells[:,0,None]]
    volume=abs(np.einsum('ij,ij->i',delta[:,0],np.cross(delta[:,1],delta[:,2])))/6
    if np.any(volume<=0):raise ValueError('REGION_NONPOSITIVE_VOLUME')
    dtype=np.dtype([('a',np.int64),('b',np.int64),('c',np.int64)])
    ordered=np.ascontiguousarray(keys).view(dtype).ravel();marker_ids={}
    for tag,f in markers.items():
        k=np.sort(f,axis=1);indices=np.searchsorted(ordered,np.ascontiguousarray(k).view(dtype).ravel())
        if np.any(indices>=len(keys)) or not np.array_equal(keys[indices],k):raise ValueError('REGION_MARKER_NOT_FOUND')
        if np.any(other[indices]>=0):raise ValueError('REGION_MARKER_INTERNAL')
        marker_ids[tag]=indices
    mapped=np.concatenate(list(marker_ids.values()))
    if len(np.unique(mapped))!=len(mapped) or not np.array_equal(np.sort(mapped),np.flatnonzero(other<0)):raise ValueError('REGION_BOUNDARY_COVERAGE')
    internal=other>=0
    graph=coo_matrix((np.ones(2*sum(internal)),(np.r_[owner[internal],other[internal]],np.r_[other[internal],owner[internal]])),shape=(n,n)).tocsr()
    return dict(faces=faces,owner=owner,other=other,normal=normal,center=center,faceCenter=fc,volume=volume,markers=marker_ids,graph=graph)


def extract_region(topo,selected):
    """模块三：完整原单元集合的边界与几何闭合；不把此结果当原生节点双控制体。"""
    selected=np.asarray(selected,dtype=bool);owner=topo['owner'];other=topo['other'];internal=other>=0
    inside=selected[owner];neighbor=np.zeros(len(other),dtype=bool);neighbor[internal]=selected[other[internal]]
    boundary=np.flatnonzero(inside!=neighbor);orientation=np.where(inside[boundary],1.,-1.)
    normal=topo['normal'][boundary]*orientation[:,None];area=np.linalg.norm(normal,axis=1)
    vector_error=float(np.linalg.norm(normal.sum(axis=0))/area.sum())
    cells=np.flatnonzero(selected);volume=float(topo['volume'][cells].sum())
    surface_volume=float(np.sum(topo['faceCenter'][boundary]*normal)/3)
    volume_error=abs(surface_volume-volume)/volume
    if vector_error>1e-11 or volume_error>1e-11:raise ValueError('REGION_GEOMETRIC_CLOSURE')
    components=int(connected_components(topo['graph'][cells][:,cells],directed=False,return_labels=False))
    return dict(cells=cells,boundary=boundary,normal=normal,interface=boundary[internal[boundary]],
      summary=dict(cells=len(cells),components=components,boundaryFaces=len(boundary),interfaceFaces=int(sum(internal[boundary])),
      volumeM3=volume,surfaceDivergenceVolumeM3=surface_volume,vectorAreaClosureRelative=vector_error,volumeClosureRelative=volume_error,
      nativeDualControlVolume=False,physicalQualification=False))


def interface_characteristics(U,faces,normal):
    """模块四：三角形顶点守恒量均值反算同气体，判定五个Euler特征的入/出射数。"""
    if U.ndim!=2 or U.shape[1]!=5 or not np.isfinite(U).all():raise ValueError('REGION_STATE_SHAPE')
    u=U[faces].mean(axis=1);rho=u[:,0]
    if np.any(rho<=0):raise ValueError('REGION_STATE_DENSITY')
    velocity=u[:,1:4]/rho[:,None];e=u[:,4]/rho-.5*np.sum(velocity**2,axis=1)-300000
    air=PhaseAwareViscousAir();t=inverse_energy(air,e);p=rho*air.R*t;state=air.evaluate(p,t)
    area=np.linalg.norm(normal,axis=1)
    if np.any(area<=0):raise ValueError('REGION_INTERFACE_NORMAL')
    un=np.einsum('ij,ij->i',velocity,normal)/area;a=np.asarray(state['soundSpeedMps'])
    eigen=np.column_stack([un-a,un,un,un,un+a]);incoming=np.sum(eigen<0,axis=1);sonic=np.any(eigen==0,axis=1)
    h=np.asarray(state['enthalpyJkg'])+.5*np.sum(velocity**2,axis=1)
    return dict(velocity=velocity,p=p,T=t,H=h,S=np.asarray(state['entropyJkgK']),normalMach=un/a,
      incoming=incoming,sonic=sonic,area=area,conservativeAverage=u,
      summary=dict(faces=len(faces),areaM2=float(area.sum()),normalMachRange=[float(min(un/a)),float(max(un/a))],
        characteristicFaceCounts={str(i):int(sum(incoming==i)) for i in range(6)},
        characteristicAreasM2={str(i):float(sum(area[incoming==i])) for i in range(6)},
        sonicFaces=int(sum(sonic)),method='conservative vertex mean / same phase Air / geometric triangular-face normal',
        nativeCharacteristicBoundary=False,nativeNumericalFlux=False,physicalQualification=False))
