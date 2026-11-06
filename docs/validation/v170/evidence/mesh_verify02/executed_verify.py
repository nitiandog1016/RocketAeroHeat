"""模块用途：独立标准形函数/面方向/实验锚点审核实际网格；不启动求解器。"""
import argparse
from pathlib import Path
import sys
import time
import itertools
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from su2ThermalRuntime import read,sha
from su2NACALaminarInterface import mesh
from naca3344Model import DATA,D,R,SHOULDER,END,wall_temperature
from validateIdealAir import save

# 独立标准单元拓扑，不复用生成器FACES。
FACE={12:((0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)),
      13:((0,2,1),(3,4,5),(0,1,4,3),(1,2,5,4),(2,0,3,5))}
SIGNS=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])


def hex_jac(q,x):
    grad=np.empty((8,3))
    for j in range(3):
        other=[k for k in range(3) if k!=j]
        grad[:,j]=SIGNS[:,j]*np.prod(1+SIGNS[:,other]*x[other],axis=1)/8
    return np.linalg.det(np.einsum('nki,kj->nij',q,grad))


def verify(folder,out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    (out/'executed_verify.py').write_bytes(Path(__file__).read_bytes())
    paths=[folder/name for name in ('mesh.su2','mesh_arrays.npz','geometry.json','report.json','wall_temperature.csv')]
    paths += [Path(__file__).resolve(),ROOT/'scripts/prepareNACAValidationMesh.py']
    locks={str(p):sha(p) for p in paths};save(out/'input_manifest.json',locks)
    p,cells,markers=mesh(folder/'mesh.su2');meta=read(folder/'geometry.json');report=read(folder/'report.json')
    with np.load(str(folder/'mesh_arrays.npz')) as a:
        assert np.array_equal(p,a['points'])
        labels=a['nodeIndices'];ix=a['columnNodeIDs']
        arrays={k:np.array([c[1:] for c in cells if c[0]==k]) for k in (12,13)}
        assert np.array_equal(arrays[12],a['hexes']) and np.array_equal(arrays[13],a['prisms'])
    s=np.array(meta['s']);y=np.array(meta['offsets']);phi=np.array(meta['azimuth'])
    assert np.array_equal(p[0],[0,0,0]) and abs(y[1]-1e-6)<1e-18 and y[-1]==2*D
    assert len(phi)==25 and max(np.diff(s))<=D/64+1e-14
    assert len(cells)<=240000 and np.all(p[:,1:]>=0)
    assert np.max(abs(p[ix[:,0,0],0]-y))<2e-16
    assert np.array_equal(ix[:,:,0][0],np.unique(ix[0],axis=0)[:,0])
    required=np.r_[np.genfromtxt(str(DATA/'segments_v170/nominal_stations.csv'),delimiter=',',names=True)['s_m'],
                   np.genfromtxt(str(DATA/'wall_temperature.csv'),delimiter=',',names=True)['s_over_D']*D,SHOULDER,END]
    anchor=max(min(abs(s-x)) for x in required);assert anchor<1e-14
    # 精确法向列和单元角度范围；张量积解析下界可保证全参考域正性，角点独立核对。
    j,i,k=labels.T;theta=np.minimum(s[i]/R,np.pi/2)
    expected=np.column_stack(((R+y[j])*np.cos(theta)-R-np.maximum(s[i]-SHOULDER,0),
                               (R+y[j])*np.sin(theta)*np.cos(phi[k]),(R+y[j])*np.sin(theta)*np.sin(phi[k])))
    coordinate_error=float(np.max(abs(p-expected)));assert coordinate_error<1e-15
    q=p[arrays[12]];j,i,k=labels[arrays[12][:,0]].T
    expected_hex=np.column_stack((ix[j,i,k],ix[j,i+1,k],ix[j,i+1,k+1],ix[j,i,k+1],
                                 ix[j+1,i,k],ix[j+1,i+1,k],ix[j+1,i+1,k+1],ix[j+1,i,k+1]))
    assert np.array_equal(arrays[12],expected_hex)
    dr=R+y[j];dy=y[j+1]-y[j];ds=s[i+1]-s[i];dp=phi[k+1]-phi[k]
    bound=np.where(s[i+1]<=SHOULDER+1e-14,dy*dr**2*np.sin(s[i]/R)*np.sin(ds/R)*np.sin(dp)/8,
                   ds*dy*dr*np.sin(dp)/8)
    corner=np.min([hex_jac(q,np.array(x,float)) for x in itertools.product((-1,1),repeat=3)],axis=0)
    bound_error=float(np.max(abs(corner/bound-1)));assert min(bound)>0 and bound_error<2e-8
    gx,gw=np.polynomial.legendre.leggauss(2);volume=0.
    for aa,bb,cc in itertools.product(range(2),repeat=3):
        value=hex_jac(q,np.array([gx[aa],gx[bb],gx[cc]]));assert value.min()>0
        volume+=float(value.sum()*gw[aa]*gw[bb]*gw[cc])
    q=p[arrays[13]];bottom=q[:,:3];delta=q[:,3:]-bottom
    j,_,k=labels[arrays[13][:,1]].T
    expected_prism=np.column_stack((ix[j,0,0],ix[j,1,k],ix[j,1,k+1],ix[j+1,0,0],ix[j+1,1,k],ix[j+1,1,k+1]))
    assert np.array_equal(arrays[13],expected_prism)
    # 同心径向伸缩使J与三角形局部坐标无关、沿层厚单调，下面为全参考域下界。
    scale=(R+y[j+1])/(R+y[j])
    homothetic=np.max(abs((q[:,3:]+np.array([R,0,0]))-(q[:,:3]+np.array([R,0,0]))*scale[:,None,None]))
    assert homothetic<1e-15
    prism_bound=(y[j+1]-y[j])*(R+y[j])**2*np.sin(s[1]/R)**2*np.sin(phi[k+1]-phi[k])
    prism_corner=np.einsum('ij,ij->i',bottom[:,1]-bottom[:,0],np.cross(bottom[:,2]-bottom[:,0],delta[:,0]))
    prism_bound_error=float(np.max(abs(prism_corner/prism_bound-1)))
    assert prism_bound.min()>0 and prism_bound_error<2e-8
    prism_min=float('inf')
    # 独立Duffy三角积分及正性：底三角方向不变，圆锥伸缩半径始终正。
    for aa,bb,cc in itertools.product(range(2),repeat=3):
        u=(gx[aa]+1)/2;v=(gx[bb]+1)/2*(1-u);z=(gx[cc]+1)/2
        ds=bottom[:,1]-bottom[:,0]+z*(delta[:,1]-delta[:,0])
        dt=bottom[:,2]-bottom[:,0]+z*(delta[:,2]-delta[:,0])
        dz=delta[:,0]*(1-u-v)+delta[:,1]*u+delta[:,2]*v
        jac=np.einsum('ij,ij->i',ds,np.cross(dt,dz));prism_min=min(prism_min,float(jac.min()))
        assert jac.min()>0
        volume+=float(jac.sum()*gw[aa]*gw[bb]*gw[cc]*(1-u)/8)
    err=abs(volume/report['fluidVolume_m3']-1);assert err<2e-10
    # 旋转到最小节点开头；每个内面必须恰有一正一反的相同环，外面必须对应marker。
    allkeys=[];allsigns=[]
    for kind,ids in arrays.items():
        for local in FACE[kind]:
            f=ids[:,local];n=f.shape[1];m=f.argmin(axis=1)
            rot=np.take_along_axis(f,(m[:,None]+np.arange(n))%n,axis=1)
            direction=np.where(rot[:,1]<rot[:,-1],1,-1)
            reverse=rot[:,np.r_[0,np.arange(n-1,0,-1)]]
            can=np.where((direction==1)[:,None],rot,reverse)
            allkeys.append(np.pad(can,((0,0),(0,4-n)),constant_values=-1));allsigns.append(direction)
    keys=np.vstack(allkeys);sign=np.concatenate(allsigns)
    unique,inv,count=np.unique(keys,axis=0,return_inverse=True,return_counts=True)
    sums=np.bincount(inv,weights=sign);assert np.all((count==1)|(count==2)) and np.all(sums[count==2]==0)
    boundary=[]
    for rows in markers.values():
        for row in rows:
            f=row[1:];m=f.index(min(f));f=f[m:]+f[:m]
            if f[1]>f[-1]:f=[f[0]]+f[:0:-1]
            boundary.append(f+[-1]*(4-len(f)))
    assert np.array_equal(unique[count==1],np.unique(np.array(boundary),axis=0))
    assert sum(count==1)==len(boundary)
    wall_ids=np.unique([v for f in markers['wall'] for v in f[1:]])
    assert np.array_equal(wall_ids,np.unique(ix[0]))
    # 独立检查原生法向邻点选择和NormalDistance的5%下限，不能由1μm输入直接推定。
    from nacaWallFlux import geometry as native_wall_geometry
    normal_geometry=native_wall_geometry(p,cells,{'wall_check':markers['wall']})['wall_check']
    ids=normal_geometry['ids'];jl,il,kl=labels[ids].T
    expected_neighbor=ix[jl+1,il,kl]
    assert np.array_equal(normal_geometry['normalNeighbor'],expected_neighbor)
    euclidean=np.linalg.norm(p[expected_neighbor]-p[ids],axis=1)
    cosine=normal_geometry['distance']/euclidean
    assert cosine.min()>.05
    table=np.loadtxt(str(folder/'wall_temperature.csv'),delimiter=',',skiprows=1)
    tw=wall_temperature(s[labels[wall_ids,1]])[0]
    error_tw=float(np.max(abs(np.interp(p[wall_ids,0],table[:,0],table[:,1])-tw)))
    assert error_tw<2e-12 and len(markers)==5 and not (folder/'request.json').exists()
    assert not meta['readyForCFD'] and not report['wallTemperature']['runtimeSpatialProfileValidated']
    stable={path:sha(Path(path))==h for path,h in locks.items()};assert all(stable.values())
    result=dict(passed=True,nodes=len(p),cells=len(cells),wallFaces=len(markers['wall']),
        completeInternalFaces=int(sum(count==2)),completeBoundaryFaces=len(boundary),
        volumeIndependentRelativeError=err,hexaAnalyticLowerBoundMinimum=float(min(bound)),
        hexaIndependentBoundMaxRelativeDifference=bound_error,prismGaussMinimum=prism_min,
        prismFullReferenceLowerBoundMinimum=float(prism_bound.min()),prismIndependentBoundRelativeError=prism_bound_error,
        tensorConnectivityExact=True,normalNeighborMatchesFirstColumn=True,normalDistanceFloorTriggered=0,
        minimumNormalDistanceCosine=float(cosine.min()),normalDistanceRangeM=[float(normal_geometry['distance'].min()),float(normal_geometry['distance'].max())],
        coordinateMaxDifferenceM=coordinate_error,requiredStationMaxDistanceM=float(anchor),
        wallNodeTemperatureMaxDifferenceK=error_tw,wallNodes=len(wall_ids),inputSHA256=locks,inputsStable=stable,
        elapsedSeconds=time.monotonic()-started,readyForCFD=False,newPhysicalQualification=False)
    save(out/'report.json',result);print({k:v for k,v in result.items() if k not in ('inputSHA256','inputsStable')})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();verify(Path(args.folder).resolve(),Path(args.output).resolve())
