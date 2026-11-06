"""模块用途：生成具名实际分辨率半球柱网格和连续壁温表；只准备，不启动CFD。"""
import argparse
from pathlib import Path
import time

import numpy as np

from analyzeSU2DualVolumes import FACES
from naca3344Model import DATA, D, R, SHOULDER, END, meridian, source_data, wall_temperature, strip_mean
from prepareSU2NACALaminar import write_mesh
from su2ThermalRuntime import ROOT, sha
from validateIdealAir import save

PROFILE = 'naca-tn3344-quarter-3d-validation-mesh-v170'


def axes():
    """模块一：按电极/温度原站确定性细分，不扫描网格参数。"""
    stations = np.genfromtxt(str(DATA/'segments_v170/nominal_stations.csv'), delimiter=',', names=True,
                            dtype=None, encoding='utf8')
    required = np.sort(np.r_[0., SHOULDER, END, stations['s_m'], source_data()[1]['s_over_D']*D])
    # 同一个肩点可能经不同的浮点表达式重复，只合并1e-14m内的点。
    required = required[np.r_[True, np.diff(required)>1e-14]]
    s = np.r_[required[0], np.concatenate([np.linspace(a,b,int(np.ceil((b-a)/(D/64)))+1)[1:]
                                           for a,b in zip(required[:-1],required[1:])])]
    n = int(np.floor(np.log1p(2*D*.2/1e-6)/np.log(1.2)))
    y = 1e-6*np.expm1(np.arange(n+1)*np.log(1.2))/.2
    if 2*D-y[-1] < y[-1]-y[-2]:
        y = y[:-1]
    y = np.r_[y, 2*D]
    return s, y, np.linspace(0,np.pi/2,25), required


def build(s,y,phi):
    """模块二：真实轴线棱柱和其余六面体，所有边界从同一原节点表生成。"""
    ix=np.empty((len(y),len(s),len(phi)),int)
    points=[];labels=[]
    for j,offset in enumerate(y):
        for i,length in enumerate(s):
            xr=meridian(length,offset)
            for k,angle in enumerate(phi):
                if i==0 and k:
                    ix[j,i,k]=ix[j,i,0];continue
                ix[j,i,k]=len(points)
                points.append([xr[0],xr[1]*np.cos(angle),xr[1]*np.sin(angle)])
                labels.append((j,i,k))
    p=np.array(points);p[np.abs(p)<1e-17]=0
    labels=np.array(labels)
    # 三角函数的π/2舍入余量不能改变实际对称面；按已知原拓扑精确落在平面。
    p[labels[:,2]==0,2]=0
    p[labels[:,2]==len(phi)-1,1]=0
    j,k=np.meshgrid(np.arange(len(y)-1),np.arange(len(phi)-1),indexing='ij')
    wedges=np.column_stack([ix[j,0,0].ravel(),ix[j,1,k].ravel(),ix[j,1,k+1].ravel(),
                           ix[j+1,0,0].ravel(),ix[j+1,1,k].ravel(),ix[j+1,1,k+1].ravel()])
    j,i,k=np.meshgrid(np.arange(len(y)-1),np.arange(1,len(s)-1),np.arange(len(phi)-1),indexing='ij')
    hexes=np.column_stack([a.ravel() for a in (ix[j,i,k],ix[j,i+1,k],ix[j,i+1,k+1],ix[j,i,k+1],
                          ix[j+1,i,k],ix[j+1,i+1,k],ix[j+1,i+1,k+1],ix[j+1,i,k+1])])
    arrays={13:wedges,12:hexes}
    keys=[];oriented=[]
    for kind,rows in arrays.items():
        for face in FACES[kind]:
            ids=rows[:,face]
            keys.append(np.pad(np.sort(ids,axis=1),((0,0),(0,4-len(face))),constant_values=-1))
            oriented.append(np.pad(ids,((0,0),(0,4-len(face))),constant_values=-1))
    keys=np.vstack(keys);oriented=np.vstack(oriented)
    _,first,counts=np.unique(keys,axis=0,return_index=True,return_counts=True)
    assert np.all((counts==1)|(counts==2))
    markers={name:[] for name in ('wall','outer','outlet','sym_y','sym_z')}
    for ids in oriented[first[counts==1]]:
        ids=ids[ids>=0];jl,il,kl=labels[ids].T
        if np.all(jl==0):tag='wall'
        elif np.all(jl==len(y)-1):tag='outer'
        elif np.all(il==len(s)-1):tag='outlet'
        elif np.all(p[ids,2]==0):tag='sym_z'
        elif np.all(p[ids,1]==0):tag='sym_y'
        else:raise AssertionError('未分类完整外边界')
        markers[tag].append([5 if len(ids)==3 else 9]+ids.tolist())
    cells=[[kind]+row.tolist() for kind,rows in arrays.items() for row in rows]
    return p,arrays,cells,markers,labels,ix


def audit(p,arrays,markers,y):
    """模块三：全单元中位子体积和闭合边界积分，另保留光滑几何离散误差。"""
    volumes=[];min_sub=float('inf')
    for kind,ids in arrays.items():
        q=p[ids];center=q.mean(axis=1);v=np.zeros(len(q))
        for local in FACES[kind]:
            f=q[:,local];fc=f.mean(axis=1)
            for a,b in zip(np.moveaxis(f,1,0),np.moveaxis(np.roll(f,-1,axis=1),1,0)):
                part=np.einsum('ij,ij->i',a-center,np.cross(b-center,fc-center))/6
                min_sub=min(min_sub,float(part.min()));v+=part
        assert min_sub>0 and v.min()>0
        volumes.append(v)
    total=float(sum(v.sum() for v in volumes));bvol=0.;area=0.
    for tag,rows in markers.items():
        for length in (4,5):
            ids=np.array([row[1:] for row in rows if len(row)==length],int)
            if not len(ids):continue
            q=p[ids];fc=q.mean(axis=1)
            for a,b in zip(np.moveaxis(q,1,0),np.moveaxis(np.roll(q,-1,axis=1),1,0)):
                normal=np.cross(b-a,fc-a)/2
                bvol+=float(np.einsum('ij,ij->i',a,normal).sum()/3)
                if tag=='wall':area+=float(np.linalg.norm(normal,axis=1).sum())
    smooth=(2*np.pi/3*((R+y[-1])**3-R**3)+np.pi*((R+y[-1])**2-R**2)*(.127-R))/4
    smooth_area=(2*np.pi*R**2+2*np.pi*R*(.127-R))/4
    error=abs(bvol/total-1);ve=abs(total/smooth-1);ae=abs(area/smooth_area-1)
    assert error<=2e-10 and ve<=.001 and ae<=.001
    return dict(nodes=len(p),cells=sum(len(a) for a in arrays.values()),prisms=len(arrays[13]),
        hexes=len(arrays[12]),wallFaces=len(markers['wall']),minimumMedianSubVolume_m3=min_sub,
        fluidVolume_m3=total,boundaryVolumeRelativeError=error,smoothVolumeRelativeError=ve,
        wallAreaQuarter_m2=area,smoothWallAreaQuarter_m2=smooth_area,smoothWallAreaRelativeError=ae,
        physicalQualified=False)


def wall_audit(s,y,profile):
    """模块四：节点空间Tw精确表；旧分带方法仅作误差尺度，不当实际新边界。"""
    tw,extension=wall_temperature(s)
    with profile.open('w',encoding='ascii',newline='\n') as f:
        f.write('RAH_AXIAL_WALL_TEMPERATURE_V1,X_m,T_K\n')
        for x,t in zip(meridian(s)[::-1,0],tw[::-1]):f.write('%.17g,%.17g\n'%(x,t))
    means=np.array([strip_mean(lambda x:wall_temperature(x)[0],a,b)[0] for a,b in zip(s[:-1],s[1:])])
    ends=np.maximum(abs(means-tw[:-1]),abs(means-tw[1:]))
    jump=np.abs(np.diff(means))
    # .03W/(m K)只是明示的放大量纲尺，不用它代替求解器实际k或真实热流误差。
    amplification=.03*jump/y[1]
    rows=np.loadtxt(str(profile),delimiter=',',skiprows=1)
    node_error=float(np.max(abs(np.interp(meridian(s)[:,0],rows[:,0],rows[:,1])-tw)))
    mid=(s[:-1]+s[1:])/2
    between=float(np.max(abs(np.interp(meridian(mid)[:,0],rows[:,0],rows[:,1])-wall_temperature(mid)[0])))
    assert node_error<2e-12
    return dict(mode='single-wall-marker-spatial-node-table-preparation-only',nodes=len(tw),
        tableNodeMaxDifferenceK=node_error,midpointXInterpolationVsSInterpolationMaxDifferenceK=between,
        givenTwRangeK=[float(tw.min()),float(tw.max())],extensionNodes=int(extension.sum()),
        legacyPiecewiseMarkerMaxEndpointDifferenceK=float(ends.max()),legacyPiecewiseMaxSharedJumpK=float(jump.max()),
        illustrativeConductivityWmK=.03,legacyKDeltaTwOverFirstDistanceMaxWm2=float(amplification.max()),
        actualHeatFluxErrorKnown=False,runtimeSpatialProfileValidated=False,
        newRuntimeRequired='reuse-reviewed-v115-file-adapter-in-exact-observed-source; not-built-or-tested')


def prepare(folder):
    """模块五：只产出几何/列/壁温证据；不写可被误启动的request或case合同。"""
    started=time.monotonic();folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    folder.mkdir(parents=True,exist_ok=False)
    inputs=[Path(__file__).resolve(),ROOT/'scripts/naca3344Model.py',ROOT/'scripts/prepareSU2NACALaminar.py',
        ROOT/'scripts/analyzeSU2DualVolumes.py',DATA/'wall_temperature.csv',DATA/'source_manifest.json',
        DATA/'segments_v170/nominal_stations.csv',DATA/'segments_v170/nominal_segments.csv',
        ROOT/'docs/validation/v170/PLAN.md']
    locks={str(path):sha(path) for path in inputs}
    (folder/'executed_prepare.py').write_bytes(Path(__file__).read_bytes())
    save(folder/'input_manifest.json',locks)
    s,y,phi,required=axes();p,arrays,cells,markers,labels,ix=build(s,y,phi)
    assert len(cells)<=240000 and np.array_equal(p[0],[0,0,0])
    report=audit(p,arrays,markers,y)
    write_mesh(folder/'mesh.su2',p,cells,markers)
    np.savez_compressed(str(folder/'mesh_arrays.npz'),points=p,prisms=arrays[13],hexes=arrays[12],
                        nodeIndices=labels,columnNodeIDs=ix)
    metadata=dict(profile=PROFILE,s=s.tolist(),offsets=y.tolist(),azimuth=phi.tolist(),requiredS=required.tolist(),
        coordinates='nose +X, flow -X, SI; actual nose=(0,0,0)',singleWallMarker=True,
        physicalEndX_m=-.127,outerNormalOffset_m=2*D,firstDistance_m=float(y[1]),
        nominalGrowth=1.2,actualLastLayerRatio=float((y[-1]-y[-2])/(y[-2]-y[-3])),
        sourceSHA256=locks,readyForCFD=False,newPhysicalQualification=False)
    save(folder/'geometry.json',metadata)
    report['wallTemperature']=wall_audit(s,y,folder/'wall_temperature.csv')
    report['inputsStable']={path:sha(Path(path))==h for path,h in locks.items()}
    assert all(report['inputsStable'].values())
    report.update(elapsedSeconds=time.monotonic()-started,geometryCandidatePrepared=True,
        nativeTasks=0,readyForCFD=False,newPhysicalQualification=False)
    save(folder/'report.json',report)
    assert sum(path.stat().st_size for path in folder.iterdir() if path.is_file())<128*1024**2
    print({key:report[key] for key in ('nodes','cells','wallFaces','smoothVolumeRelativeError',
                                     'smoothWallAreaRelativeError','elapsedSeconds','readyForCFD')})
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--folder',required=True)
    prepare(parser.parse_args().folder)
