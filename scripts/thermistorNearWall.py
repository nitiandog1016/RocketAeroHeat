"""在 B 型热敏珠的全周体网内平分靠壁区间，并按原始节点号传递启动场。"""
from pathlib import Path
import numpy as np
from cm011TunnelSector import HEX_FACES, PRISM_FACES, write_volume
from thermistorStartupSurface import columns, read, write
from thermistorBead import sha

LAYERS=48
SPLIT=8
FIELDS=["PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"]


def read_grid(path):
    """按热敏珠已有 SU2 写法读取原序单元、节点和真实边界。"""
    with path.open(encoding="ascii") as f:
        def header(key):
            name,val=next(f).strip().split("=",1)
            assert name==key
            return val.strip()
        assert header("NDIME")=="3"
        cells=[]
        for _ in range(int(header("NELEM"))):
            kind,*nodes=map(int,next(f).split())
            assert (kind,len(nodes)) in ((13,6),(12,8))
            cells.append(nodes)
        xyz=np.loadtxt([next(f) for _ in range(int(header("NPOIN")))])
        markers={}
        for _ in range(int(header("NMARK"))):
            name=header("MARKER_TAG")
            assert name not in markers
            faces=[]
            for _ in range(int(header("MARKER_ELEMS"))):
                kind,*nodes=map(int,next(f).split())
                assert (kind,len(nodes)) in ((5,3),(9,4))
                faces.append(nodes)
            markers[name]=faces
        assert set(markers)=={"wall","far"} and not f.read().strip()
    return xyz,cells,markers


def refine(xyz,ids,cells):
    """按原射线共享中点，同一个新节点服务相邻单元、两极及周向接缝。"""
    grid=ids.reshape(49,LAYERS+1,65)
    assert np.array_equal(grid[0],grid[-1])
    coarse_in_fine=np.arange(LAYERS+1)+np.minimum(np.arange(LAYERS+1),SPLIT)
    fine=np.empty((49,LAYERS+SPLIT+1,65),int);fine[:,coarse_in_fine]=grid
    mid_by_lower=np.full(len(xyz),-1,int)
    pairs=[];offset=len(xyz)
    for j in range(SPLIT):
        edges=np.unique(np.stack((grid[:,j],grid[:,j+1]),axis=-1).reshape(-1,2),axis=0)
        assert len(edges)==3026 and len(np.unique(edges[:,0]))==len(edges) and np.all(mid_by_lower[edges[:,0]]<0)
        mid_by_lower[edges[:,0]]=np.arange(offset,offset+len(edges))
        fine[:,2*j+1]=mid_by_lower[grid[:,j]]
        pairs.append(edges);offset+=len(edges)
    pairs=np.concatenate(pairs)
    new_xyz=np.concatenate((xyz,xyz[pairs].mean(axis=1)))
    assert np.array_equal(fine[0],fine[-1]) and len(np.unique(fine))==len(new_xyz)==172482
    out=[];parent=[]
    for i,cell in enumerate(cells):
        row=np.asarray(cell)
        lower,upper=([0,1,2],[3,4,5]) if len(cell)==6 else ([0,1,4,5],[3,2,7,6])
        mid=mid_by_lower[row[lower]]
        if mid[0]<0:
            assert np.all(mid<0)
            out.append(cell);parent.append(i)
        else:
            assert np.array_equal(pairs[mid-len(xyz),1],row[upper])
            inner=row.copy();outer=row.copy()
            inner[upper]=mid;outer[lower]=mid
            out.extend([inner.tolist(),outer.tolist()]);parent.extend([i,i])
    assert len(out)==172032
    return new_xyz,fine,np.asarray(parent),pairs,out


def geometry(xyz,cells,markers,shell_mask):
    """求各棱锥子体积和顶点中位对偶体积，同时验证内面的双侧定向与边界。"""
    dual=np.zeros(len(xyz));shell_dual=np.zeros(len(xyz));volume=np.zeros(len(cells));minimum=np.inf
    widths=np.fromiter(map(len,cells),int)
    faces={3:[],4:[]}
    for width,local in ((6,PRISM_FACES),(8,HEX_FACES)):
        sel=np.flatnonzero(widths==width)
        nodes=np.asarray([cells[i] for i in sel],np.int32)
        q=xyz[nodes];center=q.mean(axis=1)
        inner=shell_mask[sel]
        for face in local:
            vertex=nodes[:,face];faces[len(face)].append(vertex)
            fc=xyz[vertex].mean(axis=1)-center
            for a,b in zip(face,(*face[1:],face[0])):
                part=np.einsum("ij,ij->i",q[:,a]-center,np.cross(q[:,b]-center,fc))/6
                minimum=min(minimum,float(part.min())/2)
                assert np.all(part>0),"发现非正定向子体积"
                volume[sel]+=part
                for end in (a,b):
                    np.add.at(dual,nodes[:,end],part*.5)
                    np.add.at(shell_dual,nodes[inner,end],part[inner]*.5)
    boundary=sum(markers.values(),[])
    interior_count=0
    for width,blocks in faces.items():
        face=np.vstack(blocks)
        tags=np.asarray([f for f in boundary if len(f)==width],dtype=np.int32).reshape(-1,width)
        # 以同一顶点为起点检查两侧环绕方向；真实边界的方向也须与标记一致。
        both=np.vstack([face,tags]);idx=np.argmin(both,axis=1);line=np.arange(len(both))
        sign=np.where(both[line,(idx+1)%width]<both[line,(idx-1)%width],1,-1)
        sign[len(face):]*=-1
        _,key,count=np.unique(np.sort(both,axis=1),axis=0,return_counts=True,return_inverse=True)
        assert np.all(count==2) and np.all(np.bincount(key,weights=sign)==0)
        assert len(np.unique(key[len(face):]))==len(tags)
        interior_count+=len(count)-len(tags)
    assert abs(dual.sum()/volume.sum()-1)<2e-12 and dual.min()>0
    return dual,shell_dual,volume,dict(minimumSignedSubvolumeD3=minimum,interiorFaces=interior_count,
                            boundaryFaces={k:len(v) for k,v in markers.items()},volumeD3=float(volume.sum()))


def prepare(source,native_file,comparison,folder):
    original,cells,markers=read_grid(source/"mesh.su2")
    with np.load(source/"nodes.npz") as z:
        assert np.array_equal(z["xyz"],original)
        ids=z["ids"].reshape(49,49,65)
    xyz,fine,parent,pairs,new_cells=refine(original,ids,cells)
    write_volume(folder,xyz,new_cells,markers)
    stored,stored_cells,stored_markers=read_grid(folder/"mesh.su2")
    assert np.array_equal(stored,xyz) and stored_cells==new_cells and stored_markers==markers
    shell=np.bincount(parent,minlength=len(cells))==2
    assert shell.sum()==SPLIT*48*64
    base_v,base_shell,base_cell,base_geom=geometry(original,cells,markers,shell)
    new_v,new_shell,new_cell,new_geom=geometry(stored,stored_cells,stored_markers,shell[parent])
    volume_error=float(np.max(abs(np.bincount(parent,weights=new_cell)/base_cell-1)))
    native=columns(native_file,["id","volume"]);native=native[np.argsort(native[:,0])]
    assert np.array_equal(native[:,0],np.arange(len(original)))
    dual_error=float(max(abs(base_v/native[:,1]-1)))
    assert volume_error<2e-10 and dual_error<2e-10 and abs(new_v.sum()/base_v.sum()-1)<2e-12
    old=columns(source/"seed_00000.csv",FIELDS)
    assert np.array_equal(old[:,0],np.arange(len(original)))
    np.testing.assert_allclose(old[:,1:4],original,rtol=1e-14,atol=2e-14)
    values=np.vstack((old[:,4:],old[pairs,4:].mean(axis=1)))
    np.savetxt(folder/"seed_00000.csv",np.column_stack((np.arange(len(xyz)),xyz,values)),delimiter=",",
               header=",".join('"'+key+'"' for key in FIELDS),comments="",fmt=["%d"]+["%.17g"]*9)
    check=columns(folder/"seed_00000.csv",FIELDS)
    assert np.array_equal(check[:,0],np.arange(len(xyz))) and np.array_equal(check[:len(old),4:],old[:,4:])
    assert np.array_equal(check[:,1:4],xyz) and np.array_equal(check[:,4:],values)
    assert np.isfinite(values).all() and values[:,-1].min()>=0 and values[:,-1].max()<=1
    assert np.all(values[fine[:,0],1:4]==0) and np.all(values[fine[:,0],-1]==1)
    assert np.array_equal(fine[:,-1],ids[:,-1]) and np.array_equal(fine[:,2],ids[:,1])
    g=read(source/"geometry.json");diam=g["transverseDiameterM"]
    wall=np.unique(ids[:,0])
    # 按实际角、母线点排列，包括两极各一次。
    _,i=np.unique(ids[:,0].ravel(),return_index=True)
    first=ids[:,1].ravel()[i];middle=fine[:,1].ravel()[i]
    edge=(xyz[first]-xyz[wall])*diam
    thick=np.linalg.norm(edge,axis=1)
    p=xyz[wall];a=np.max(abs(p[:,0]));b=np.max(np.hypot(p[:,1],p[:,2]))
    normal=p/np.array([a*a,b*b,b*b]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    height=np.sum(normal*edge,axis=1)
    height_new=np.sum(normal*(xyz[middle]-p),axis=1)*diam
    np.testing.assert_allclose(height_new,height/2,rtol=1e-11,atol=0)
    dt=g["nominalWallK"]-g["airK"]
    slope=dt*(values[first,-1]-values[wall,-1])/thick
    half_slope=dt*(values[middle,-1]-values[wall,-1])/np.linalg.norm((xyz[middle]-p)*diam,axis=1)
    slope_err=float(max(abs(half_slope-slope))/max(abs(slope)))
    scale=g["rho"]*g["cp"]*dt*diam**3
    old_energy=float(base_v@old[:,-1])*scale;new_energy=float(new_v@values[:,-1])*scale
    shell_old=float(base_shell@old[:,-1])*scale;shell_new=float(new_shell@values[:,-1])*scale
    assert abs((new_energy-old_energy)-(shell_new-shell_old))<1e-16
    reference=read(comparison);interval_heat=reference["secondHalfEnergy"]["wallReleaseW"]*reference["elapsedPhysicalTimeS"]
    err=(new_energy-old_energy)/interval_heat
    admitted=bool(slope_err<1e-10 and abs(err)<.01)
    np.savez_compressed(folder/"nodes.npz",xyz=xyz,ids=fine.reshape(49,-1),originalFirst=first,refinedFirst=middle,
                        oldGridNodes=len(original),dualVolume=new_v)
    write(folder/"grid.json",dict(admitted=admitted,oldNodes=len(original),nodes=len(xyz),oldCells=len(cells),cells=len(new_cells),
          oldLayers=48,newLayers=56,splitIntervals=8,originalNodesPreserved=True,matchedOriginalFirstLayer=2,
          oldGeometry=base_geom,newGeometry=new_geom,parentVolumeMaxRelativeDifference=volume_error,
          originalVsNativeDualMaxRelativeDifference=dual_error,originalRaySpacingM=[float(thick.min()),float(thick.max())],
          originalNormalSpacingM=[float(height.min()),float(height.max())],newNormalSpacingM=[float(height_new.min()),float(height_new.max())],
          temperatureRaySlopeMaxRelativeChange=slope_err,originalExcessEnergyJ=old_energy,refinedExcessEnergyJ=new_energy,
          excessEnergyChangeJ=new_energy-old_energy,originalInnerShellEnergyJ=shell_old,refinedInnerShellEnergyJ=shell_new,
          referenceIntervalHeatJ=interval_heat,energyChangeFractionOfWallHeat=err,
          scalarRange=[float(values[:,-1].min()),float(values[:,-1].max())],sourceSeedSHA256=sha(source/"seed_00000.csv"),
          gridSHA256=sha(folder/"mesh.su2"),seedSHA256=sha(folder/"seed_00000.csv"),nodeMapSHA256=sha(folder/"nodes.npz")))
