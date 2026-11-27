"""固定低速风洞的全周显热输运：轴对称速度、48 个周期方位、单值轴线。"""
import hashlib
import json
from io import StringIO
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu

from simpleLauncherHeatMesh import sampler
from simpleLauncherTemperature import meridian_operator, transfer
from simpleLauncherTunnel import CP, RHO, SPEED, T_AIR, T_WALL, geometry
from runEmbeddedSU2 import sha


SECTORS = 48
ANGLE = 2*np.pi/SECTORS
PHI = (np.arange(SECTORS)+0.5)*ANGLE


def domain(job, grid, reference):
    """在已登记的热网重建边界及速度投影，原温度只供独立对照。"""
    flow = job/grid
    heat = job/("heat_"+grid)
    file = flow/"restart.csv"
    receipt = json.loads((flow/"engine_process.json").read_bytes())
    assert receipt["status"]=="terminal" and receipt["exitCode"]==0 and all(receipt["inputsStable"].values())
    assert not (receipt["cancelled"] or receipt["timedOut"])
    assert sha(file)==receipt["outputSHA256"]["restart.csv"]==reference["flowSHA256"]
    assert sha(heat/"temperature.npz")==reference["heatMeshSHA256"]
    info = json.loads((flow/"condition.json").read_bytes())
    definition, corners, diameter, wall = geometry()
    assert (info["speedMps"],info["wallK"],info["airK"],info["surfaceSHA256"]) == (
        SPEED,T_WALL,T_AIR,definition["fileSHA256"])
    assert info["grid"]==grid and info["diameterM"]==diameter
    state = np.loadtxt(file,delimiter=",",skiprows=1)
    ids = np.arange(len(state)).reshape(info["layers"]+1,len(info["wall"]))
    np.testing.assert_array_equal(state[:,0],np.arange(len(state)))
    np.testing.assert_allclose(state[ids[0],1:3],info["wall"],atol=1e-13,rtol=0)
    assert np.isfinite(state).all() and np.max(abs(state[ids[0],4:6]))<1e-9
    with np.load(heat/"temperature.npz") as model:
        xy,tr,wet = (model[k] for k in ("xy","triangles","wall"))
    axis = np.flatnonzero(abs(xy[:,1])<1e-12)
    sides = np.vstack((tr[:,[0,1]],tr[:,[1,2]],tr[:,[2,0]]))
    edges, counts = np.unique(np.sort(sides,axis=1),axis=0,return_counts=True)
    border = edges[counts==1]
    far_edges = border[~(np.isin(border,axis).all(axis=1) | np.isin(border,wet).all(axis=1))]
    neighbours = {}
    for a,b in far_edges:
        neighbours.setdefault(a,[]).append(b)
        neighbours.setdefault(b,[]).append(a)
    ends = [k for k,values in neighbours.items() if len(values)==1]
    assert len(ends)==2 and all(k in axis for k in ends)
    chain = [min(ends,key=lambda k:xy[k,0])]
    while len(chain)<=len(far_edges):
        choices = [n for n in neighbours[chain[-1]] if len(chain)==1 or n!=chain[-2]]
        assert len(choices)==1
        chain.append(choices[0])
    far = np.asarray(chain)
    assert len(np.unique(far))==len(far) and far[-1] in ends
    quads = np.stack((ids[:-1,:-1],ids[:-1,1:],ids[1:,1:],ids[1:,:-1]),axis=-1).reshape(-1,4)
    flow_tri = np.vstack((quads[:,[0,1,2]],quads[:,[0,2,3]]))
    centroid = xy[tr].mean(axis=1)
    velocity = sampler(state[:,1:3],flow_tri,centroid)(state[:,4:6])
    pe = info["pecletD"]
    matrix,outward,inlet,_,projection = meridian_operator(xy,tr,wet,axis,far,velocity*centroid[:,1,None],pe)
    triangle_area = np.cross(xy[tr[:,1]]-xy[tr[:,0]],xy[tr[:,2]]-xy[tr[:,0]])/2
    assert np.all(triangle_area>0) and np.intersect1d(wet,inlet).size==0
    # 以三角面重心分配子午控制面积。角向导热系数为 A_i / (Pe r_i Δφ²)。
    dual_area = np.bincount(tr.ravel(),weights=np.repeat(triangle_area/3,3),minlength=len(xy))
    conductance = np.zeros(len(xy))
    off = np.ones(len(xy),bool);off[axis]=False
    conductance[off] = dual_area[off]/(pe*xy[off,1]*ANGLE**2)
    arc = np.r_[0.,np.cumsum(np.linalg.norm(np.diff(wall,axis=0),axis=1))]
    wet_arc = np.r_[0.,np.cumsum(np.linalg.norm(np.diff(xy[wet],axis=0),axis=1))]
    assert abs(arc[-1]-wet_arc[-1])<1e-8
    stations = np.interp(wet_arc,arc,np.arange(len(wall)))
    # 在原母线转折处保留确切的参数端点，以免求积误跨相邻的面带。
    rounded = np.rint(stations)
    same = abs(stations-rounded)<1e-9
    stations[same]=rounded[same]
    face_map = transfer(corners,wall,stations,diameter)
    return dict(xy=xy,tri=tr,wet=wet,axis=axis,far=far,inlet=inlet,outward=outward,matrix=matrix,
                conductance=conductance,pe=pe,volumePerRadian=triangle_area@centroid[:,1],
                projection=projection,stations=stations,corners=corners,wall=wall,faceMap=face_map,
                diameter=diameter,definition=definition,grid=grid,flowSHA256=reference["flowSHA256"])


def steady_fields(model, walls):
    """按角向变换分解周期线性方程，轴线只保留零模态。一次分解可处理多组壁温。"""
    n = len(model["xy"])
    wet,axis,inlet = (model[k] for k in ("wet","axis","inlet"))
    spectral_wall = np.fft.rfft(walls,axis=1)
    if np.any(np.ptp(walls[:,:,np.isin(wet,axis)],axis=1)>1e-13):
        raise ValueError("轴线上的给定壁温必须单值")
    coefficients = np.zeros((len(walls),SECTORS//2+1,n),dtype=complex)
    for mode in range(SECTORS//2+1):
        if not np.any(spectral_wall[:,mode]):
            continue
        fixed = np.union1d(wet,inlet)
        if mode:
            fixed = np.union1d(fixed,axis)
        free = np.ones(n,bool);free[fixed]=False
        matrix = model["matrix"]+diags(4*np.sin(mode*ANGLE/2)**2*model["conductance"])
        current = coefficients[:,mode,:]
        current[:,wet] = spectral_wall[:,mode]
        if mode:
            current[:,axis]=0.
        forcing = -(matrix[free][:,fixed]@current[:,fixed].T)
        factor = splu(matrix[free][:,free].tocsc())
        columns = np.column_stack((forcing.real,forcing.imag))
        solved = factor.solve(columns)
        current[:,free] = (solved[:,:len(walls)]+1j*solved[:,len(walls):]).T
    return np.fft.irfft(coefficients,n=SECTORS,axis=1)


def heat_rows(model, theta):
    """在物理空间逐扇区回算反力、轴线方程、相邻扇区交换及远界能量。"""
    conductance = model["conductance"]
    outgoing = (theta-np.roll(theta,-1,axis=0))*conductance
    rows = (model["matrix"]@theta.T).T + outgoing-np.roll(outgoing,1,axis=0)
    wet,inlet,far,axis = (model[k] for k in ("wet","inlet","far","axis"))
    free = np.ones(theta.shape[1],bool);free[np.union1d(wet,inlet)]=False
    shared = axis[free[axis]]
    free[axis]=False
    # 轴线只有一个自由度，检验相连的所有扇区贡献之和；壁面反力保留每块面积的分配。
    free_error = max(np.max(abs(rows[:,free])),np.max(abs(rows[:,shared].sum(axis=0))))
    fixed = np.union1d(wet,inlet)
    # 固定温度在自由行产生的右端项是残差标度。轴线按同一全局方程求和。
    boundary_rows = model["matrix"][:,fixed]@theta[:,fixed].T
    scale = max(np.max(abs(boundary_rows[free])),np.max(abs(boundary_rows[shared].sum(axis=1))))
    wall = rows[:,wet]
    remote = np.sum((theta[:,far[:-1]]+theta[:,far[1:]])*model["outward"]*.5)-rows[:,inlet].sum()
    return wall,dict(relativeFreeRowResidual=float(free_error/scale),maxFreeRowPerRadian=float(free_error),
        energyImbalance=float(abs(remote-wall.sum())/abs(wall.sum())),
        wallReleasePerRadianSum=float(wall.sum()),remoteReleasePerRadianSum=float(remote),
        minTheta=float(theta.min()),maxTheta=float(theta.max()),maxAxisSpread=float(np.ptp(theta[:,axis],axis=0).max()),
        maxFreeAngularExchangePerRadian=float(np.max(abs(outgoing[:,free]))))


def wall_results(model, reaction):
    """按扇区和原面面积保守映射；半径四等分的底盘另积分到原扇片范围。"""
    tri = model["corners"]
    area = np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    center = tri.mean(axis=1)
    sector = np.floor(np.mod(np.arctan2(center[:,2],center[:,1]),2*np.pi)/ANGLE).astype(int)
    scale = RHO*SPEED*CP*(T_WALL-T_AIR)*model["diameter"]**2*ANGLE
    face_map = model["faceMap"]
    knots = model["stations"]
    base_end = np.flatnonzero(np.diff(model["wall"][:,0])>0)[0]
    radial = knots[knots<=base_end]/base_end
    ring_weights = np.zeros((4,len(radial)))
    # 底面扇片横截宽随半径线性变化，分段积分 2r N_i(r)，r 以底半径归一化。
    for band in range(4):
        for i,(a,b) in enumerate(zip(radial[:-1],radial[1:])):
            left,right=max(a,band/4),min(b,(band+1)/4)
            if right<=left:
                continue
            for z in (-1/np.sqrt(3),1/np.sqrt(3)):
                r=(left+right+z*(right-left))/2
                fraction=(r-a)/(b-a)
                ring_weights[band,i:i+2] += (right-left)*r*np.array([1-fraction,fraction])
    assert abs(ring_weights.sum()-1)<1e-12
    q = np.empty(len(area))
    rings = np.zeros(4)
    for j in range(SECTORS):
        which = sector==j
        mapping = face_map[which]
        patch_area = np.asarray(mapping.T@area[which]).ravel()
        assert np.all(patch_area>0)
        flux = -reaction[j]*scale/patch_area
        q[which] = mapping@flux
        base_area = area[which&(tri[:,:,0]==0).all(axis=1)].sum()
        rings -= base_area*(ring_weights@flux[:len(radial)])
    outward = -area@q
    np.testing.assert_allclose(outward,reaction.sum()*scale,rtol=1e-10,atol=0)
    base = (tri[:,:,0]==0).all(axis=1)
    np.testing.assert_allclose(rings.sum(),-area[base]@q[base],rtol=1e-9,atol=1e-15)
    regions={}
    for group in model["definition"]["regions"]:
        face=slice(group["firstFace"]-1,group["lastFace"])
        regions[group["name"]]=float(-area[face]@q[face])
    return q,dict(outwardW=float(outward),regionW=regions,baseAnnulusW=rings.tolist())


def write_uniform_comparison(job, destination, primary, flow, heat):
    """在本次已通过的扩域流场上新算等温复核，并与该次原壁面表一起绑定。"""
    assert flow.resolve()==(job/"outer").resolve() and heat.resolve()==(job/"heat_outer").resolve()
    source=primary["runs"]["outer"]
    assert primary["accepted"] and sha(heat/"wall.csv")==source["wallSHA256"]==primary["wallSHA256"]
    model=domain(job,"outer",dict(flowSHA256=source["restartSHA256"],heatMeshSHA256=sha(heat/"temperature.npz")))
    boundary=np.ones((1,SECTORS,len(model["wet"])))
    field=steady_fields(model,boundary)[0]
    rows,account=heat_rows(model,field)
    flux,integral=wall_results(model,rows)
    wall=np.loadtxt(heat/"wall.csv",delimiter=",",skiprows=1)
    area=wall[:,1]
    rms=float(np.sqrt((area@(flux-wall[:,3])**2)/(area@wall[:,3]**2)))
    assert np.isfinite(flux).all() and rms<.02
    assert account["relativeFreeRowResidual"]<1e-8 and account["energyImbalance"]<.005
    assert account["minTheta"]>=-1e-12 and account["maxTheta"]<=1+1e-12 and account["maxAxisSpread"]<1e-12
    assert abs(integral["outwardW"]/primary["outwardW"]-1)<.01
    for name,power in integral["regionW"].items():
        assert abs(power/primary["regions"][name]["outwardW"]-1)<.02
    local=destination.parent/"angular_uniform"
    local.mkdir(exist_ok=False)
    np.savez_compressed(local/"field.npz",theta=field,wallReactions=rows,sectorCentreRad=PHI,wet=model["wet"])
    text=StringIO()
    np.savetxt(text,np.column_stack((wall[:,0],area,flux)),fmt=["%d","%.17g","%.17g"],delimiter=",",
               header="faceID,areaM2,heatFluxWm2",comments="")
    csv_data=text.getvalue().encode("utf8")
    summary=dict(sectors=SECTORS,faces=len(area),accepted=True,condition=primary["condition"],
                 surfaceSHA256=primary["surfaceSHA256"],flowSHA256=source["restartSHA256"],
                 primaryWallSHA256=source["wallSHA256"],heatGridSHA256=sha(heat/"temperature.npz"),
                 fieldSHA256=sha(local/"field.npz"),csvSHA256=hashlib.sha256(csv_data).hexdigest(),
                 freeAndBoundaryDOF=SECTORS*(len(model["xy"])-len(model["axis"]))+len(model["axis"]),
                 rmsDifference=rms,**account,**integral)
    with ZipFile(destination,"x",ZIP_DEFLATED) as bundle:
        bundle.writestr("comparison.json",json.dumps(summary,indent=2,allow_nan=False))
        bundle.writestr("angular_wall.csv",csv_data)
    return summary
