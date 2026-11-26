"""低速、常物性、零攻角风洞算例的守恒显热输运与原面投影。"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

from cm011Tunnel import volume, write_json
from runEmbeddedSU2 import sha
from simpleLauncherTunnel import CP, RHO, SPEED, T_AIR, T_WALL, P_REF, GRIDS, ROOT, geometry


def triangulate(xy, ids):
    """在同一速度节点上作约束剖分，保留边界段和原面对应。"""
    sys.path.insert(0, str(ROOT / "tools/runtime/gmsh/4.15.2"))
    import gmsh
    boundary = np.r_[ids[0], ids[1:, -1], ids[-1, -2::-1], ids[-2:0:-1, 0]]
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("tunnel-air")
        geo = gmsh.model.geo
        for i, (x, y) in enumerate(xy):
            geo.addPoint(x, y, 0, 100, i+1)
        lines = [geo.addLine(int(a+1), int(b+1)) for a,b in zip(boundary, np.roll(boundary,-1))]
        surface = geo.addPlaneSurface([geo.addCurveLoop(lines)])
        geo.synchronize()
        gmsh.model.mesh.embed(0, (np.setdiff1d(np.arange(len(xy)), boundary)+1).tolist(), 2, surface)
        for line in lines:
            gmsh.model.mesh.setTransfiniteCurve(line, 2)
        gmsh.option.setNumber("Mesh.Algorithm", 3)  # 仅作约束 Delaunay 初始网。
        gmsh.option.setNumber("Mesh.Smoothing", 0)
        gmsh.model.mesh.generate(2)
        tags, points, _ = gmsh.model.mesh.getNodes()
        delta, original = cKDTree(xy).query(points.reshape(-1,3)[:,:2])
        assert len(tags)==len(xy) and max(delta)<1e-11 and len(np.unique(original))==len(xy)
        lookup = np.full(int(max(tags))+1,-1,int)
        lookup[tags.astype(int)] = original
        tr = lookup[gmsh.model.mesh.getElementsByType(2)[1].reshape(-1,3)]
    finally:
        gmsh.finalize()
    edges, count = np.unique(np.sort(np.vstack((tr[:,[0,1]], tr[:,[1,2]], tr[:,[2,0]])),axis=1),axis=0,return_counts=True)
    expected = set(map(tuple,np.sort(np.column_stack((boundary,np.roll(boundary,-1))),axis=1)))
    assert np.all(count<=2) and set(map(tuple,edges[count==1]))==expected
    return tr


def dirichlet(matrix, rhs, fixed, values):
    free=np.ones(len(rhs),bool)
    free[fixed]=False
    u=np.zeros(len(rhs))
    u[fixed]=values
    u[free]=spsolve(matrix[free][:,free],(rhs-matrix@u)[free])
    residual=(matrix@u-rhs)[free]
    scale=max(np.linalg.norm((matrix[free][:,free]@u[free]),np.inf),1e-30)
    assert np.max(abs(residual))/scale<1e-8
    return u,float(np.max(abs(residual))/scale)


def transport(xy,tr,wet,axis,far,target,pe):
    """轴对称 r 权重；外域入流强给温，出流自然扩散，绝热轴不加通量。"""
    v=xy[tr]
    edge=np.roll(v,-1,axis=1)-np.roll(v,1,axis=1)
    twice=np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])
    assert np.all(twice>0)
    gradient=np.stack((edge[:,:,1],-edge[:,:,0]),axis=2)/twice[:,None,None]
    area=twice/2
    dot=np.einsum('tid,tjd->tij',gradient,gradient)
    ii=np.broadcast_to(tr[:,:,None],dot.shape).ravel()
    jj=np.broadcast_to(tr[:,None,:],dot.shape).ravel()
    def assemble(element):
        return coo_matrix((element.ravel(),(ii,jj)),shape=(len(xy),len(xy))).tocsr()
    # 轴线和实体是同一条不可穿透流线。对速度作最小二乘流函数投影。
    psi_rhs=area[:,None]*np.einsum('tid,td->ti',gradient,np.column_stack((-target[:,1],target[:,0])))
    fixed=np.unique(np.r_[wet,axis])
    psi,psi_res=dirichlet(assemble(dot*area[:,None,None]),np.bincount(tr.ravel(),weights=psi_rhs.ravel()),fixed,0.)
    gradpsi=np.einsum('ti,tid->td',psi[tr],gradient)
    stream=np.column_stack((gradpsi[:,1],-gradpsi[:,0]))
    projection=np.sqrt(np.sum(area*np.sum((target-stream)**2,axis=1))/np.sum(area*np.sum(target**2,axis=1)))
    diffusion=dot*(area*v[:,:,1].mean(axis=1)/pe)[:,None,None]
    advection=(area[:,None]*np.einsum('td,tjd->tj',stream,gradient)/3)[:,None,:]
    raw=assemble(diffusion+advection)
    # 守恒低阶代数扩散：成对反向通量，零行和，非对角项非正。温度值不裁剪。
    up=raw.maximum(raw.T).maximum(0.)
    up.setdiag(0);up.eliminate_zeros()
    system=raw-up+diags(np.asarray(up.sum(axis=1)).ravel())
    outward=np.diff(psi[far])  # 此处 far 的点序从尾侧到鼻侧，单位外流为 -X。
    # 外域半圆从左到右为顺时针，外向质量流量为 -dψ。
    outward=-outward
    inlet=np.unique(np.column_stack((far[:-1],far[1:]))[outward<0])
    bound=np.r_[wet,inlet]
    theta,res=dirichlet(system,np.zeros(len(xy)),bound,np.r_[np.ones(len(wet)),np.zeros(len(inlet))])
    row=system@theta
    wall=row[wet]
    adv=float(np.dot(outward,(theta[far[:-1]]+theta[far[1:]])/2))
    net=adv-row[inlet].sum()
    balance=abs(net/wall.sum()-1)
    assert np.isfinite(theta).all() and theta.min()>=-1e-12 and theta.max()<=1+1e-12
    assert balance<.005 and wall.sum()>0
    return theta,wall,dict(linearResidual=res,streamfunctionResidual=psi_res,
            weightedVelocityProjectionError=float(projection),temperatureMin=float(theta.min()),
            temperatureMax=float(theta.max()),wallRelease=float(wall.sum()),
            farTotal=float(net),balance=float(balance))


def transfer(tri,wall,knots,diameter):
    """对原三角面内线性母线坐标作面积积分，底面包含轴心至边缘的全部子段。"""
    base_end=int(np.where(np.diff(wall[:,0])>0)[0][0])
    x=tri[:,:,0]/diameter
    row=np.interp(x,wall[base_end:,0],np.arange(base_end,len(wall)))
    bottom=(x==0).all(axis=1)
    row[bottom]=np.hypot(tri[bottom,:,1],tri[bottom,:,2])/(diameter*wall[base_end,1])*base_end
    row=np.rint(row).astype(int)
    lo=row.min(axis=1);hi=row.max(axis=1)
    assert np.all(hi>lo)
    double_high=np.sum(row==hi[:,None],axis=1)==2
    rr=[];cc=[];weights=[]
    bounds=np.unique(np.column_stack((lo,hi)),axis=0)
    for left,right in bounds:
        faces=np.flatnonzero((lo==left)&(hi==right))
        cuts=np.flatnonzero((knots>=left-1e-8)&(knots<=right+1e-8))
        assert len(cuts)>=2 and abs(knots[cuts[0]]-left)<1e-8 and abs(knots[cuts[-1]]-right)<1e-8
        for k in cuts[:-1]:
            width=knots[k+1]-knots[k]
            for q in (.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)):
                t=(knots[k]+q*width-left)/(right-left)
                weight=np.where(double_high[faces],t,1-t)*width/(right-left)
                for offset,factor in ((0,1-q),(1,q)):
                    rr.append(faces);cc.append(np.full(len(faces),k+offset));weights.append(weight*factor)
    matrix=coo_matrix((np.concatenate(weights),(np.concatenate(rr),np.concatenate(cc))),
                      shape=(len(tri),len(knots))).tocsr()
    assert np.max(abs(np.asarray(matrix.sum(axis=1))-1))<1e-12
    return matrix


def evaluate(folder,out,spacing=None):
    folder=folder.resolve()
    info=json.loads((folder/"condition.json").read_bytes())
    rec=json.loads((folder/"engine_process.json").read_bytes())
    assert rec["status"]=="terminal" and rec["exitCode"]==0 and not(rec["cancelled"] or rec["timedOut"])
    assert rec["runtimeFlavor"]=="official" and rec["executableVersion"]=="8.6.0" and all(rec["inputsStable"].values())
    for name in ("restart.csv","history.csv"):
        assert sha(folder/name)==rec["outputSHA256"][name]
    for name in ("mesh.su2","case.cfg","condition.json"):
        assert sha(folder/name)==rec["inputSHA256"][str(folder/name)]
    definition,corners,d,wall=geometry()
    assert definition["fileSHA256"]==info["surfaceSHA256"] and d==info["diameterM"]
    level,layers,extent=GRIDS[info["grid"]]
    w=np.array(info["wall"]);xy,_,_=volume(w,layers,extent,level)
    ids=np.arange(len(xy)).reshape(layers+1,-1)
    from simpleLauncherHeatMesh import heat_grid, sampler
    knots=np.arange(len(w))/level
    pmat=transfer(corners,wall,knots,d)
    if spacing is None:
        points=xy;tr=triangulate(xy,ids);wet=ids[0];axis=np.r_[ids[:,0],ids[:,-1]];far=ids[-1]
        def target_flow(u): return (u[tr]*points[tr,1,None]).mean(axis=1)
        def near_field(t): return t[ids[1]]
        matrix=pmat
    else:
        points,tr,wet,far,axis,wet_stations=heat_grid(xy,ids,spacing)
        quads=np.stack((ids[:-1,:-1],ids[:-1,1:],ids[1:,1:],ids[1:,:-1]),axis=-1).reshape(-1,4)
        flow_tri=np.vstack((quads[:,[0,1,2]],quads[:,[0,2,3]]))
        centroid=points[tr].mean(axis=1)
        interpolate=sampler(xy,flow_tri,centroid)
        def target_flow(u): return interpolate(u)*centroid[:,1,None]
        near_field=sampler(points,tr,xy[ids[1]])
        matrix=transfer(corners,wall,wet_stations/level,d)
    with (folder/"history.csv").open() as f:
        reader=csv.DictReader(f,skipinitialspace=True)
        reader.fieldnames=[s.strip() for s in reader.fieldnames]
        history=list(reader)
    keys=("rms[P]","rms[U]","rms[V]","rms[h]")
    maxima={key:max(float(r[key]) for r in history[-20:]) for key in keys}
    assert max(maxima.values())< -8
    final=int(history[-1]["Inner_Iter"])
    assert sha(folder/f"restart_{final:06d}.csv")==sha(folder/"restart.csv")
    prior=max(int(f.stem.split('_')[-1]) for f in folder.glob('restart_*.csv') if int(f.stem.split('_')[-1])<=final-20)
    power_unit=2*np.pi*RHO*SPEED*CP*(T_WALL-T_AIR)*d**2
    saved=[]
    for step in (prior,final):
        state=np.loadtxt(folder/f"restart_{step:06d}.csv",delimiter=',',skiprows=1)
        assert np.array_equal(state[:,0],np.arange(len(xy))) and np.all(np.isfinite(state))
        np.testing.assert_allclose(state[:,1:3],xy,atol=1e-12,rtol=0)
        assert abs(state[ids[0],4:6]).max()<1e-9 and np.min(P_REF+RHO*SPEED**2*state[:,3])>0
        theta,release,thermal=transport(points,tr,wet,axis,far,target_flow(state[:,4:6]),info["pecletD"])
        saved.append(release.sum()*power_unit)
    change=abs(saved[0]/saved[1]-1);assert change<.002
    area=.5*np.linalg.norm(np.cross(corners[:,1]-corners[:,0],corners[:,2]-corners[:,0]),axis=1)
    dual=np.asarray(matrix.T@area).ravel();assert np.all(dual>0)
    flux=matrix@(-release*power_unit/dual)
    pressure=P_REF+RHO*SPEED**2*(pmat@state[ids[0],3])
    near=T_AIR+(T_WALL-T_AIR)*(pmat@near_field(theta))
    power=float(-(area@flux));assert abs(power/saved[-1]-1)<1e-10
    assert np.all(np.isfinite(flux)) and np.all(flux<0)
    regions={}
    for region in definition['regions']:
        sl=slice(region['firstFace']-1,region['lastFace'])
        regions[region['name']]=dict(areaM2=float(area[sl].sum()),outwardW=float(-(area[sl]@flux[sl])),
                 minWallFluxWm2=float(flux[sl].min()),maxWallFluxWm2=float(flux[sl].max()))
    out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'temperature.npz',xy=points,triangles=tr,theta=theta,wall=wet,reaction=release)
    np.savetxt(out/'wall.csv',np.column_stack((np.arange(1,len(area)+1),area,pressure,flux,near)),
               delimiter=',',fmt=['%d']+['%.17g']*4,comments='',header='faceID,areaM2,pressurePa,heatFluxWm2,nearWallK')
    write_json(out/'summary.json',dict(grid=info['grid'],faces=len(area),areaM2=float(area.sum()),outwardW=power,
               regions=regions,thermal=thermal,heatSpacingD=spacing,heatNodes=len(points),residuals=maxima,heatStability=change,priorStep=prior,finalStep=final,
               sourceRestartSHA256=sha(folder/'restart.csv'),surfaceSHA256=definition['fileSHA256'],wallSHA256=sha(out/'wall.csv')))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('flow',type=Path);parser.add_argument('out',type=Path)
    parser.add_argument('--spacing',type=float)
    args=parser.parse_args();evaluate(args.flow,args.out,args.spacing)
