"""在固定 0.05 m/s 风洞热网传递三分量速度，按共享轴点和开放远界投影质量通量。"""
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib.tri as mtri
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

DOC=Path(__file__).resolve().parent
ROOT=DOC.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from simpleLauncherAngularHeat import ANGLE,SECTORS,domain
from simpleLauncherHeatMesh import sampler
from runEmbeddedSU2 import sha


def read(p):
    return json.loads(p.read_bytes())


def boundary_sampler(points,triangles,targets):
    """流体内沿用原三角形插值，落在真实边界段上的舍入点按该段两端取值。"""
    miss=mtri.Triangulation(*points.T,triangles).get_trifinder()(*targets.T)<0
    interior=sampler(points,triangles,targets[~miss])
    all_edges=np.sort(np.vstack((triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]])),axis=1)
    unique,counts=np.unique(all_edges,axis=0,return_counts=True)
    border=unique[counts==1];a,b=points[border].transpose(1,0,2)
    delta=b-a;squared=np.sum(delta**2,axis=1);ends=[];weights=[]
    for target in targets[miss]:
        t=np.clip(np.sum((target-a)*delta,axis=1)/squared,0.,1.)
        distance=np.linalg.norm(target-a-t[:,None]*delta,axis=1)
        index=int(distance.argmin())
        assert distance[index]<1e-11
        ends.append(border[index]);weights.append([1-t[index],t[index]])
    ends=np.asarray(ends,int);weights=np.asarray(weights)
    def sample(field):
        out=np.empty((len(targets),)+field.shape[1:]);out[~miss]=interior(field)
        if miss.any():out[miss]=np.einsum("ij,ij...->i...",weights,field[ends])
        return out
    return sample


def build_model(flow_job,prepared):
    ref=read(DOC.parent/"v352/baseline.json")["coarse"]
    m=domain(flow_job,"coarse",ref)
    n=len(m["xy"]);xy=m["xy"];tr=m["tri"];p=xy[tr]
    area=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])/2
    ctr=p.mean(axis=1);radius=p[:,:,1];m["cellArea"]=area;m["centroid"]=ctr
    # ∫ r N_i dA 的精确线性三角形积分用作每弧度的集总流体热容量。
    m["mass"]=np.bincount(tr.ravel(),weights=(area[:,None]*(radius.sum(axis=1)[:,None]+radius)/12).ravel(),minlength=n)
    np.testing.assert_allclose(m["mass"].sum(),m["volumePerRadian"],rtol=1e-13)
    dual_area=np.bincount(tr.ravel(),weights=np.repeat(area/3,3),minlength=n)
    m["angularArea"]=dual_area/ANGLE
    edge=np.roll(p,-1,axis=1)-np.roll(p,1,axis=1)
    grad=np.stack((edge[:,:,1],-edge[:,:,0]),axis=2)/(2*area[:,None,None])
    next_tr=np.roll(tr,-1,axis=1);pairs=np.stack((tr,next_tr),axis=2)
    edges,indices=np.unique(np.sort(pairs,axis=2).reshape(-1,2),axis=0,return_inverse=True)
    m["edges"]=edges;m["edgeIndices"]=indices
    signs=np.sign(next_tr-tr)
    m["edgeVectors"]=(np.roll(grad,-1,axis=1)-grad)*(signs*area[:,None]*ctr[:,1,None]/3)[:,:,None]
    mid=(p+np.roll(p,-1,axis=1))/2
    lengths=np.linalg.norm(np.roll(p,-1,axis=1)-p,axis=2)
    w=np.linalg.norm(mid-ctr[:,None,:],axis=2)/lengths*(mid[:,:,1]+ctr[:,1,None])/2
    m["weights"]=np.bincount(indices,weights=w.ravel())
    i,j=edges.T;ix=np.arange(len(i))
    B=coo_matrix((np.r_[np.ones(len(i)),-np.ones(len(i))],(np.r_[i,j],np.r_[ix,ix])),shape=(n,len(i))).tocsr()
    m["incidence"]=B
    far=m["far"];e=np.diff(xy[far],axis=0);far_mid=(xy[far[:-1]]+xy[far[1:]])/2
    # 远界子午线由尾侧向鼻侧沿顺时针方向排列，外法向为左法向。
    normal=np.column_stack((-e[:,1],e[:,0]))
    m["farNormalArea"]=normal*far_mid[:,1,None]
    far_weights=np.linalg.norm(normal,axis=1)*far_mid[:,1]/np.linalg.norm(e,axis=1)
    open_w=np.zeros(n);np.add.at(open_w,far[:-1],far_weights*.5);np.add.at(open_w,far[1:],far_weights*.5)
    m["openWeight"]=open_w
    m["projectionMatrix"]=(B@diags(m["weights"])@B.T+diags(open_w)).tocsr()
    m["angularWeight"]=m["conductance"]*m["pe"]
    base_col=np.asarray(m["matrix"].sum(axis=0)).ravel()
    base_b=np.zeros(n);np.add.at(base_b,far[:-1],m["outward"]*.5);np.add.at(base_b,far[1:],m["outward"]*.5)
    assert max(abs(base_col-base_b))<1e-10
    m["baseBoundary"]=base_b
    m["upper"]=np.asarray(m["matrix"][i,j]).ravel();m["lower"]=np.asarray(m["matrix"][j,i]).ravel()
    assert max(m["upper"].max(),m["lower"].max())<1e-12
    state=np.loadtxt(flow_job/"coarse/restart.csv",skiprows=1,delimiter=",")
    info=read(flow_job/"coarse/condition.json")
    ids=np.arange(len(state)).reshape(info["layers"]+1,len(info["wall"]))
    quads=np.stack((ids[:-1,:-1],ids[:-1,1:],ids[1:,1:],ids[1:,:-1]),axis=-1).reshape(-1,4)
    flow_tri=np.vstack((quads[:,[0,1,2]],quads[:,[0,2,3]]))
    m["atCentroid"]=sampler(state[:,1:3],flow_tri,ctr)
    m["atNode"]=boundary_sampler(state[:,1:3],flow_tri,xy)
    m["atFar"]=boundary_sampler(state[:,1:3],flow_tri,far_mid)
    on_axis=state[:,2]==0
    lookup=np.empty((SECTORS,len(state)),int);lookup[:,on_axis]=np.arange(on_axis.sum())
    lookup[:,~on_axis]=on_axis.sum()+np.arange(SECTORS*(~on_axis).sum()).reshape(SECTORS,-1)
    m["flowIds"]=lookup;m["flowWallCount"]=len(info["wall"])
    with np.load(prepared/"stations.npz") as native:
        np.testing.assert_array_equal(lookup[:,:len(info["wall"])],native["wall"])
        coord=native["xyz"][lookup]
        np.testing.assert_allclose(coord[:,:,0],np.broadcast_to(state[:,1],lookup.shape),rtol=0,atol=1e-12)
        np.testing.assert_allclose(np.hypot(coord[:,:,1],coord[:,:,2]),np.broadcast_to(state[:,2],lookup.shape),rtol=0,atol=1e-12)
    m["seedVelocity"]=np.loadtxt(prepared/"control.csv",delimiter=",",skiprows=1,usecols=(5,6,7))
    return m


def project(m,native):
    """对相对已守恒基流的通量增量做加权最小修正；远场可交换质量，实体壁始终不可穿透。"""
    ids=m["flowIds"];i,j=m["edges"].T
    delta=native[:,5:8]-m["seedVelocity"]
    v=delta[ids];phase=np.arange(SECTORS)*ANGLE
    c,s=np.cos(phase)[:,None],np.sin(phase)[:,None]
    vr=v[:,:,1]*c+v[:,:,2]*s;vt=-v[:,:,1]*s+v[:,:,2]*c
    uv=np.stack((v[:,:,0],vr),axis=2);uv=.5*(uv+np.roll(uv,-1,axis=0))
    edge=np.empty((SECTORS,len(i)));angular=np.empty((SECTORS,len(m["xy"])))
    boundary=np.zeros_like(angular)
    for k in range(SECTORS):
        cell=m["atCentroid"](uv[k])
        per_tri=np.einsum('tjd,td->tj',m["edgeVectors"],cell)
        edge[k]=np.bincount(m["edgeIndices"],weights=per_tri.ravel(),minlength=len(i))
        angular[k]=m["atNode"](vt[(k+1)%SECTORS])*m["angularArea"]
        angular[k,m["axis"]]=0.
        fb=np.einsum('ij,ij->i',m["atFar"](uv[k]),m["farNormalArea"])*.5
        np.add.at(boundary[k],m["far"][:-1],fb);np.add.at(boundary[k],m["far"][1:],fb)
    B=m["incidence"]
    def divergence():
        return (B@edge.T).T+angular-np.roll(angular,1,axis=0)+boundary
    before=divergence()
    source=np.fft.rfft(before,axis=0);coeff=np.zeros_like(source)
    off=np.ones(len(m["xy"]),bool);off[m["axis"]]=False
    for mode,rhs in enumerate(source):
        if not np.any(rhs):continue
        matrix=m["projectionMatrix"]+diags(4*np.sin(mode*ANGLE/2)**2*m["angularWeight"])
        free=off if mode else np.ones_like(off)
        parts=spsolve(matrix[free][:,free],np.column_stack((rhs.real[free],rhs.imag[free])))
        coeff[mode,free]=parts[:,0]+1j*parts[:,1]
    phi=np.fft.irfft(coeff,n=SECTORS,axis=0)
    base_flux=m["upper"]-m["lower"]
    full_norm=increment_norm=correction_norm=0.
    for k in range(SECTORS):
        for old,base,correction in (
            (edge[k],base_flux,(phi[k,i]-phi[k,j])*m["weights"]),
            (angular[k],0.,(phi[k]-phi[(k+1)%SECTORS])*m["angularWeight"]),
            (boundary[k],m["baseBoundary"],phi[k]*m["openWeight"])):
            full_norm+=np.sum((old+base)**2);increment_norm+=np.sum(old**2);correction_norm+=np.sum(correction**2)
            old-=correction
    after=divergence()
    def max_unique(a):
        return float(max(np.abs(a[:,off]).max(),np.abs(a[:,~off].sum(axis=0)).max()))
    scale=max(np.abs(base_flux).max(),np.abs(m["baseBoundary"]).max())
    measure=dict(relativeProjection=float(np.sqrt(correction_norm/full_norm)),
        projectionRelativeToIncrement=float(np.sqrt(correction_norm/increment_norm)) if increment_norm else 0.,
        initialMaxDivergence=max_unique(before),projectedMaxDivergence=max_unique(after),
        normalizedDivergence=max_unique(after)/scale,maxAngularMassFlux=float(np.abs(angular).max()),
        netFarFlux=float(np.sum(boundary+m["baseBoundary"])))
    assert measure["relativeProjection"]<.05 and measure["normalizedDivergence"]<1e-8,measure
    # 守恒中心通量的增量成对装配，再成对加入仅够保持非正非对角项的代数扩散。
    matrices=[]
    for k in range(SECTORS):
        half=.5*edge[k]
        viscosity=np.maximum(np.maximum(m["upper"]+half,m["lower"]-half),0.)
        a=half-viscosity;b=-half-viscosity
        diag=-np.bincount(i,weights=a,minlength=len(off))-np.bincount(j,weights=b,minlength=len(off))
        change=coo_matrix((np.r_[a,b],(np.r_[i,j],np.r_[j,i])),shape=m["matrix"].shape).tocsr()
        matrices.append(m["matrix"]+change+diags(diag))
    half=.5*angular
    loss=np.maximum(abs(half)-m["conductance"],0.)
    right=-m["conductance"]+half-loss
    left=np.roll(-m["conductance"]-half-loss,1,axis=0)
    assert max(right.max(),left.max())<1e-12
    return matrices,right,left,boundary+m["baseBoundary"],measure
