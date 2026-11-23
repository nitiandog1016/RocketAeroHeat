"""通过各分区共用的母线位置与方位参数，将原面和求解壁片按交叠配对。"""
import numpy as np
from shortSupportMesh import COUNTS, PHASE, REGIONS, SECTORS


def clip(poly,axis,bound,side):
    result=[]
    for i,end in enumerate(poly):
        start=poly[i-1];a=side*(start[axis]-bound);b=side*(end[axis]-bound)
        if (a<0)!=(b<0):
            w=a/(a-b)
            result.append((start[0]+w*(end[0]-start[0]),start[1]+w*(end[1]-start[1])))
        if b>=0:result.append(end)
    return result


def polygon_area(poly):
    x,y=poly[0]
    return .5*abs(sum((a-x)*(d-y)-(b-y)*(c-x) for (a,b),(c,d) in zip(poly[1:],poly[2:])))


def sphere_sectors(v,data,patch):
    """中心投影把同一扇区中的原面直边映为直线，保留接头极窄片的方向。"""
    half=np.tan(np.pi/SECTORS);sectors=[]
    for k,strip in enumerate(patch.reshape(SECTORS,sum(COUNTS),4,3)):
        angle=PHASE+(k+.5)*2*np.pi/SECTORS
        radial=np.array([np.cos(angle),0,-np.sin(angle)])
        tangent=np.array([-np.sin(angle),0,-np.cos(angle)])
        band=strip[23:71]
        row=np.vstack((band[:,0],band[-1,3]))
        levels=row[:,1]/(row@radial)
        measure=2*half*np.diff(levels)
        for name,edge in (("joinNegative",0),("joinPositive",-1)):
            seam=v[np.asarray(data[name])[:,0]-1]
            depth=seam@radial
            use=(depth>0)&(abs(seam@tangent)<=half*depth+1e-12)
            uv=np.c_[seam[use]@tangent,seam[use,1]]/depth[use,None]
            uv=uv[np.argsort(uv[:,0])]
            np.testing.assert_allclose(uv[[0,-1],0],[-half,half],rtol=0,atol=1e-10)
            offset=np.trapz(uv[:,1]-levels[edge],uv[:,0])
            measure[edge]+=(-offset if edge==0 else offset)
        assert min(measure)>0
        sectors.append((radial,tangent,levels,measure))
    return sectors


def split_sphere_sectors(v,data,patch):
    """两个子壁片共用父扇区的中心投影图，在图中沿 u=0 精确切开接缝。"""
    halves=patch.reshape(SECTORS,2,sum(COUNTS),4,3)
    parent=np.stack((halves[:,0,:,0],halves[:,1,:,1],halves[:,1,:,2],halves[:,0,:,3]),axis=2)
    parent_charts=sphere_sectors(v,data,parent)
    half=np.tan(np.pi/SECTORS);charts=[]
    for radial,tangent,levels,whole in parent_charts:
        measures=np.tile(half*np.diff(levels),(2,1))
        for name,edge,sign in (("joinNegative",0,-1),("joinPositive",-1,1)):
            seam=v[np.array(data[name])[:,0]-1];depth=seam@radial
            keep=(depth>0)&(abs(seam@tangent)<=half*depth+1e-12)
            uv=np.c_[seam[keep]@tangent,seam[keep,1]]/depth[keep,None];uv=uv[np.argsort(uv[:,0])]
            for j,(lo,hi) in enumerate(((-half,0),(0,half))):
                x=np.r_[lo,uv[(uv[:,0]>lo)&(uv[:,0]<hi),0],hi]
                measures[j,edge]+=sign*np.trapz(np.interp(x,uv[:,0],uv[:,1])-levels[edge],x)
        np.testing.assert_allclose(measures.sum(0),whole,atol=1e-12,rtol=1e-12)
        assert measures.min()>0
        charts.extend((radial,tangent,levels,measure) for measure in measures)
    return charts


def sphere_pairs(tri,normal,angular,sectors,patch,subdivisions=1):
    half=np.tan(np.pi/SECTORS)
    for k in range(int(np.floor(min(angular))),int(np.ceil(max(angular)))):
        index=k%(SECTORS*subdivisions)
        radial,tangent,levels,measure=sectors[index]
        depth=tri@radial
        assert min(depth)>0
        uv=np.c_[tri@tangent,tri[:,1]]/depth[:,None]
        lo_u=0 if subdivisions==2 and k%2 else -half
        hi_u=0 if subdivisions==2 and k%2==0 else half
        slab=clip(clip(uv.tolist(),0,lo_u,1),0,hi_u,-1)
        if len(slab)<3:continue
        heights=[p[1] for p in slab]
        lo=max(0,np.searchsorted(levels,min(heights))-1)
        hi=min(COUNTS[2],np.searchsorted(levels,max(heights))+1)
        for l in range(lo,hi):
            piece=slab
            if l>0:piece=clip(piece,1,levels[l],1)
            if l<COUNTS[2]-1:piece=clip(piece,1,levels[l+1],-1)
            if len(piece)<3:continue
            projected=polygon_area(piece)
            if projected<=1e-19:continue
            poly=np.asarray(piece)
            direction=radial+poly[:,0,None]*tangent
            direction[:,1]=poly[:,1]
            original=direction*((tri[0]@normal)/(direction@normal))[:,None]
            rel=original-original[0]
            fraction=np.linalg.norm(np.cross(rel,np.roll(rel,-1,axis=0)).sum(0))/np.linalg.norm(normal)
            j=index*sum(COUNTS)+23+l
            pt=original.mean(0)
            ray=pt/(pt@radial);ray[1]=np.clip(ray[1],levels[l],levels[l+1])
            face=patch[j];n=np.cross(face[2]-face[0],face[3]-face[1])
            native=ray*((face[0]@n)/(ray@n))
            yield j,fraction,projected/measure[l],pt,native


def transfer(folder,data,condition,stats):
    v=np.array(data["V"])/condition["diameterM"];f=np.array(data["F"])-1
    region=np.array(data["regionID"]);area=np.array(data["areaM2"])
    with np.load(folder/"mapping.npz") as z:
        xyz=z["xyz"];corners=z["corners"];source_area=z["nativeArea"];groups=z["group"]
    patch=xyz[corners];src_vec=np.cross(patch[:,2]-patch[:,0],patch[:,3]-patch[:,1])
    tri=v[f];dst_vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    phi=(np.arctan2(-v[:,2],v[:,0])-PHASE)%(2*np.pi)
    count=stats.get("angularSectors",SECTORS)
    assert count in (SECTORS,2*SECTORS)
    subdivisions=count//SECTORS
    angular=phi*count/(2*np.pi)
    lat=np.arctan2(v[:,1],np.hypot(v[:,0],v[:,2]))
    joints=[]
    max_joint_error=0.
    for name,j in (("joinNegative",23),("joinPositive",71)):
        loop=np.array(data[name])[:,0]-1
        native_ring=patch.reshape(count,-1,4,3)[:,j,0]
        # 每个原接缝节点落在本方位对应的原生接缝线段附近。
        a=native_ring[np.floor(angular[loop]).astype(int)%count]
        b=native_ring[(np.floor(angular[loop]).astype(int)+1)%count]
        e=b-a;w=np.clip(((v[loop]-a)*e).sum(1)/(e*e).sum(1),0,1)
        max_joint_error=max(max_joint_error,np.linalg.norm(v[loop]-a-w[:,None]*e,axis=1).max())
        joints.append((np.interp(phi,phi[loop],v[loop,1],period=2*np.pi),
                       np.interp(phi,phi[loop],lat[loop],period=2*np.pi)))
    assert max_joint_error<.01,max_joint_error
    height={3:(v[:,1]+1.5)/(joints[0][0]+1.5),2:(v[:,1]-joints[1][0])/(1.5-joints[1][0]),
            1:(lat-joints[0][1])/(joints[1][1]-joints[0][1])}
    offsets=dict(zip(REGIONS,np.r_[0,np.cumsum(COUNTS)[:-1]]))
    counts=dict(zip(REGIONS,COUNTS));stride=sum(COUNTS)
    target=[];source=[];share_t=[];share_s=[]
    max_distance=max_joint_error;min_cos=1.
    sectors=sphere_sectors(v,data,patch) if subdivisions==1 else split_sphere_sectors(v,data,patch)
    def pair(i,j,wt,ws,point,position):
        nonlocal max_distance,min_cos
        max_distance=max(max_distance,np.linalg.norm(point-position))
        cosine=dst_vec[i]@src_vec[j]/(np.linalg.norm(dst_vec[i])*np.linalg.norm(src_vec[j]))
        min_cos=min(min_cos,cosine)
        target.append(i);source.append(j);share_t.append(wt);share_s.append(ws)
    for i,(nodes,g) in enumerate(zip(f,region)):
        start=offsets[g];n=counts[g]
        if g in (4,5):
            # 两套端盖共用同一个 64 边形，原三角扇形由该角区的三排原生壁片完整填充。
            centroid=tri[i].mean(axis=0)
            k=int(((np.arctan2(-centroid[2],centroid[0])-PHASE)%(2*np.pi))*SECTORS/(2*np.pi))
            for strip in range(k*subdivisions,(k+1)*subdivisions):
                for j in range(strip*stride+start,strip*stride+start+n):
                    point=patch[j].mean(axis=0)
                    assert abs(point[1]-centroid[1])<1e-13
                    pair(i,j,source_area[j]/area[i],1.,point,point)
            continue
        a=angular[nodes];a=a[0]+(a-a[0]+count/2)%count-count/2
        if g==1:
            for j,wt,ws,pt,native in sphere_pairs(tri[i],dst_vec[i],a,sectors,patch,subdivisions):
                pair(i,j,wt,ws,pt,native)
            continue
        b=height[g][nodes]*n
        chart=np.c_[a,b]
        near=np.abs(chart-np.rint(chart))<1e-11
        chart[near]=np.rint(chart[near])
        poly=chart.tolist();total=polygon_area(poly);assert total>0,(i,g)
        basis=(chart[1:]-chart[0]).T
        lo=np.floor(chart.min(0)).astype(int);hi=np.ceil(chart.max(0)).astype(int)
        assert lo[1]>=0 and hi[1]<=n
        for k in range(lo[0],hi[0]):
            slab=clip(clip(poly,0,k,1),0,k+1,-1)
            for l in range(lo[1],hi[1]):
                piece=clip(clip(slab,1,l,1),1,l+1,-1)
                if len(piece)<3:continue
                size=polygon_area(piece)
                if size<=1e-18:continue
                mid=np.mean(piece,axis=0);u,w=mid-[k,l]
                bary=np.linalg.solve(basis,mid-chart[0])
                point=tri[i,0]+bary@(tri[i,1:]-tri[i,0])
                j=(k%count)*stride+start+l
                position=np.array([(1-u)*(1-w),u*(1-w),u*w,(1-u)*w])@patch[j]
                pair(i,j,size/total,size,point,position)
    target=np.asarray(target);source=np.asarray(source)
    st=np.bincount(target,share_t,minlength=len(f));ss=np.bincount(source,share_s,minlength=len(corners))
    err_t=max(abs(st-1));err_s=max(abs(ss-1))
    group_area=np.bincount(region,area)[1:];relative=np.array(stats["nativeAreasM2"])/group_area-1
    assert (region[target]==groups[source]).all() and err_t<1e-7 and err_s<1e-7,(err_t,err_s)
    assert min_cos>0 and max_distance<.012 and max(abs(relative))<.01,(min_cos,max_distance,relative)
    np.savez_compressed(folder/"transfer.npz",target=target,source=source,targetFraction=share_t,
                        sourceFraction=share_s,targetArea=area,region=region)
    return dict(rows=len(target),targetFaces=len(f),sourceFaces=len(corners),targetCoverageError=err_t,
                sourceCoverageError=err_s,minNormalCos=min_cos,maxSampleDistanceD=max_distance,
                maxJoinDistanceD=max_joint_error,targetAreasM2=group_area.tolist(),areaRelativeDifference=relative.tolist())
