"""把每个加密条带的分片线性读数按面积积分回原三角面，平底另按半径积分。"""
from functools import lru_cache
import numpy as np
from scipy.sparse import coo_matrix


def cross2(a,b):
    return a[0]*b[1]-a[1]*b[0]


def clip(poly,a,b):
    """将含形函数值的凸多边形裁至直线左边，交点处线性分配两端权重。"""
    result=[]
    for v,w in zip(poly,poly[1:]+poly[:1]):
        dv=cross2(b-a,v[:2]-a);dw=cross2(b-a,w[:2]-a)
        if dv>=0:result.append(v)
        if (dv<0<dw) or (dw<0<dv):result.append(v+(w-v)*(dv/(dv-dw)))
    return result


@lru_cache(maxsize=2048)
def side_weights(r0,r2,corners):
    """同一扇片的两个半带各自三角剖分。裁到原面后对线性形函数精确求积。"""
    radius=np.repeat([r0,(r0+r2)/2,r2],2)
    uv=np.column_stack((np.repeat([0.,.5,1.],2),radius*np.tile([-1,1],3)))
    keep=uv[list(corners)]
    twice=cross2(keep[1]-keep[0],keep[2]-keep[0])
    if twice<0:keep=keep[::-1]
    area=abs(twice)/2
    assert area>0
    diag=next(((i,j) for i in corners for j in corners if i<2 and j>=4 and i%2!=j%2),None)
    rising=diag is None or diag==(0,5)
    total=np.zeros(6)
    for left in (0,2):
        elements=((0,2,3),(0,3,1)) if rising else ((0,2,1),(1,2,3))
        for indices in elements:
            ids=np.array(indices)+left;poly=list(np.column_stack((uv[ids],np.eye(3))))
            for i in range(3):
                poly=clip(poly,keep[i],keep[(i+1)%3])
            for k in range(1,len(poly)-1):
                a,b,c=poly[0],poly[k],poly[k+1]
                piece=abs(cross2(b[:2]-a[:2],c[:2]-a[:2]))/2
                total[ids]+=piece*(a[2:]+b[2:]+c[2:])/(3*area)
    assert total.min()>=-1e-13 and abs(total.sum()-1)<1e-12
    return total


def original_face_map(tri,d,old_wall,fine_wall,ids):
    """计算每个原面上的原生壁节点权重，坐标轴与 SU2 壁节点编号保持一致。"""
    sectors,width=ids.shape
    assert sectors==48 and width==2*len(old_wall)-1 and np.allclose(fine_wall[::2],old_wall)
    base_end=np.flatnonzero(np.diff(old_wall[:,0])>0)[0]
    mid=tri.mean(axis=1);angle=2*np.pi/sectors
    sector=np.floor(np.mod(np.arctan2(mid[:,2],mid[:,1]),2*np.pi)/angle).astype(int)
    base=(tri[:,:,0]==0).all(axis=1)
    row=np.interp(tri[:,:,0]/d,old_wall[base_end:,0],np.arange(base_end,len(old_wall)))
    ids_old=np.rint(row).astype(int)
    assert abs(row[~base]-ids_old[~base]).max()<1e-8
    index=[];nodes=[];weights=[]
    for face,j in enumerate(sector):
        k=(j+1)%sectors
        if base[face]:
            r=fine_wall[:2*base_end+1,1]/old_wall[base_end,1]
            a,b=r[:-1],r[1:];ray=np.zeros(len(r))
            ray[:-1]+=(b-a)*(2*a+b)/3;ray[1:]+=(b-a)*(a+2*b)/3
            ids_face=np.r_[ids[j,:len(r)],ids[k,:len(r)]];w=np.tile(ray/2,2)
        else:
            band=ids_old[face].min()
            assert ids_old[face].max()==band+1
            phi=(j+.5)*angle
            plus=(-tri[face,:,1]*np.sin(phi)+tri[face,:,2]*np.cos(phi))>0
            corners=tuple(4*(ids_old[face]-band)+plus.astype(int))
            w=side_weights(old_wall[band,1],old_wall[band+1,1],corners)
            a=2*band
            ids_face=ids[[j,k],a:a+3].T.ravel()
        use=w!=0
        nodes.extend(ids_face[use]);weights.extend(w[use]);index.extend([face]*int(use.sum()))
    matrix=coo_matrix((weights,(index,nodes)),shape=(len(tri),int(ids.max())+1)).tocsr()
    np.testing.assert_allclose(np.asarray(matrix.sum(axis=1)).ravel(),1.,atol=1e-12,rtol=0)
    return matrix
