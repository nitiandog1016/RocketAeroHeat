"""独立核对 CM-011 新剖分的几何投影、重排不变性及错误输入拒绝。"""
import argparse
from pathlib import Path

import numpy as np

from compare import ROOT, binary_stl, read
from cm011Imported import match_surface
from cm011Tunnel import D, meridian, write_json
from runEmbeddedSU2 import sha


def direct_intersection(triangles, origin, ray):
    """对全部三角面求 Möller–Trumbore 交点，不使用生产代码的球面筛选。"""
    e1=triangles[:,1]-triangles[:,0]
    e2=triangles[:,2]-triangles[:,0]
    h=np.cross(ray,e2)
    det=np.einsum("ij,ij->i",e1,h)
    offset=origin-triangles[:,0]
    q=np.cross(offset,e1)
    keep=abs(det)>1e-22
    u=np.einsum("ij,ij->i",offset[keep],h[keep])/det[keep]
    v=(q[keep]@ray)/det[keep]
    t=np.einsum("ij,ij->i",e2[keep],q[keep])/det[keep]
    hits=t[(u>=-1e-11)&(v>=-1e-11)&(u+v<=1+1e-11)&(t>0)]
    assert len(hits) and np.ptp(hits)<1e-11
    return hits.mean()


def check(gui):
    stl=ROOT/"cases/CM011-alternate-SI.stl"
    tri=binary_stl(stl)
    vertices,ids=np.unique(tri.reshape(-1,3),axis=0,return_inverse=True)
    faces=ids.reshape(-1,3)+1
    mesh=Path((gui/"mesh_folder.txt").read_text("utf-8-sig").strip())
    base,reg=match_surface(vertices,faces)
    assert np.allclose(read(mesh/"body.json")["meridian"],base["meridian"],rtol=0,atol=1e-14)
    rng=np.random.RandomState(276)
    vp=rng.permutation(len(vertices))
    fp=rng.permutation(len(faces))
    shuffled=np.roll(np.argsort(vp)[faces[fp]-1]+1,1,axis=1)
    other,_=match_surface(vertices[vp],shuffled)
    np.testing.assert_array_equal(other["F"],shuffled)
    np.testing.assert_array_equal(other["V"],vertices[vp])
    np.testing.assert_array_equal(other["regionCode"],np.array(base["regionCode"])[fp])
    np.testing.assert_allclose(other["meridian"],base["meridian"],rtol=0,atol=1e-14)
    np.testing.assert_allclose(other["wallStation"],np.array(base["wallStation"])[vp],rtol=0,atol=2e-12)

    # 给定四个内部站。对每站的 48 个方位与原面直接求交，再检查平均位置。
    nominal=meridian(2)[0]*D
    origin=np.array([.5*(nominal[0,0]+nominal[-1,0]),0,0])
    gaps=[]
    for coarse_station in (3,14,29,41):
        j=2*coarse_station
        ray=nominal[j]-origin[:2]
        ray/=np.linalg.norm(ray)
        lengths=[direct_intersection(tri,origin,[ray[0],ray[1]*np.cos(phi),ray[1]*np.sin(phi)])
                 for phi in np.arange(48)*2*np.pi/48]
        xr=origin[:2]+ray*np.mean(lengths)
        gaps.append(float(np.linalg.norm(xr-D*np.asarray(base["meridian"])[j])))
    assert max(gaps)<1e-13

    # 保留名义外形上的顶点，构建闭合粗面，单独验证面片偏差拒绝路径。
    x,r=nominal[54]
    coarse_v=np.array([[nominal[0,0],0,0],[x,r,0],[x,0,r],[x,-r,0],[x,0,-r],[nominal[-1,0],0,0]])
    coarse_f=np.array([[1,2+(i+1)%4,2+i] for i in range(4)]+
                      [[2+i,2+(i+1)%4,6] for i in range(4)])
    damaged=vertices.copy()
    k=np.argmax(np.linalg.norm(vertices[:,1:],axis=1))
    damaged[k,1:]*=1.01
    reverse=faces.copy()
    reverse[0]=reverse[0,::-1]
    duplicate=faces.copy()
    duplicate[1]=duplicate[0]
    bad=(("hole",vertices,faces[1:],"开口"),("duplicate",vertices,duplicate,"非流形"),
         ("oneReversed",vertices,reverse,"方向"),("allReversed",vertices,faces[:,::-1],"内向"),
         ("wrongScale",vertices*.001,faces,"名义轮廓"),("localBump",damaged,faces,"名义轮廓"),
         ("coarseShell",coarse_v,coarse_f,"边中点"))
    rejected={}
    for name,v,f,reason in bad:
        try:
            match_surface(v,f)
        except ValueError as err:
            assert reason in str(err), str(err)
            rejected[name]=str(err)
        else:
            raise AssertionError("应当拒绝的试件："+name)
    return dict(stlSHA256=sha(stl),vertices=len(vertices),faces=len(faces),registration=reg,
                shuffledAndCyclic=True,independentRays=192,meanPositionDifferenceM=gaps,rejections=rejected)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gui",type=Path)
    parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    write_json(args.new_output,check(args.gui))
