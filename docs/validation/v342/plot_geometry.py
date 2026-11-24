"""从工作台的冻结 STL 绘制两个等比例正交视图，供几何入口的报告复核。"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def view(ax, triangles, direction, horizontal, title, bar):
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    right = np.array(horizontal) - np.dot(horizontal, direction)*direction
    right /= np.linalg.norm(right)
    up = np.cross(direction, right)
    coords = triangles @ np.array([right, up]).T
    normal = np.cross(triangles[:,1]-triangles[:,0], triangles[:,2]-triangles[:,0])
    normal /= np.linalg.norm(normal,axis=1)[:,None]
    light = .64 + .3*np.clip(normal @ (.4*up+.8*direction-.2*right),-1,1)
    color = light[:,None]*[.68,.78,.85]
    lo, hi = coords.min(axis=(0,1)), coords.max(axis=(0,1))
    # 用逐像素深度判定处理薄壁的内外表面，避免按三角形重心排序产生穿透伪影。
    scale = (1600 if bar==5 else 750)/max(hi-lo)
    screen=(coords-lo)*scale+.5
    width,height=np.ceil((hi-lo)*scale+1).astype(int)
    nearest=np.full(width*height,-np.inf)
    rgb=np.full((width*height,3),255,dtype=np.uint8)
    depths=triangles @ direction
    for offset in range(0,len(triangles),512):
        block=screen[offset:offset+512]
        xy0=np.ceil(block.min(axis=1)).astype(int)
        size=np.maximum(0,np.floor(block.max(axis=1)).astype(int)-xy0+1)
        counts=size.prod(axis=1)
        face=np.repeat(np.arange(len(block)),counts)
        if not len(face):continue
        idx=np.arange(len(face))-np.repeat(np.cumsum(counts)-counts,counts)
        x=xy0[face,0]+idx%size[face,0];y=xy0[face,1]+idx//size[face,0]
        a,b,c=block[face].transpose(1,0,2)
        ab=b-a;ac=c-a;px=x-a[:,0];py=y-a[:,1]
        den=ab[:,0]*ac[:,1]-ab[:,1]*ac[:,0]
        valid=np.abs(den)>1e-14
        u=np.divide(px*ac[:,1]-py*ac[:,0],den,out=np.zeros_like(px),where=valid)
        v=np.divide(ab[:,0]*py-ab[:,1]*px,den,out=np.zeros_like(px),where=valid)
        hit=valid & (u>=-1e-10) & (v>=-1e-10) & (u+v<=1+1e-10)
        face=face[hit]+offset;u=u[hit];v=v[hit];pixel=y[hit]*width+x[hit]
        if not len(pixel):continue
        z=(1-u-v)*depths[face,0]+u*depths[face,1]+v*depths[face,2]
        order=np.lexsort((z,pixel))
        front=order[np.r_[pixel[order][:-1]!=pixel[order][1:],True]]
        front=front[z[front]>nearest[pixel[front]]]
        nearest[pixel[front]]=z[front];rgb[pixel[front]]=(255*color[face[front]]).astype(np.uint8)
    ax.imshow(rgb.reshape(height,width,3),origin="lower",extent=[lo[0]-1/scale,lo[0]+(width-1)/scale,
              lo[1]-1/scale,lo[1]+(height-1)/scale],interpolation="hanning")
    pad = max(hi-lo)*.05
    ax.set(xlim=(lo[0]-pad,hi[0]+pad),ylim=(lo[1]-2.6*pad,hi[1]+pad),aspect="equal")
    ax.axis("off")
    ax.plot([lo[0],lo[0]+bar],[lo[1]-1.65*pad]*2,color="#253d4a",lw=2)
    ax.text(lo[0]+bar/2,lo[1]-1.45*pad,f"{bar} m",ha="center",color="#253d4a",fontsize=9)
    ax.set_title(title,loc="left",color="#253d4a",fontsize=12)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stl",type=Path);p.add_argument("png",type=Path);a=p.parse_args()
    spec=json.loads((ROOT/"cases/rfz_rec/geometry.json").read_bytes())
    assert hashlib.sha256(a.stl.read_bytes()).hexdigest()==spec["stlSHA256"]
    dt=[("n","<f4",3),("v","<f4",(3,3)),("tag","<u2")]
    triangles=np.fromfile(a.stl,dtype=dt,offset=84)["v"].astype(float)
    plt.rcParams["font.family"]="Microsoft YaHei"
    fig,(ax,bx)=plt.subplots(1,2,figsize=(11.8,3.8),gridspec_kw={"width_ratios":[2.45,1]},constrained_layout=True)
    view(ax,triangles,[-.35,-1,.52],[1.,0,0],"RFZ-REC  ·  CAD 边界",5)
    view(bx,triangles,[-1.,0,0],[0.,-1,0],"尾部正视",2)
    fig.savefig(a.png,dpi=155,facecolor="white")
