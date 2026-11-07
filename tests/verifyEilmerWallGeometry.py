"""补充独立读取原生边界块：逐三角面面积和外向关系，不重复单元/初场审核。"""
from pathlib import Path
import argparse
import hashlib
import json
import struct
import time
import numpy as np


def verify(case,out):
    start=time.monotonic();case=Path(case);dest=Path(out);assert not dest.exists()
    m=json.loads((case/'case_manifest.json').read_text());data=(case/'lmrsim/grid/grid-0000').read_bytes()
    assert data.startswith(b'unstructured_grid 1.0');pos=21
    # 固定二进制布局，只读取所需壁面；跳过单元对象避免重复全域检查。
    def i():
        nonlocal pos
        v=struct.unpack_from('<i',data,pos)[0];pos+=4;return v
    def arr():
        nonlocal pos
        n=i();a=np.frombuffer(data,dtype='<i4',count=n,offset=pos);pos+=4*n;return a
    def string():
        nonlocal pos
        n=i();s=data[pos:pos+n].decode();pos+=n;return s
    string();dim,nv,nf,nc,nb=[i() for _ in range(5)];assert dim==3
    points=np.frombuffer(data,dtype='<f8',count=3*nv,offset=pos).reshape(nv,3);pos+=24*nv
    offsets=np.empty(nf,dtype=np.int64)
    for k in range(nf):offsets[k]=pos;pos+=4*(1+struct.unpack_from('<i',data,pos)[0])
    for k in range(nc):
        pos+=4
        for _ in range(3):pos+=4*(1+struct.unpack_from('<i',data,pos)[0])
    wall=None
    for k in range(nb):
        name=string();ids,signs=arr(),arr()
        if name==m['condition']['boundaryMarkers']['wall']:wall=(ids,signs)
    assert pos==len(data) and wall is not None
    sums=np.zeros(3);area=0.;minarea=float('inf');perimeter_keys={}
    for fid,sign in zip(*wall):
        at=offsets[fid];n=struct.unpack_from('<i',data,at)[0];assert n==3
        ids=np.frombuffer(data,dtype='<i4',count=3,offset=at+4);p=points[ids]
        av=np.cross(p[1]-p[0],p[2]-p[0])/2
        a=float(np.linalg.norm(av));assert a>0;area+=a;minarea=min(minarea,a);sums+=sign*av
        for x,y in zip(ids,np.roll(ids,-1)):
            key=tuple(sorted((int(x),int(y))));perimeter_keys[key]=perimeter_keys.get(key,0)+1
    assert len(wall[0])==19752 and abs(area-1.9646292565715551)<2e-12
    assert abs(area-m['mesh']['wallAreaM2'])<2e-12
    assert max(abs(sums))<2e-12 and set(perimeter_keys.values())=={2}
    r=dict(passed=True,wallFaces=len(wall[0]),nativeWallAreaM2=area,minimumFaceAreaM2=minarea,
           signedAreaClosureVectorM2=sums.tolist(),wallEdges=len(perimeter_keys),allWallEdgesTwoFaces=True,
           nativeGridSHA256=hashlib.sha256(data).hexdigest(),elapsedSeconds=time.monotonic()-start,
           newQualifiedAreaM2=0)
    dest.write_text(json.dumps(r,indent=2)+'\n');return r


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(verify(a.case,a.out)))
