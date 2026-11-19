"""从同一 CM-011 解析外形另作一套顶点和三角面；不调用导入配准代码。"""
import argparse
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from cm011Tunnel import D, meridian


def generate(path):
    dense,_ = meridian(400)
    azimuths = 72
    j = np.arange(azimuths)
    theta = 2*np.pi*j/azimuths
    points = [[dense[0,0],0.,0.]]
    # 各区单独设置剖分密度。内圈顶点沿两维错开，保留四区实际边界。
    rings = []
    for start,end,count in ((0,5,8),(5,23,24),(23,27,7),(27,34,10),(34,48,22)):
        for k in range(1,count+1):
            if end==48 and k==count: continue
            row = len(rings)+1
            shift = .19*np.sin(3*theta+.61*row) if k<count else np.zeros_like(theta)
            s = start+(k+shift)*(end-start)/count
            x,r = dense[np.rint(400*s).astype(int)].T
            angle = theta+(.31+.16*np.sin(row*.85)+.09*np.sin(2*theta+row))*2*np.pi/azimuths
            rings.append(len(points)+j)
            points.extend(np.column_stack((x,r*np.cos(angle),r*np.sin(angle))))
    pole = len(points);points.append([dense[-1,0],0.,0.])
    faces = [(0,rings[0][(t+1)%azimuths],rings[0][t]) for t in j]
    for row,(a,b) in enumerate(zip(rings,rings[1:])):
        for t in j:
            n = (t+1)%azimuths
            if (row+t)%2:
                faces.extend([(a[t],a[n],b[n]),(a[t],b[n],b[t])])
            else:
                faces.extend([(a[t],a[n],b[t]),(a[n],b[n],b[t])])
    faces.extend((rings[-1][t],rings[-1][(t+1)%azimuths],pole) for t in j)
    v=np.array(points)*D
    f=np.asarray(faces)[np.random.RandomState(276).permutation(len(faces))]
    tri=v[f]
    normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    records=np.zeros(len(tri),dtype=[("normal","<f4",3),("xyz","<f4",(3,3)),("attr","<u2")])
    records["xyz"]=tri;records["normal"]=normals
    with path.open("xb") as out:
        out.write(b"CM-011 low-speed model; independent triangulation; coordinates in metres".ljust(80,b" "))
        out.write(len(tri).to_bytes(4,"little"));out.write(records.tobytes())
    print(len(v),len(f))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("new_stl",type=Path)
    generate(parser.parse_args().new_stl)
