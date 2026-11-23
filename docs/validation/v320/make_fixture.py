"""按标准面坐标生成两处重剖的毫米制 B 珠样件和独立的文件面号真值。"""
import json
import sys
from pathlib import Path
from zipfile import ZipFile

import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/"v319"))
from make_fixture import FACET


def make(out):
    out.mkdir(parents=True,exist_ok=False)
    with ZipFile(ROOT/"cases/thermistor_b/leads3.zip") as z:body=json.loads(z.read("body.json"))
    triangles=np.asarray(body["V"])[np.asarray(body["F"])-1];targets=[];parents=[];pieces=[]
    for i,t in enumerate(triangles,1):
        if i in (183,184):
            center=np.array([.2,.3,.5])@t
            for k in range(3):
                targets.append([t[k],t[(k+1)%3],center]);parents.append(i);pieces.append(k+1)
        else:targets.append(t);parents.append(i);pieces.append(0)
    order=np.random.RandomState(320).permutation(len(parents))
    records=np.zeros(len(parents),dtype=FACET);records["xyz"]=1000*np.asarray(targets)[order]
    xyz=records["xyz"].astype(float);v=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
    records["normal"]=v/np.linalg.norm(v,axis=1)[:,None]
    with (out/"probe-partition-mm.stl").open("xb") as f:
        f.write(b"B thermistor with split lead and joint facets | mm".ljust(80));f.write(len(order).to_bytes(4,"little"));records.tofile(f)
    np.savetxt(out/"face_truth.csv",np.column_stack((np.asarray(parents)[order],np.asarray(pieces)[order])),delimiter=",",fmt="%d")
    info=[dict(stlFaceID=k+1,parent=int(parents[i]),piece=int(pieces[i]),sourceM=triangles[parents[i]-1].tolist(),
               inputM=(xyz[k]*.001).tolist(),nativeFaces=([363,365,366] if parents[i]==183 else [367,368,370]))
          for k,i in enumerate(order) if pieces[i]]
    (out/"geometry.json").write_text(json.dumps(info,indent=2),encoding="utf8")


if __name__=="__main__":make(Path(sys.argv[1]))
