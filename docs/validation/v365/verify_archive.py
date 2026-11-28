"""独立用三角形顶点均值和底面径向积分复核两组映射、圆柱热流及反向面积。"""
import hashlib
import json
from pathlib import Path
import numpy as np

DOC=Path(__file__).resolve().parent
CASE=DOC.parents[2]/"cases/simple_launcher"


def verify():
    folder=DOC/"research"
    s=json.loads((folder/"surface.json").read_bytes())
    geom=json.loads((CASE/"geometry.json").read_bytes())
    stl=CASE/geom["file"]
    assert hashlib.sha256(stl.read_bytes()).hexdigest()==s["surfaceSHA256"]==geom["fileSHA256"]
    dtype=np.dtype([("normal","<f4",3),("v","<f4",(3,3)),("att","<u2")])
    tri=np.fromfile(stl,dtype=dtype,offset=84)["v"].astype(float)
    area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    table=np.genfromtxt(folder/"wall.csv",delimiter=",",names=True)
    x=table["xM"];r=table["rM"];base=x==0
    xside=x[np.flatnonzero(base)[-1]:]
    flat=np.max(tri[:,:,0],axis=1)==0
    # 原始壁以直径为尺度做过归一化，米制还原最多可能使端点移动一个浮点位。
    xcorners=tri[~flat,:,0]
    node=np.abs(xcorners[:,:,None]-xside).argmin(axis=2)
    np.testing.assert_allclose(xside[node],xcorners,rtol=0,atol=1e-12)
    sel=np.arange(2016,6912)
    assert not flat[:sel[-1]+1].any() and abs(area[sel].sum()-geom["regions"][2]["areaM2"])<1e-12
    a=area[sel];face_group=np.array(s["faceGroup"])-1
    ratio=r[base]/r[base].max()
    rad0,rad1=ratio[:-1],ratio[1:]
    # 对底面的线性节点插值按 r·dr 积分。
    w0=(rad1-rad0)*(2*rad0+rad1)/3
    w1=(rad1-rad0)*(rad0+2*rad1)/3
    q_fields=[];rev_both=np.ones(len(sel),bool)
    for grid,tag in (("normal4","4um"),("normal2","2um")):
        raw=np.column_stack([table[f"{key}_{tag}"] for key in s["columns"]])
        values=np.empty((len(area),len(s["columns"])))
        values[~flat]=raw[np.flatnonzero(base)[-1]+node].mean(axis=1)
        values[flat]=w0@raw[base][:-1]+w1@raw[base][1:]
        expected=np.array(s[grid])[face_group]
        np.testing.assert_allclose(values,expected,rtol=1e-11,atol=1e-9)
        q=values[sel,1];q_fields.append(q)
        np.testing.assert_allclose(a@q,s["cylinder"][grid]["cylinderPowerW"],atol=1e-9)
        idx=node[sel]+np.flatnonzero(base)[-1]
        reverse=(raw[idx,2]<0).all(axis=1)&(raw[idx,3]<0).all(axis=1)
        assert int(reverse.sum())==s["cylinder"][grid]["reverseCylinderFaces"]
        rev_both &= reverse
    rms=np.sqrt(a@(q_fields[0]-q_fields[1])**2/(a@q_fields[1]**2))
    np.testing.assert_allclose(rms,s["gridDifference"]["heatRelativeRms"])
    np.testing.assert_allclose(a[rev_both].sum(),s["commonReverseCylinderAreaM2"])
    ranges=np.concatenate([np.arange(lo,hi+1) for lo,hi in s["commonReverseCylinderFaceRanges"]]) if rev_both.any() else []
    np.testing.assert_array_equal(ranges,sel[rev_both]+1)
    assert s["wholeBodyPowerW"] is None and not s["physicalQualification"]
    # 圆柱轴向两端的独立贡献帮助区分局部缺口与整段平均值。
    spans,bands=np.unique(np.sort(tri[sel,:,0],axis=1)[:,[0,-1]],axis=0,return_inverse=True)
    square=a*(q_fields[0]-q_fields[1])**2
    fractions=np.bincount(bands,square)/square.sum()
    return dict(checkedFaces=len(area),checkedFields=len(s["columns"]),grids=2,commonReverseFaces=int(rev_both.sum()),
                shoulderEndSquaredDifferenceShare=float(fractions[0]),noseEndSquaredDifferenceShare=float(fractions[-1]),
                shoulderEndXM=spans[0].tolist(),noseEndXM=spans[-1].tolist())


if __name__=="__main__":
    print(json.dumps(verify()))
