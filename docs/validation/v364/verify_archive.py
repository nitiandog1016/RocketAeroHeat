"""从 Git 内的三角面和壁段表独立核对面积映射、绝对压正性与六区外域差。"""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[2]/"cases/simple_launcher"


def verify(folder, positive):
    source = json.loads((folder/"surface.json").read_bytes())
    geom = json.loads((MODEL/"geometry.json").read_bytes())
    assert source["surfaceSHA256"] == geom["fileSHA256"]
    dtype = np.dtype([("normal","<f4",3),("v","<f4",(3,3)),("att","<u2")])
    tri = np.fromfile(MODEL/geom["file"],dtype=dtype,offset=84)["v"].astype(float)
    area = np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    table = np.loadtxt(folder/"bands.csv",delimiter=",",skiprows=1)
    assert table.shape == (216,6)
    base = table[:,2] == 0
    xb = table[~base,:][:,[0,2]]
    r = table[base,:][:,[1,3]]/table[base,3].max()
    # 同一环的两轴向半带，在一张三角面内分别占 1/4、3/4（或反之）。
    expected = np.empty((len(area),2))
    for i,vertices in enumerate(tri):
        x = vertices[:,0]
        if max(x)==0:
            expected[i] = np.diff(r**2,axis=1).ravel() @ table[base,4:]
        else:
            segments = (xb[:,0] >= min(x)-1e-12) & (xb[:,1] <= max(x)+1e-12)
            assert segments.sum() == 2
            first = .75 if np.count_nonzero(x==min(x))==2 else .25
            expected[i] = np.array([first,1-first]) @ table[~base,4:][segments]
    ids = np.array(source["faceGroup"])-1
    actual = np.column_stack((source["groupCp14D"],source["groupCp22D"]))[ids]
    np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-11)
    pressure = source["freePressurePa"] + actual*source["dynamicPressurePa"]
    assert np.isfinite(pressure).all() and bool((pressure>0).all()) == positive
    for reg in [*geom["regions"],dict(name="all",firstFace=1,lastFace=len(area))]:
        s = slice(reg["firstFace"]-1,reg["lastFace"])
        rec = source["regions"][reg["name"]]
        rms = np.sqrt(np.average(np.diff(expected[s],axis=1).ravel()**2,weights=area[s]))
        np.testing.assert_allclose([rms,area[s].sum(),expected[s,1].min(),expected[s,1].max()],
            [rec["rmsDeltaCp"],rec["areaM2"],rec["cpMin"],rec["cpMax"]],atol=1e-11)
    assert source["qualifiedPressureFaces"]==source["heatFaces"]==0 and source["wholeBodyPowerW"] is None
    return dict(nonPositiveFaces=(pressure<=0).sum(axis=0).tolist(),
                minimumAbsolutePa=pressure.min(axis=0).tolist(),checkedFaces=len(area))


if __name__=="__main__":
    print(json.dumps({key:verify(HERE/key,positive) for key,positive in (("fields",False),("research",True))},indent=2))
