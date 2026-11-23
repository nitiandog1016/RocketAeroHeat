"""重读两份原生网，检验每个旧近壁单元的两个子单元及全部外区连接。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from shortSupportMesh import COUNTS, SECTORS
from shortSphere65 import CASE, sha
from thermistorOuterDomain import load_mesh

COARSE=ROOT/"outputs/v325_validation/prep02/prepared"
COARSE_KEY="ca9c622e056924cb0482e06e5f3561b402d7c40df079da93dd62cd61a98587b0"
FINE_KEY="72fc2c60e2ea9aaa4a9bbf924dd54c8ae9a930c0a5b4961005cea25e0de15902"


def inspect(fine):
    specs=[]
    for path,key in ((COARSE,COARSE_KEY),(fine,FINE_KEY)):
        assert sha(path/"prepared.json")==key
        spec=json.loads((path/"prepared.json").read_bytes());specs.append(spec)
        assert sha(CASE)==spec["caseSHA256"] and all(sha(path/n)==h for n,h in spec["files"].items())
    coarseSpec,fineSpec=specs
    assert coarseSpec["condition"]==fineSpec["condition"] and coarseSpec["supportedSHA256"]==fineSpec["supportedSHA256"]
    assert coarseSpec["files"]["transfer.npz"]==fineSpec["files"]["transfer.npz"]
    cc,xc,bc=load_mesh(COARSE/"mesh.su2",("wall","far"));fc,xf,bf=load_mesh(fine/"mesh.su2",("wall","far"))
    dist,ids=cKDTree(xf).query(xc);assert max(dist)<1e-12 and len(set(ids))==len(xc)
    for marker in ("wall","far"):
        assert [tuple(ids[list(face)]) for face in bc[marker]]==bf[marker]
    original=set(ids.tolist());new=set(range(len(xf)))-original;seen=set()
    width=sum(COUNTS);ncoarse=len(cc)//SECTORS;nfine=len(fc)//SECTORS
    maxMid=0.;outerCount=0
    for k in range(SECTORS):
        c=cc[k*ncoarse:(k+1)*ncoarse];f=fc[k*nfine:(k+1)*nfine]
        assert len(f)==len(c)+5*width
        for j,cell in enumerate(c):
            mapped=ids[list(cell)]
            if j>=5*width:
                assert tuple(mapped)==f[j+5*width];outerCount+=1;continue
            layer,col=divmod(j,width);idx=2*layer*width+col
            lower,upper=set(f[idx]),set(f[idx+width])
            middle=lower & upper;ends_a=sorted(lower-middle);ends_b=sorted(upper-middle)
            assert set(mapped)==set(ends_a+ends_b) and middle<=new and len(middle)==len(ends_a)==len(ends_b)
            # 由实际顶点坐标反推中点的两端，检验每一条连线及一一对应。
            pairs=(xf[ends_a,None,:]+xf[None,ends_b,:])/2
            distances=np.linalg.norm(xf[sorted(middle),None,None,:]-pairs[None,:,:,:],axis=3)
            ai,bi=np.unravel_index(distances.reshape(len(middle),-1).argmin(1),pairs.shape[:2])
            assert len(set(ai))==len(middle)==len(set(bi))
            maxMid=max(maxMid,float(distances[np.arange(len(middle)),ai,bi].max()));seen.update(middle)
    assert seen==new and maxMid<1e-12
    with np.load(COARSE/"mapping.npz") as c,np.load(fine/"mapping.npz") as f:
        np.testing.assert_array_equal(ids[c["corners"]],f["corners"])
        np.testing.assert_array_equal(ids[c["wet"]],f["wet"])
        np.testing.assert_array_equal(c["theta"],f["theta"]);np.testing.assert_array_equal(c["nativeArea"],f["nativeArea"])
        np.testing.assert_allclose(xf[f["first"]],(xc[c["wet"]]+xc[c["first"]])/2,atol=1e-14,rtol=0)
        wet=f["wet"];cold=np.loadtxt(fine/"initial.csv",skiprows=1,delimiter=",")
        np.testing.assert_array_equal(cold[:,1:4],xf);mask=np.ones(len(xf),bool);mask[wet]=False
        assert np.all(cold[mask,9]==0) and np.all(cold[mask,5]==-1) and np.all(cold[wet,5]==0)
        np.testing.assert_array_equal(cold[wet,9],f["theta"])
        np.testing.assert_array_equal(np.loadtxt(fine/"profile.dat")[:,0],wet)
    geom=fineSpec["geometry"];assert geom["minConeVolumeD3"]>0 and geom["minCellVolumeD3"]>0
    return dict(preparedSHA256=FINE_KEY,meshSHA256=sha(fine/"mesh.su2"),coldSHA256=sha(fine/"initial.csv"),
                coarseNodes=len(xc),fineNodes=len(xf),coarseCells=len(cc),fineCells=len(fc),
                newMidpoints=len(new),bisectedCells=len(cc)-outerCount,untouchedOuterCells=outerCount,
                maxOldCoordinateErrorD=float(dist.max()),maxMidpointErrorD=maxMid,geometry=geom,transferSHA256=sha(fine/"transfer.npz"))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("prepared",type=Path);p.add_argument("output",type=Path)
    a=p.parse_args();result=inspect(a.prepared)
    with a.output.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2);f.write("\n")
    print(result["bisectedCells"],result["untouchedOuterCells"],result["maxMidpointErrorD"])
