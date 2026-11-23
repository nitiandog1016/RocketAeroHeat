"""独立读回两套银球体网，核对实际父子壁片、原面分数、温界和测线。"""
import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from shortSphere65 import CASE, kelvin, read, sha, write_json
from thermistorOuterDomain import load_mesh

BASE=ROOT/"outputs/v330_validation/prepared01"
BASE_KEY="72fc2c60e2ea9aaa4a9bbf924dd54c8ae9a930c0a5b4961005cea25e0de15902"
NEW_KEY="22ef3bf27abf1a093fae6bfb03145fac1b9f401229f11a3d12042be19cd9a04c"


def inspect(prep):
    for path,key in ((BASE,BASE_KEY),(prep,NEW_KEY)):
        assert sha(path/"prepared.json")==key
        spec=read(path/"prepared.json")
        assert all(sha(path/n)==h for n,h in spec["files"].items()) and sha(CASE)==spec["caseSHA256"]
    old=read(BASE/"prepared.json");g=spec["geometry"];s=spec["condition"]
    assert old["condition"]==s and old["supportedSHA256"]==spec["supportedSHA256"]
    cc,xc,bc=load_mesh(BASE/"mesh.su2",("wall","far"));fc,xf,bf=load_mesh(prep/"mesh.su2",("wall","far"))
    dist,match=cKDTree(xf).query(xc)
    assert len(set(match))==len(xc) and max(dist)<1e-13
    old_set=set(match.tolist());nw=len(cc)//64
    assert len(fc)==2*len(cc)==397824 and len(xf)==261708
    for j,cell in enumerate(cc):
        sector,loc=divmod(j,nw)
        a=set(fc[2*sector*nw+loc]);b=set(fc[(2*sector+1)*nw+loc])
        assert ((a|b)&old_set)==set(match[list(cell)]) and len(a&b)>=3
    old_far=set(i for f in bc["far"] for i in f);new_far=set(i for f in bf["far"] for i in f)
    assert set(match[list(old_far)])<=new_far and len(new_far)==g["farNodes"]==2690
    np.testing.assert_allclose(np.linalg.norm(xf[list(new_far)],axis=1),14,rtol=0,atol=1e-11)
    with np.load(BASE/"mapping.npz") as m,np.load(prep/"mapping.npz") as f:
        wet=f["wet"];wall=xf[wet];corners=f["corners"];groups=f["group"];first=f["first"]
        np.testing.assert_array_equal(xf,f["xyz"]);np.testing.assert_array_equal(xc,m["xyz"])
        coarse=xc[m["corners"]];coarse_first=xc[m["first"]][np.searchsorted(m["wet"],m["corners"])]
        def compare_children(parent,children):
            a=(parent[:,0]+parent[:,1])/2;b=(parent[:,2]+parent[:,3])/2
            expected=np.stack((parent[:,0],a,b,parent[:,3],a,parent[:,1],parent[:,2],b),axis=1)
            actual=children.reshape(64,2,94,4,3).transpose(0,2,1,3,4).reshape(-1,8,3)
            delta=np.linalg.norm(actual-expected,axis=2).max()
            assert delta<1e-13
            return delta
        wall_diff=compare_children(coarse,xf[corners])
        first_diff=compare_children(coarse_first,xf[first][np.searchsorted(wet,corners)])
        native=f["nativeArea"];parent=np.repeat(np.arange(6016).reshape(64,94),2,axis=0).ravel()
        np.testing.assert_allclose(native,m["nativeArea"][parent]/2,rtol=2e-12,atol=1e-20)
        np.testing.assert_array_equal(groups,m["group"][parent]);np.testing.assert_array_equal(match[m["first"]],first[np.searchsorted(wet,match[m["wet"]])])
        np.testing.assert_allclose(m["theta"],f["theta"][np.searchsorted(wet,match[m["wet"]])],rtol=0,atol=2e-14)
        sphere=np.unique(corners[groups==1])
        angle=np.rad2deg(np.arctan2(np.linalg.norm(xf[sphere,1:],axis=1),xf[sphere,0]))
        case=read(CASE);ref=np.full(len(wet),spec["equatorK"])
        ref[np.searchsorted(wet,sphere)]=np.interp(angle,case["angleDeg"],kelvin(case["surfaceF"]))
        np.testing.assert_allclose(ref,s["airK"]+s["scalarUnitK"]*f["theta"],atol=1e-12,rtol=0)
        line=sphere[abs(xf[sphere,1])<1e-12];assert len(line)==128
        phi=np.rad2deg(np.arctan2(xf[line,2],xf[line,0]))
        for a in (143.4375,146.25,149.0625):assert min(abs(phi-a))<1e-11
    with np.load(prep/"transfer.npz") as f,np.load(BASE/"transfer.npz") as m:
        ids,src,wt,ws=[f[k] for k in ("target","source","targetFraction","sourceFraction")]
        np.testing.assert_array_equal(f["region"][ids],groups[src]);np.testing.assert_array_equal(f["targetArea"],m["targetArea"])
        assert wt.min()>0 and ws.min()>0 and len(set(zip(ids,src)))==len(src)
        for ix,share in ((ids,wt),(src,ws)):np.testing.assert_allclose(np.bincount(ix,share),1,rtol=0,atol=1e-7)
        shape=(len(f["targetArea"]),len(coarse))
        grouped=csr_matrix((wt,(ids,parent[src])),shape)
        baseline=csr_matrix((m["targetFraction"],(m["target"],m["source"])),shape)
        gap=abs(grouped-baseline).max();assert gap<1e-7
        delta=[]
        for i in range(1,6):
            # 用异号、非均匀的源壁片功率检查五区逐原面传递。
            watts=native*(np.arange(len(native))%11-5)
            recv=np.bincount(ids,ws*watts[src])
            delta.append(abs(recv[f["region"]==i].sum()-watts[groups==i].sum()))
        assert max(delta)<1e-14
    assert g["minCellVolumeD3"]>0 and g["minConeVolumeD3"]>0
    return dict(preparedSHA256=NEW_KEY,oldPreparedSHA256=BASE_KEY,meshSHA256=spec["files"]["mesh.su2"],
                coldSHA256=spec["files"]["initial.csv"],geometry=g,oldNodeMaxErrorD=max(dist),oldCells=len(cc),
                wallMidpointMaxErrorD=wall_diff,firstMidpointMaxErrorD=first_diff,
                maxParentTargetFractionDifference=gap,maxRegionSignedTrialPowerErrorW=max(delta),
                outerDomainRelativeVolumeChange=g["volumeD3"]/old["geometry"]["volumeD3"]-1)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("prepared",type=Path);parser.add_argument("output",type=Path)
    a=parser.parse_args();result=inspect(a.prepared);write_json(a.output,result);print(result)
