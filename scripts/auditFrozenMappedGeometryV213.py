"""模块职责：按固定原生分解复算有向体积，判断节点约定，不运行CFD。"""
from pathlib import Path
import hashlib,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/v213_validation'
PERM=np.array([1,0,2,4,3,5])
FACES=((3,5,2,0),(1,2,5,4),(0,1,4,3),(0,2,1),(3,4,5))

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def tet(a,b,c,d):return np.einsum('ij,ij->i',d-a,np.cross(b-a,c-a))/6

def volume(p):
    # 原生wedge_properties：三四棱锥和两四面体；每四棱锥分成四四面体。
    center=p.mean(axis=1);v=np.zeros(len(p))
    for f in FACES:
        if len(f)==3:v+=tet(p[:,f[0]],p[:,f[1]],p[:,f[2]],center)
        else:
            mid=p[:,f].mean(axis=1)
            for i in range(4):v+=tet(p[:,f[i]],p[:,f[(i+1)%4]],mid,center)
    return v

def main():
    start=time.monotonic();case=BASE/'control01';z=np.load(str(case/'mapping.npz'))
    donors=json.loads((case/'donors.json').read_text())
    sources=[ROOT/'outputs/v209_validation/source04/source/src/geom/elements/properties.d',
             BASE/'source02/source/src/geom/grid/usgrid.d',case/'mapping.npz',case/'donors.json',Path(__file__).resolve()]
    results={}
    for name,p in [('receptor',z['points'][z['cells']]),('donor',np.array([r['vertices'] for r in donors['donors']]))]:
        a=volume(p);b=volume(p[:,PERM]);assert np.all(a<0) and np.all(b>0)
        err=float(np.max(np.abs(a+b)/b));assert err<1e-12
        faces=lambda row:{tuple(sorted(row[list(f)])) for f in FACES}
        assert faces(np.arange(6))==faces(PERM)
        results[name]=dict(count=len(p),originalNegative=int((a<0).sum()),
            convertedPositive=int((b>0).sum()),minimumConvertedM3=float(b.min()),
            maximumConvertedM3=float(b.max()),relativeMagnitudeDifference=err)
    result=dict(passed=True,results=results,permutation=PERM.tolist(),
        sourceSHA256={str(p.relative_to(ROOT)):sha(p) for p in sources},
        physicalGeometryChanged=False,CFDRuns=0,seconds=time.monotonic()-start,
        limit='offline reproduction proves invalid native orientation; no captured failing native cell ID')
    out=BASE/'geometry_order_audit.json';assert not out.exists()
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
