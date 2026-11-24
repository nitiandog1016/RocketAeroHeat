"""用删去测量列的输入验证物性，并比较独立生成的等温冷网与已登记的原几何。"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/"scripts"))
from shortSupportFilmProperties import film_state,native_text,BASE_CONFIG,CASE,read
from shortSphere65 import sha


def check(prepared=None):
    case=read(CASE);keep=("airF","diameterIn","speedFtS","pressureLbfFt2","waterMassFraction","reynoldsFree")
    physical={key:case[key] for key in keep}
    decoy=dict(case,meanSurfaceF=-5000,surfaceF=[9999]*7,localNuSurface=[-100]*7)
    configs={}
    for excess in (20.,30.,40.):
        c,f=film_state(physical,excess)
        assert (c,f)==film_state(case,excess)==film_state(decoy,excess)
        data=dict(condition=c,filmTransport=f,uniformExcessK=excess)
        configs[str(excess)]={}
        for steps,release in BASE_CONFIG:
            txt=native_text(data,steps,release)
            configs[str(excess)][f"{steps}-{'release' if release else 'solve'}"]=hashlib.sha256(txt.encode("ascii")).hexdigest()
        try:native_text(dict(data,condition=dict(c,pe=1.)),200,False)
        except AssertionError:pass
        else:raise AssertionError("错误 Pe 未被拒绝")
    proof=dict(inputs=list(keep),candidateConfigHashes=configs,changedMeasuredValuesHaveNoEffect=True,wrongPeRejected=True)
    if prepared:
        p=read(prepared/"prepared.json");ref=ROOT/"outputs/v332_validation/prepared01"
        c,f=film_state(physical,20.)
        assert p["condition"]==c and p["equatorK"]==f["wallK"]
        for name in ("mesh.su2","transfer.npz"):
            assert sha(prepared/name)==sha(ref/name)==p["files"][name]
        with np.load(prepared/"mapping.npz") as cold,np.load(ref/"mapping.npz") as orig:
            for key in cold.files:
                if key!="theta":np.testing.assert_array_equal(cold[key],orig[key])
            np.testing.assert_allclose(cold["theta"],20/90,atol=1e-15,rtol=0)
            state=np.loadtxt(prepared/"initial.csv",skiprows=1,delimiter=",");wet=cold["wet"]
            scalar=np.zeros(len(state));scalar[wet]=20/90
            np.testing.assert_allclose(state[:,9],scalar,atol=2e-15,rtol=0)
            np.testing.assert_allclose(state[:,8],1-298.15+scalar,atol=6e-14,rtol=0)
            velocity=np.zeros((len(state),3));velocity[:,0]=-1;velocity[wet]=0
            np.testing.assert_array_equal(state[:,5:8],velocity)
        proof["coldGrid"]=dict(nodes=len(state),wetNodes=len(wet),mesh=p["files"]["mesh.su2"],transfer=p["files"]["transfer.npz"],
                               initial=p["files"]["initial.csv"],properties=f)
    return proof


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--prepared",type=Path);parser.add_argument("--record",type=Path);a=parser.parse_args()
    result=check(a.prepared)
    if a.record:
        a.record.parent.mkdir(exist_ok=True)
        with a.record.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2);f.write("\n")
    print("Thermal properties verified with held-out experimental temperatures.")
