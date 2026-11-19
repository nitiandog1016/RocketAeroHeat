"""求解前核对两份网格的共同节点、冷启动和 STL 面元，绑定本次输入。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"scripts"))
from cm011Tunnel import write_json
from runEmbeddedSU2 import sha


def check(mesh, out):
    prior=ROOT/"outputs/v276_validation/gui01/run_folder.txt"
    coarse=Path(prior.read_text("utf-8-sig").strip())/"flow"
    receipt=json.loads((coarse/"engine_process.json").read_bytes())
    assert receipt["status"]=="terminal" and receipt["exitCode"]==0 and not receipt["timedOut"]
    assert sha(coarse/"restart.csv")==receipt["outputSHA256"]["restart.csv"]=="01d0309309e29ab791c800e74544eb31b71a71415aefb65eb372c0350d8089f4"
    before,after=[json.loads((d/"geometry.json").read_bytes()) for d in (coarse,mesh)]
    assert sha(coarse/"body.json")==sha(mesh/"body.json")=="de44b9ac4b3a74ac530f7c339dda8d31c62fb4f872a9183616188f6d6a65d2b2"
    assert sha(coarse/"case.cfg")==sha(mesh/"case.cfg")=="f93b190bf00a68979c02d56c59ccd2fdb51a3856356ad0c00971f79bd5660b34"
    stl=ROOT/"cases/CM011-alternate-SI.stl"
    assert sha(stl)==after["importedSurface"]["sourceSHA256"]==before["importedSurface"]["sourceSHA256"]
    assert after["rayRadians"][::2]==before["rayRadians"] and len(after["rayRadians"])==25
    assert after["nodes"]==29425 and after["hexahedra"]+after["prisms"]==27648
    assert after["meridianNodes"]==before["meridianNodes"]==49 and after["radialLayers"]==before["radialLayers"]==24
    old,new=[np.load(d/"nodes.npz") for d in (coarse,mesh)]
    np.testing.assert_array_equal(new["xyz"][new["ids"][::2]],old["xyz"][old["ids"]])
    initial=[np.loadtxt(d/"initial_flow.csv",delimiter=",",skiprows=1) for d in (coarse,mesh)]
    for grid,field in zip((old,new),initial):
        np.testing.assert_array_equal(field[:,0],np.arange(len(field)))
        np.testing.assert_array_equal(field[:,1:4],grid["xyz"])
    np.testing.assert_array_equal(initial[0][old["ids"],4:],initial[1][new["ids"][::2],4:])
    np.testing.assert_array_equal(initial[1][new["ids"][1::2],4:],initial[0][old["ids"][:-1],4:])
    write_json(out,dict(originalFaces=after["originalFaces"],baselineRestart=sha(coarse/"restart.csv"),
                        commonNodes=before["nodes"],refinedNodes=after["nodes"],identicalCommonCoordinates=True,
                        initialValuesVerified=True,stlSHA256=sha(stl),preparedSHA256=sha(mesh/"prepared.json"),
                        inputs=json.loads((mesh/"prepared.json").read_bytes())))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh",type=Path)
    parser.add_argument("receipt",type=Path)
    args=parser.parse_args()
    check(args.mesh,args.receipt)
