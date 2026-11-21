"""从前引线细网终帧提取共同位置温度，与原面数据一起生成离线工作台参照包。"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha
from thermistorReaction import terminal


def package(prepared,run,dest):
    info=read(prepared/"mesh_audit.json")
    key=read(ROOT/"cases/thermistor_b/sidewind3.json")
    rec=terminal(run)
    s=read(run/"surface/summary.json")
    c=read(run/"case.json")
    assert info["baseMeshSHA256"]==key["meshSHA256"] and info["commonPositionNodesExact"]
    assert all(sha(prepared/n)==v for n,v in info["files"].items())
    for name,field in (("case.json","caseSHA256"),("mesh.su2","meshSHA256"),("body.json","bodySHA256")):
        assert sha(run/name)==info["files"][name]==s[field]
    assert s["bodySHA256"]==key["bodySHA256"] and s["wallSHA256"]==sha(run/"surface/wall.csv")
    assert s["sourceFrameSHA256"]==rec["outputSHA256"]["restart.csv"]==s["energy"]["inputFrameSHA256"]
    assert s["sourceReactionSHA256"]==terminal(run/"reaction0")["outputSHA256"]["restart.csv"]
    assert s["energy"]==read(run/"ledger/energy.json") and abs(s["energy"]["relativeDifference"])<.005
    assert s["relativeChange20"]<.002 and max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8
    assert c["frontCapFirstLayerRatio"]==.5 and c["velocityMps"]==[0,-3,0] and c["outerRadiusOverD"]==20
    with np.load(run/"nodes.npz") as z:
        xyz,wall,near,common=(z[k] for k in ("xyz","wallMap","firstMap","commonMap"))
    body=read(run/"body.json")
    tri=np.asarray(body["V"])[np.asarray(body["F"])-1]
    normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    state=np.loadtxt(run/"restart.csv",delimiter=",",skiprows=1,usecols=(9,))
    def height(ids):
        return c["diameterM"]*np.einsum("fvi,fi->f",xyz[ids]-xyz[wall],normal)/3
    def temperature(ids):
        return c["airK"]+(c["wallK"]-c["airK"])*state[ids].mean(axis=1)
    matched=temperature(common)
    rows=np.genfromtxt(run/"surface/wall.csv",delimiter=",",names=True)
    np.testing.assert_allclose(temperature(near),rows["nearWallCellTemperatureK"],rtol=0,atol=1e-12)
    np.testing.assert_array_equal(rows["faceID"],np.arange(10272))
    assert min(matched)>=c["airK"] and max(matched)<=c["wallK"]
    capture=io.StringIO()
    np.savetxt(capture,np.column_stack((rows["faceID"],matched,height(common),height(near))),
               delimiter=",",fmt=["%d","%.17g","%.17g","%.17g"],comments="",
               header="faceID,matchedK,baseDistanceM,fineDistanceM")
    sample=capture.getvalue().encode("ascii")
    dest.mkdir(parents=True,exist_ok=True)
    manifest=dest/"sidewindTipSpacing.json"
    assert not manifest.exists()
    archive=dest/"sidewindTipSpacing.zip"
    with ZipFile(archive,"x") as z:
        for name,blob in (("matched_temperature.csv",sample),("case.json",(run/"case.json").read_bytes()),
                          ("body.json",(run/"body.json").read_bytes()),("wall.csv",(run/"surface/wall.csv").read_bytes()),
                          ("summary.json",(run/"surface/summary.json").read_bytes())):
            z.writestr(ZipInfo(name,(2026,1,1,0,0,0)),blob,compress_type=ZIP_DEFLATED)
    identity=dict(baseCaseSHA256=key["caseSHA256"],baseMeshSHA256=key["meshSHA256"],seedFrameSHA256=info["seedSHA256"],
                  sampleSHA256=hashlib.sha256(sample).hexdigest(),archiveSHA256=sha(archive))
    identity.update({k:s[k] for k in ("caseSHA256","meshSHA256","bodySHA256","sourceFrameSHA256","wallSHA256","sourceReactionSHA256")})
    manifest.write_bytes((json.dumps(identity,indent=2)+"\n").encode("utf8"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("prepared","run","dest"):p.add_argument(k,type=Path)
    a=p.parse_args();package(a.prepared.resolve(),a.run.resolve(),a.dest.resolve())
