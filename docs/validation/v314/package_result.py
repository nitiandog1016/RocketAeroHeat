"""绑定带线热敏珠 96 档终场、几何映射和 48 档参照，发布逐原面数据。"""
import argparse
import json
from pathlib import Path
import sys
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha
from thermistorReaction import terminal


def package(prepared,run,dest):
    audit=read(prepared/"mesh_audit.json");base=read(ROOT/"cases/thermistor_b/sidewindTipSpacing.json")
    s=read(run/"surface/summary.json");rec=terminal(run)
    assert audit["commonPositionsExact"] and audit["baseMeshSHA256"]==base["meshSHA256"]
    assert all(sha(prepared/n)==h for n,h in audit["files"].items())
    for n,k in (("nodes.npz","nodeMapSHA256"),("mesh.su2","meshSHA256"),("case.json","caseSHA256"),("body.json","bodySHA256")):
        assert sha(run/n)==audit["files"][n]==s[k]
    assert s["bodySHA256"]==base["bodySHA256"] and s["wallSHA256"]==sha(run/"surface/wall.csv")
    assert s["sourceFrameSHA256"]==rec["outputSHA256"]["restart.csv"]==s["energy"]["inputFrameSHA256"]
    assert s["sourceReactionSHA256"]==terminal(run/"reaction0")["outputSHA256"]["restart.csv"]
    assert s["plus20FrameSHA256"]==terminal(run/"plus20")["outputSHA256"]["restart.csv"]
    assert s["reaction20SHA256"]==terminal(run/"reaction20")["outputSHA256"]["restart.csv"]
    assert all(sha(ROOT/k)==h for k,h in s["postprocessSHA256"].items())
    assert s["relativeChange20"]<.002 and max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8
    assert s["energy"]==read(run/"ledger/energy.json") and abs(s["energy"]["relativeDifference"])<.005
    assert s["angularCells"]==96 and s["nativeWallTriangles"]==audit["nativeWallTriangles"]==20544
    assert s["faces"]==audit["reportFaces"]==10272 and s["jointExperimentVerifiedM2"]==0
    dest.mkdir(parents=True,exist_ok=True)
    key=dest/"sidewindAngular96.json";archive=key.with_suffix(".zip");assert not key.exists()
    with ZipFile(archive,"x") as z:
        for name,path in (("body.json",run/"body.json"),("case.json",run/"case.json"),("wall.csv",run/"surface/wall.csv"),
                          ("summary.json",run/"surface/summary.json"),("comparison.json",prepared/"mesh_audit.json")):
            z.writestr(ZipInfo(name,(2026,1,1,0,0,0)),path.read_bytes(),compress_type=ZIP_DEFLATED)
    ident={k:s[k] for k in ("caseSHA256","meshSHA256","bodySHA256","nodeMapSHA256","sourceFrameSHA256","wallSHA256","sourceReactionSHA256")}
    ident.update(baseMeshSHA256=base["meshSHA256"],baseFrameSHA256=base["sourceFrameSHA256"],
                 baseWallSHA256=base["wallSHA256"],archiveSHA256=sha(archive),mappingAuditSHA256=sha(prepared/"mesh_audit.json"))
    key.write_bytes((json.dumps(ident,indent=2)+"\n").encode("utf8"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("prepared","run","dest"):p.add_argument(k,type=Path)
    a=p.parse_args();package(a.prepared.resolve(),a.run.resolve(),a.dest.resolve())
