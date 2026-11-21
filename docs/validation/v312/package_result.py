"""核验 20D 与 30D 的原生帧、输入和七区信息，打包工作台使用的全周对照数据。"""
import argparse
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha
from thermistorReaction import terminal


def package(prepared,run,out):
    info=read(prepared/"mesh_audit.json");rec=terminal(run);s=read(run/"surface/summary.json");source=read(ROOT/"cases/thermistor_b/sidewind3.json")
    assert all(sha(prepared/k)==h for k,h in info["files"].items())
    assert info["baseMeshSHA256"]==source["meshSHA256"] and info["seedSHA256"]==source["sourceFrameSHA256"]
    assert sha(run/"mesh.su2")==info["files"]["mesh.su2"]==s["meshSHA256"]
    assert s["bodySHA256"]==source["bodySHA256"]==sha(run/"body.json")
    assert s["sourceFrameSHA256"]==rec["outputSHA256"]["restart.csv"]==s["energy"]["inputFrameSHA256"]
    assert sha(run/"surface/wall.csv")==s["wallSHA256"] and read(run/"case.json")["outerRadiusOverD"]==30
    assert s["relativeChange20"]<.002 and abs(s["energy"]["relativeDifference"])<.005 and s["faces"]==10272
    assert max(s["residualWindow"].values())<-8 and max(s["plus20Window"].values())<-8
    target=out/"sidewindOuter30.zip";manifest=out/"sidewindOuter30.json";assert not manifest.exists()
    with ZipFile(target,"x") as z:
        for path in (run/"body.json",run/"case.json",run/"surface/wall.csv",run/"surface/summary.json"):
            z.writestr(ZipInfo(path.name,(2026,1,1,0,0,0)),path.read_bytes(),compress_type=ZIP_DEFLATED)
    key=dict(baseMeshSHA256=source["meshSHA256"],baseCaseSHA256=source["caseSHA256"],seedFrameSHA256=info["seedSHA256"],
             bodySHA256=s["bodySHA256"],meshSHA256=s["meshSHA256"],caseSHA256=s["caseSHA256"],
             sourceFrameSHA256=s["sourceFrameSHA256"],wallSHA256=s["wallSHA256"],archiveSHA256=sha(target))
    manifest.write_bytes((json.dumps(key,indent=2)+"\n").encode("utf8"))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("prepared","run","out"):p.add_argument(key,type=Path)
    a=p.parse_args();package(a.prepared.resolve(),a.run.resolve(),a.out.resolve())
