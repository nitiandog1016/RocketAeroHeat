"""为 B 珠风洞测点准备实际的 96 档体网；没有本机缓存时从名义几何重建。"""
import argparse
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
from thermistorBead import ROOT,read,sha,write_json
from thermistorLeadTipSpacing import build_geometry as refine_tip
from thermistorAngularSpacing import build_geometry as split_angles
from thermistorSidewindWorkbench import prepare as prepare48


def prepare(out,rebuild=False):
    case_dir=ROOT/"cases/thermistor_b";key=read(case_dir/"sidewindAngular96.json")
    cache=ROOT/"outputs/v314_validation/prepared01"
    files={"mesh.su2":"meshSHA256","nodes.npz":"nodeMapSHA256","case.json":"caseSHA256","body.json":"bodySHA256"}
    reuse=not rebuild and all((cache/n).is_file() and sha(cache/n)==key[k] for n,k in files.items())
    if reuse:
        out.mkdir(parents=True,exist_ok=False)
        for n in files:
            try:os.link(cache/n,out/n)
            except OSError:shutil.copyfile(cache/n,out/n)
    else:
        out.parent.mkdir(parents=True,exist_ok=True)
        with TemporaryDirectory(prefix="bead96_",dir=str(out.parent)) as work:
            tmp=Path(work);prepare48(tmp/"base");refine_tip(tmp/"base",tmp/"tip");split_angles(tmp/"tip",out)
    assert all(sha(out/n)==key[k] for n,k in files.items())
    archive=case_dir/"sidewindAngular96.zip";assert sha(archive)==key["archiveSHA256"]
    with ZipFile(archive) as z:(out/"mesh_audit.json").write_bytes(z.read("comparison.json"))
    assert sha(out/"mesh_audit.json")==key["mappingAuditSHA256"]
    audit=read(out/"mesh_audit.json");body=read(out/"body.json")
    assert audit["reportFaces"]==len(body["F"])==10272 and audit["nativeWallTriangles"]==20544
    with np.load(out/"nodes.npz") as z:
        wall=z["wallMap"];vertices,index=np.unique(wall,return_inverse=True)
        preview=dict(V=(z["xyz"][vertices]*read(out/"case.json")["diameterM"]).tolist(),F=(index.reshape(-1,3)+1).tolist())
    write_json(out/"prepared.json",dict(directory=str(out),nodes=audit["nodes"],cells=audit["cells"],angularCells=96,
               wallFaces=audit["nativeWall"],meshSHA256=key["meshSHA256"],bodySHA256=key["bodySHA256"],
               nodeMapSHA256=key["nodeMapSHA256"],meshAuditSHA256=key["mappingAuditSHA256"],reusedGrid=reuse,
               preview=preview))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("action",choices=["prepare"])
    p.add_argument("out",type=Path);p.add_argument("--rebuild",action="store_true")
    a=p.parse_args();out=a.out.resolve();out.relative_to(ROOT/"outputs");prepare(out,a.rebuild)
