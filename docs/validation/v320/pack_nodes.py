"""从已验收的 B 珠原生终帧取出紧凑的壁面及首层节点场，保留输入指纹。"""
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/"scripts"))
from thermistorSubfaces import PACK,extract,read,sha


if __name__=="__main__":
    loc=read(ROOT/"outputs/v319_validation/native01/ui.json");prep=Path(loc["prepared"])
    manifest=read(prep/"prepared.json");assert sha(prep/"nodes.npz")==manifest["nodeMapSHA256"]
    cases=ROOT/"cases/thermistor_b";wp=read(cases/"windPoint12.json")
    job=Path(loc["native"]["directory"])/"run";cold=ROOT/"outputs/v316_validation/run01"
    long=read(cases/"longWire12.json");res=read(cases/"resistance12.json")
    sources=[(ROOT/"outputs/v315_validation/six01",wp["sixKFrame"]),
             (ROOT/"outputs/v315_validation/native01",wp["sourceFrameSHA256"]),
             (cold/read(cold/"completed.json")["stage"],read(job/"surface/summary.json")["baseFrame"]),
             (job/read(job/"completed.json")["stage"],long["frameSHA256"]),
             (ROOT/"outputs/v318_validation/check01/steady01",res["frameSHA256"])]
    PACK.mkdir(exist_ok=True);metadata={}
    with np.load(prep/"nodes.npz") as grid:
        for stage,frame in sources:
            src=stage if (stage/"surface/summary.json").is_file() else stage.parent
            s=read(src/"surface/summary.json");spec=stage/"case.json"
            if not spec.is_file():spec=stage.parent/"prepared/case.json"
            reaction=s.get("sourceReactionSHA256",s.get("reactionSHA256"))
            assert s["sourceFrameSHA256"]==frame and sha(spec)==s["caseSHA256"]
            assert s["meshSHA256"]==manifest["meshSHA256"] and s["nodeMapSHA256"]==manifest["nodeMapSHA256"]
            vals=extract(stage,frame,reaction,grid,read(spec))
            # 用冻结的原面表反查缩放、反力符号及压力参考；这些判断独立于子面。
            tab=np.genfromtxt(src/"surface/wall.csv",delimiter=",",names=True)
            wi=np.searchsorted(np.unique(grid["wallMap"]),grid["wallMap"])
            for key,col in (("gauge","pressureGaugePa"),("flux","heatFluxWm2")):
                cs=vals[key][wi[grid["overlapChild"]]]
                total=np.einsum("ij,ij->i",grid["shapeAreaM2"],cs)
                coarse=np.bincount(grid["overlapParent"],total)/grid["originalArea"]
                np.testing.assert_allclose(coarse,tab[col],rtol=1e-12,atol=1e-8)
            path=PACK/(frame+".npz");np.savez_compressed(path,**vals)
            metadata[frame]=dict(sha256=sha(path),reaction=reaction,mesh=s["meshSHA256"],nodeMap=s["nodeMapSHA256"],
                                 case=s["caseSHA256"],wall=s["wallSHA256"],wallNodes=len(vals["flux"]),firstNodes=len(vals["first"]))
    (PACK/"manifest.json").write_text(json.dumps(metadata,indent=2)+"\n",encoding="utf8")
