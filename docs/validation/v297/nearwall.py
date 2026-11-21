"""配对两套网格的近壁流向并保存细网逐面诊断。"""
import gzip
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2];sys.path.insert(0,str(ROOT/"scripts"))
from thermistorLeadNearWall import assess,write_faces

if __name__=="__main__":
    directory=ROOT/"outputs/v297_validation";report={};masks=[]
    for label in ("base","fine"):
        summary,rows=assess(directory/(label+"01"))
        report[label]=summary;masks.append(rows[:,-1].astype(bool))
    report["changedReverseClassificationFaces"]=int(np.count_nonzero(masks[0]!=masks[1]))
    with gzip.open(HERE/"evidence/nearwall_faces.csv.gz","xt",encoding="ascii") as target:
        write_faces(target,rows)
    with (HERE/"evidence/nearwall.json").open("x",encoding="utf8") as f:json.dump(report,f,ensure_ascii=False,indent=2)
