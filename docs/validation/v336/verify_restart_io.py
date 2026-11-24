"""在相同节点上逐值比较原生观察到的温差和实际磁盘保存值，定位文本保存的差别。"""
import argparse
import json
from pathlib import Path
import numpy as np


def compare(steady):
    result={}
    for mark in ("early","late"):
        folder=steady/mark
        nodes=np.fromfile(folder/"nodes.bin","<f8").reshape(-1,16)
        nodes=nodes[np.argsort(nodes[:,0])]
        np.testing.assert_array_equal(nodes[:,0],np.arange(len(nodes)))
        for file,col in (("initial_flow.csv",1),("restart.csv",3)):
            values=np.loadtxt(folder/file,delimiter=",",skiprows=1,usecols=9)
            mismatch=values!=nodes[:,col]
            idx=np.flatnonzero(mismatch)
            data=dict(nodes=len(values),mismatches=len(idx),maxAbsDifference=float(max(abs(values-nodes[:,col]))))
            if len(idx):
                i=idx[0]
                data.update(firstNode=int(i),fromFile=values[i].hex(),fromNative=nodes[i,col].hex())
            result[mark+"/"+file]=data
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original",type=Path)
    parser.add_argument("corrected",type=Path)
    parser.add_argument("report",type=Path)
    a=parser.parse_args()
    result={k:compare(v) for k,v in (("original",a.original),("corrected",a.corrected))}
    assert result["original"]["early/initial_flow.csv"]["mismatches"]==8351
    assert not any(r["mismatches"] for r in result["corrected"].values())
    with a.report.open("x",encoding="utf8") as f:
        json.dump(result,f,indent=2);f.write("\n")
