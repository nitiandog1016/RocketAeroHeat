"""从同一初始动力场保持全壁 12 K，推进到引线分支实际迭代位置并独立读取原面。"""
import argparse
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from thermistorBead import read,sha
from thermistorCrossflowSurface import collect
from thermistorMeasuredCold import advance,audit


def run(wire,out):
    done=read(wire/"completed.json");spec=read(wire/"prepared/profile.json");base=Path(spec["base"])
    assert done["frame"]==sha(wire/done["stage"]/"restart.csv") and spec["baseFrame"]==sha(base/"restart.csv")
    out.mkdir(parents=True,exist_ok=False)
    start=base/"restart.csv";previous=0
    for w in sorted(wire.glob("steady[0-9][0-9]")):
        dest=out/w.name;steps=read(w/"native.json")["steps"]
        advance(base.parent,dest,start,steps)
        receipt=audit(base.parent,dest,start,previous,steps)
        assert receipt["steps"]==steps,"两条分支原生停止位置不同"
        start=dest/"restart.csv";previous=receipt["cumulativeSteps"]
        if w.name==done["stage"]:break
    assert receipt["allResidualsBelowMinus8"]
    collect(dest)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("wire",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();out=a.out.resolve();out.relative_to(ROOT/"outputs");run(a.wire.resolve(),out)
