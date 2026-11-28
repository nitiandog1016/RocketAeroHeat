"""仅用 Git 中的原生边和阶段小表，独立重算局部方程账与试算资格。"""
import csv
import json
import math
from pathlib import Path

DATA=Path(__file__).resolve().parent/"evidence"


def rows(name):
    with (DATA/name).open(encoding="ascii",newline="") as f:
        return list(csv.DictReader(f))


def verify():
    study=json.loads((DATA/"assembly.json").read_bytes())
    status=json.loads((DATA/"status.json").read_bytes())
    native={int(r["point"]):r for r in rows("matched_stages.csv")}
    supplied={int(r["point"]):r for r in rows("matched_nodes.csv")}
    ids={109*i+j for i in range(1,88) for j in (11,12,13,14)}
    assert set(supplied)==set(native)==ids
    edges=rows("matched_edges.csv");around={p:[] for p in ids}
    assert len({int(e["edge"]) for e in edges})==len(edges)==787
    for edge in edges:
        a,b=int(edge["left"]),int(edge["right"])
        if a in ids:around[a].append((b,1,edge))
        if b in ids:around[b].append((a,-1,edge))
    terms=("central","roeDissipation","viscous","axisSource","physicalTime")
    squared=[];dot={k:[] for k in terms};peak=(-1.,None);max_error=0.
    for p,connections in around.items():
        row={k:float(v) for k,v in native[p].items()};expect=supplied[p]
        assert {q for q,_,_ in connections}=={p-109,p+109,p-1,p+1}
        def flux(key):
            return math.fsum(sign*float(edge[key]) for _,sign,edge in connections)
        central,roe,total=[flux(k) for k in ("central","roe","total")]
        actual=(central,roe-central,total-roe,row["source"]-row["afterViscousCall"],
                (row["energy"]-row["previousEnergy"])*row["volume"]/row["dtND"])
        for key,value in zip(terms,actual):assert abs(value-float(expect[key]))<1e-12
        residual=math.fsum(actual)
        max_error=max(max_error,abs(residual+row["rhs"]))
        assert abs(total-row["transport"])<1e-12 and row["start"]==0.
        assert row["afterViscousCall"]==row["transport"] and row["boundary"]==row["time"]
        assert abs(residual-float(expect["physicalResidual"]))<1e-12
        peak=max(peak,(abs(row["rhs"]),p))
        if p%109 in (12,13):
            squared.append(residual**2)
            for key,value in zip(terms,actual):dot[key].append(value*residual)
    assert peak[1]==status["peakPoint"]==3828 and max_error<1e-12
    for key in terms:
        ratio=math.fsum(dot[key])/math.fsum(squared)
        assert abs(ratio-study["matched"]["directionalShare"][key])<1e-9
    assert max(study["matchedLogDifference"].values())<.001 and study["prefix"]["maxPrefixLogDifference"]<1e-9
    assert max(study["toggle"]["maxStateDifference"].values())<1e-9
    assert max(study["toggle"]["maxLogHistoryDifference"].values())<1e-9
    for r in study["refinedSteps"].values():
        assert r["passed"]==(max(r["last20Max"].values())< -8 and r["inner"]>=20)
        assert r["originalFaces"]==status["faces"]==9216
    assert status["sameTimePairAvailable"] is False and status["wholeBodyPowerW"] is None
    print(f"{len(ids)} nodes, {len(edges)} edges, peak {peak[1]}, max signed closure error {max_error:.3g}.")


if __name__=="__main__":
    verify()
