"""核对同次原生组装的所有边符号、隐式右端、官方基态及局部加密试算。"""
import argparse
from pathlib import Path
import numpy as np
from match_state import PARENT,BASE,DOC,RES,numbers,read,save,sha,terminal

STAGES=("start","transport","afterViscousCall","source","time","boundary","rhs")
TERMS=("central","roeDissipation","viscous","axisSource","physicalTime")
IDS=np.array([109*row+col for row in range(1,88) for col in range(11,15)])


def table(path):
    a=np.genfromtxt(path,delimiter=",",names=True)
    return np.sort(a,order=a.dtype.names[0])


def reduce_edges(edges,values):
    """沿原生边方向向左节点加、向右节点减，保留相邻对照列的全部外侧邻边。"""
    a=np.zeros(9701)
    np.add.at(a,edges["left"].astype(int),values)
    np.add.at(a,edges["right"].astype(int),-values)
    return a[IDS]


def assembly(folder,i):
    stage={name:table(folder/f"{i}_{name}.csv") for name in STAGES}
    start=stage["start"];val={k:v["energyResidual"] for k,v in stage.items()}
    for name,a in stage.items():
        np.testing.assert_array_equal(a["point"],IDS)
        for key in ("x","r","energy","previousEnergy","volume","dtND"):
            np.testing.assert_array_equal(a[key],start[key])
    edges=table(folder/f"{i}_edges.csv");ids=edges["edge"].astype(int)
    assert len(np.unique(ids))==len(edges)
    left,right=edges["left"].astype(int),edges["right"].astype(int)
    neighbours=[set() for _ in range(9701)]
    for a,b in zip(left,right):
        assert 0<=a<9701 and 0<=b<9701
        neighbours[a].add(b);neighbours[b].add(a)
    assert all(neighbours[n]=={n-1,n+1,n-109,n+109} for n in IDS)
    central=reduce_edges(edges,edges["central"]);roe=reduce_edges(edges,edges["roe"])
    total=reduce_edges(edges,edges["total"])
    components=np.column_stack((central,roe-central,total-roe,val["source"]-val["afterViscousCall"],val["time"]-val["source"]))
    verified_time=(start["energy"]-start["previousEnergy"])*start["volume"]/start["dtND"]
    rhs=val["rhs"];physical=-rhs
    errors={"edgeAssembly":max(abs(total-(val["transport"]-val["start"]))),
            "separateViscousCall":max(abs(val["afterViscousCall"]-val["transport"])),
            "timeReadback":max(abs(components[:,4]-verified_time)),
            "boundaryOnInterior":max(abs(val["boundary"]-val["time"])),
            "implicitSign":max(abs(rhs+val["boundary"])),
            "signedReconstruction":max(abs(rhs+components.sum(axis=1)))}
    assert max(errors.values())<1e-12 and not np.any(val["start"])
    core=np.isin(IDS%109,[12,13]);denom=physical[core]@physical[core]
    shares=components[core].T@physical[core]/denom
    np.testing.assert_allclose(shares.sum(),1.,atol=1e-10,rtol=0)
    peak=int(np.argmax(abs(physical)))
    around=abs(left-right)==1
    report=dict(nodes=len(IDS),edges=len(edges),checks={k:float(v) for k,v in errors.items()},
        peakPoint=int(IDS[peak]),peakPhysicalResidual=float(physical[peak]),peakImplicitRHS=float(rhs[peak]),
        peakComponents=dict(zip(TERMS,components[peak].tolist())),
        centralBetweenColumnsAtPeak=float(reduce_edges(edges,edges["central"]*around)[peak]),
        centralBetweenLayersAtPeak=float(reduce_edges(edges,edges["central"]*~around)[peak]),
        directionalShare=dict(zip(TERMS,shares.tolist())))
    return report,np.column_stack((IDS,start["x"],start["r"],components,physical,rhs)),edges


def compare_states(a,b):
    keys=("Density","Momentum_x","Momentum_y","Energy")
    x,y=[numbers(p/"restart_00001.csv") for p in (a,b)]
    assert np.array_equal(x["PointID"],y["PointID"])
    diff={k:float(np.max(abs(x[k]-y[k]))) for k in keys}
    h,j=[numbers(p/"history_00001.csv") for p in (a,b)]
    np.testing.assert_array_equal(h["Inner_Iter"],j["Inner_Iter"])
    return dict(maxStateDifference=diff,maxLogHistoryDifference={k:float(max(abs(h[k]-j[k]))) for k in RES})


def export(out):
    pair=BASE/"pair01";on=pair/"on";off=pair/"off";matched=BASE/"matched01";prefix=BASE/"official_prefix01"
    for folder in (on,off,matched,prefix):
        rc=terminal(folder)
        assert all(sha(Path(k))==v for k,v in rc["inputSHA256"].items())
    toggled=compare_states(on,off)
    assert max(toggled["maxStateDifference"].values())<1e-9 and max(toggled["maxLogHistoryDifference"].values())<1e-9
    for folder in (on,off):
        r=read(folder/"replay.json")
        assert all(sha(folder/k)==v for k,v in r["outputs"].items())
    local={str(i):assembly(on,i)[0] for i in (0,1,249)}
    current=read(prefix/"prefix.json");single=read(matched/"sample.json")
    for r,folder in ((current,prefix),(single,matched)):assert sha(folder/"restart_00001.csv")==r["frameSHA256"]
    assert current["maxPrefixLogDifference"]<1e-9 and max(single["logDifferences"].values())<.001
    binary=np.fromfile(matched/"bound_state.bin",dtype="<f8").reshape(9701,8)
    for k,file in enumerate((prefix/"restart_00001.csv",PARENT/"seed_00000.csv")):
        state=numbers(file)
        np.testing.assert_array_equal(binary[:,4*k:4*k+4],np.column_stack([state[j] for j in ("Density","Momentum_x","Momentum_y","Energy")]))
    summary,nodes,edges=assembly(matched,0)
    start=table(matched/"0_start.csv")
    np.testing.assert_array_equal(binary[IDS,3],start["energy"])
    np.testing.assert_array_equal(binary[IDS,7],start["previousEnergy"])
    official=numbers(PARENT/"restart_00001.csv")["Residual_Energy"][IDS]
    saved_difference=float(max(abs(official-nodes[:,-2])))
    initial=read(BASE/"refined_mesh01/prepared.json")
    assert initial["initialSHA256"]==sha(BASE/"refined_mesh01/initial.csv") and max(initial["relativeIntegralDifferences"])<1e-12
    refined={}
    for kind in ("full","half1","half2"):
        folder=BASE/"refined01"/kind
        if kind=="half2" and not folder.exists():continue
        rc=terminal(folder);assert all(sha(Path(k))==v for k,v in rc["inputSHA256"].items())
        r=read(folder/"step.json")
        i=2 if kind=="half2" else 1
        assert r["meshSHA256"]==initial["meshSHA256"] and r["frameSHA256"]==sha(folder/f"restart_{i:05d}.csv")
        h=numbers(folder/f"history_{i:05d}.csv")
        assert all(max(h[k][-20:])==r["last20Max"][k] for k in RES)
        refined[kind]=r
    status=dict(mach=.3,altitudeM=5000,wallGivenK=300,faces=9216,surfaceSHA256=read(DOC.parent/"v366/research/status.json")["surfaceSHA256"],
        modelXM=float(read(PARENT.parents[2]/"v365_validation/warm2/condition.json")["lengthM"]-nodes[np.argmax(abs(nodes[:,-1])),1]),
        peakPoint=summary["peakPoint"],peakCentral=summary["peakComponents"]["central"],
        peakPhysicalResidual=summary["peakPhysicalResidual"],reconstructionMax=summary["checks"]["signedReconstruction"],
        refinedFullPassed=refined["full"]["passed"],refinedHalfPassed=refined["half1"]["passed"],
        sameTimePairAvailable="half2" in refined and refined["half2"]["passed"] and refined["full"]["passed"],
        wholeBodyPowerW=None)
    out.mkdir(parents=True,exist_ok=False)
    save(out/"status.json",status)
    save(out/"assembly.json",dict(matched=summary,prefix=current,matchedLogDifference=single["logDifferences"],
        officialSavedFieldVsBasePhysicalResidualMax=saved_difference,localStages=local,toggle=toggled,localVsOfficial=compare_states(on,PARENT),
        seedSHA256=single["seedSHA256"],refinedMesh=initial,refinedSteps=refined))
    with (out/"matched_nodes.csv").open("x",encoding="ascii",newline="\n") as stream:
        np.savetxt(stream,nodes,delimiter=",",comments="",fmt=["%d"]+["%.17g"]*9,
            header=",".join(("point","solverX","radius",*TERMS,"physicalResidual","implicitRHS")))
    (out/"matched_edges.csv").write_bytes((matched/"0_edges.csv").read_bytes())
    state_cols=("energy","previousEnergy","volume","dtND")
    raw=np.column_stack([start["point"],*[start[k] for k in state_cols],
                         *[table(matched/f"0_{k}.csv")["energyResidual"] for k in STAGES]])
    with (out/"matched_stages.csv").open("x",encoding="ascii",newline="\n") as stream:
        np.savetxt(stream,raw,delimiter=",",comments="",
            header=",".join(("point",*state_cols,*STAGES)),fmt=["%d"]+["%.17g"]*(raw.shape[1]-1))
    print(status,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("out",type=Path);export(p.parse_args().out)
