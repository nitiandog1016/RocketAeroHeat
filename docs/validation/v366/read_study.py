"""复核原生物理步、旋转域收支和残差空间分布，把可证实的状态接入同形档案。"""
import argparse
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from physical_step import BASE, DOC, OLD, ORIGIN, RES, STEPS, inspect, numbers, read, save, sha, vtk_tokens
from read_wall import mapping
from energy_budget import budget

METHODS=(("time01","adaptive"),("time02","fixed"),("time03","newton"))
RHS=("Residual_Density","Residual_Momentum_x","Residual_Momentum_y","Residual_Energy")


def spatial_residual(a,h,condition):
    """先从节点重建原生日志的 RMS 和峰值，再定位主导节点及其所属子午线。"""
    rows,cols=condition["shape"]
    wall=np.array(condition["meridianM"])
    assert wall[12,0]==0 and wall[12,1]>0 and .0149<wall[13,0]<.0151
    report={}
    for key,column in zip(RES,RHS):
        v=a[column].reshape(rows,cols)
        square=v**2
        lognorm=float(np.log10(np.sqrt(square.mean())))
        peak=float(np.log10(np.max(abs(v))))
        native_max=key.replace("rms[","max[")
        # 原生 NK 同一轮的右端取值在有限差分/输出阶段会有末位差；不把它作为热通量。
        assert abs(lognorm-h[key][-1])<.001 and abs(peak-h[native_max][-1])<.001
        row,col=map(int,np.unravel_index(abs(v).argmax(),v.shape));idx=row*cols+col
        report[column]=dict(logRMS=lognorm,historyLogRMS=float(h[key][-1]),logMaximum=peak,
            historyLogMaximum=float(h[native_max][-1]),peakPointID=idx,layerZeroBased=row,columnZeroBased=col,
            xM=condition["lengthM"]-float(a["x"][idx]),rM=float(a["y"][idx]),
            baseRimAndFirstTailColumnShare=float(square[:,12:14].sum()/square.sum()),
            firstFourLayerShare=float(square[:4].sum()/square.sum()))
    return report


def read_case(folder,kind,method,seed,c,matrix,areas,scopes):
    dt,i=STEPS[kind]
    checked=inspect(folder,kind,method)
    rec=read(folder/"step.json")
    assert all(checked[k]==value for k,value in rec.items())
    assert rec["seedSHA256"]==sha(seed) and rec["meshSHA256"]==sha(OLD/"mesh.su2")
    a=numbers(folder/f"restart_{i:05d}.csv");h=numbers(folder/f"history_{i:05d}.csv")
    assert min(abs(rec["nativeCurTimeSeconds"]-t) for t in (dt,i*dt))<1e-12
    rec["savedFieldBudget"]=budget(numbers(seed),a,dt,c)
    # 独立的壁面 VTK 和体场前 109 个点在求解坐标下逐点匹配，再按原 STL 映射。
    xyz,f=vtk_tokens(folder/f"surface_{i:05d}.vtk")
    distance,ix=cKDTree(xyz[:,:2]).query(np.column_stack((a["x"][:109],a["y"][:109])))
    assert len(np.unique(ix))==109 and max(distance)<1e-6
    for k in ("Pressure","Heat_Flux","Temperature"):
        np.testing.assert_allclose(f[k][ix,0],a[k][:109],rtol=6e-6,atol=1e-9)
    np.testing.assert_allclose(f["Skin_Friction_Coefficient"][ix,0],a["Skin_Friction_Coefficient_x"][:109],rtol=6e-6,atol=1e-9)
    p0=c["freePressurePa"];rho0=c["densityKgM3"];v0=np.sqrt(p0/rho0)
    dp=.5*rho0*c["speedMps"]**2
    tau=dp*a["Skin_Friction_Coefficient_x"][:109]
    nodes=np.column_stack((p0*a["Pressure"][:109],a["Heat_Flux"][:109],tau))
    qfirst=a["Temperature"][109:218]*c["freeTemperatureK"]
    tnode=a["Temperature"][:109]*c["freeTemperatureK"]
    xy=np.column_stack((c["lengthM"]-a["x"],a["y"]))
    edge=np.diff(xy[:109],axis=0);n=np.column_stack((-edge[:,1],edge[:,0]))
    dual=.5*np.vstack((n[0],n[:-1]+n[1:],n[-1]));length=np.linalg.norm(dual,axis=1)
    dist=np.einsum("ij,ij->i",xy[109:218]-xy[:109],dual/length[:,None])
    k=1.716e-5*(tnode/273.15)**1.5*(273.15+110.4)/(tnode+110.4)*1.4*287.05287/.4/.72
    np.testing.assert_allclose(nodes[:,1],k*(qfirst-c["givenWallK"])/dist,rtol=2e-8,atol=2e-5)
    np.testing.assert_allclose(nodes[:,1]@length,h["HF"][-1],rtol=2e-8,atol=1e-5)
    rev_nodes=(tau<0)&(v0*a["Velocity_x"][109:218]<0)
    reverse=np.isclose(matrix@rev_nodes.astype(float),1.,atol=1e-12,rtol=0)
    faces=matrix@nodes
    regions={}
    for name,s in scopes.items():
        p,q,shear=faces[s].T
        regions[name]=dict(faces=len(q),areaM2=float(areas[s].sum()),
            candidatePowerW=float(areas[s]@q),qRangeWm2=[float(min(q)),float(max(q))],
            pRangePa=[float(min(p)),float(max(p))],shearRangePa=[float(min(shear)),float(max(shear))],
            fullyReverseFaces=int(reverse[s].sum()))
    rec["stepFrameRegions" if rec["innerPassed"] else "rejectedFrameRegions"]=regions
    if "Residual_Energy" in a:
        rec["spatialResidual"]=spatial_residual(a,h,c)
    return rec,a,faces,reverse


def compare_pair(full,half,area,scopes,dynamic):
    """只对来源与时刻均已核实的原面场计算差额；零参照保持未定义。"""
    a,rev_a=full;b,rev_b=half
    comparison={}
    for name,s in scopes.items():
        w=area[s];delta=a[s]-b[s]
        rms=np.sqrt(np.average(delta**2,axis=0,weights=w))
        scale=np.sqrt(np.average(b[s,1]**2,weights=w))
        power=w@b[s,1]
        worst=np.argmax(abs(delta),axis=0)
        comparison[name]=dict(heatRmsFraction=float(rms[1]/scale) if scale>0 else None,
            powerDifferenceFraction=float(abs(w@delta[:,1]/power)) if power else None,
            pressureRmsOverDynamic=float(rms[0]/dynamic),shearRmsPa=float(rms[2]),
            maxAbsDifference=np.max(abs(delta),axis=0).tolist(),
            worstOriginalFaceIDs=(np.arange(1,len(area)+1)[s][worst]).tolist(),
            reverseStatusDifferentFaces=int(np.count_nonzero(rev_a[s]!=rev_b[s])))
    return comparison


def export(out):
    model,_,area,matrix,_=mapping();c=read(OLD/"condition.json")
    scopes={"all":slice(None),"cylinder":slice(2016,6912),"contractionEnd":slice(2016,2112),"noseEnd":slice(6816,6912)}
    out.mkdir(parents=True,exist_ok=False)
    records={};sums={};pairs={}
    for directory,name in METHODS:
        branch=BASE/directory
        admitted={}
        for kind in ("full","half1","half2"):
            seed=ORIGIN
            if kind=="half2":
                if not (branch/kind).exists():
                    continue
                if "half1" not in admitted:
                    raise ValueError("第二半步缺少通过的第一半步")
                seed=branch/"half1/restart_00001.csv"
            key=name+"_"+kind
            rec,a,faces,reverse=read_case(branch/kind,kind,name,seed,c,matrix,area,scopes)
            if rec["innerPassed"]:
                admitted[kind]=(faces,reverse)
            records[key]=rec
            sums[key]={"iterations":rec["innerIterations"],"maxLast20":max(rec["last20Max"].values()),
                       "innerPassed":rec["innerPassed"]}
            if "spatialResidual" in rec:
                # 各列及各层独立平方和可用小型归档复算 RMS 与集中程度。
                z=np.stack([a[col].reshape(*c["shape"]) for col in RHS],axis=-1)
                np.savetxt(out/f"{kind}_columns.csv",(z*z).sum(axis=0),delimiter=",",comments="",header=",".join(RHS))
                np.savetxt(out/f"{kind}_layers.csv",(z*z).sum(axis=1),delimiter=",",comments="",header=",".join(RHS))
        if "full" in admitted and "half2" in admitted:
            assert records[name+"_full"]["elapsedFromSeedSeconds"]==records[name+"_half2"]["elapsedFromSeedSeconds"]
            pairs[name]=compare_pair(admitted["full"],admitted["half2"],area,scopes,.5*c["densityKgM3"]*c["speedMps"]**2)
    full,half=[records["newton_"+k]["spatialResidual"]["Residual_Energy"] for k in ("full","half1")]
    status=dict(surfaceSHA256=model["fileSHA256"],faces=len(area),areaM2=float(area.sum()),mach=.30,
        altitudeM=5000,wallGivenK=c["givenWallK"],originSHA256=sha(ORIGIN),requestedEndpointUs=50,
        requestedHalfStepUs=25,stepReceipts=sums,acceptedPhysicalSteps=sum(r["innerPassed"] for r in records.values()),sameTimePairAvailable=bool(pairs),
        sameTimeHeatRmsFraction=None,wholeBodyPowerW=None,physicalQualification=False,
        fullEnergyRimShare=full["baseRimAndFirstTailColumnShare"],halfEnergyRimShare=half["baseRimAndFirstTailColumnShare"])
    if pairs:
        status["sameTimePairs"]=pairs
        fractions=[r["all"]["heatRmsFraction"] for r in pairs.values()]
        status["sameTimeHeatRmsFraction"]=max(fractions) if all(v is not None for v in fractions) else None
    save(out/"runs.json",records);save(out/"status.json",status)
    print(status,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("out",type=Path)
    export(p.parse_args().out)
