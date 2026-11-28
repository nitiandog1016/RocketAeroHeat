"""把两网的同刻压力、热流与原生守恒项归到相同的原面和父网径向分带。"""
import argparse
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
import numpy as np
from outer_observer import BASE,DOC,SHAPE,read,save,sha,terminal,numbers
from outer_step import COARSE,COARSE_FRAME,OLD,LAYERS,RES,CFG_HASH,residual_window,mesh_points
from read_ledger import audit_state,save_csv
from read_pair import FIELDS,GROUP_SOURCE,wall_field
from read_wall import mapping
from read_study import compare_pair
from energy_budget import budget
from far_readback import calculate

ARCHIVE=DOC.parent/"v369/evidence"
BASE_BALANCE=DOC.parent/"v370/evidence"


def check_receipts(step,manifest):
    r=read(step/"step.json")
    proc=terminal(step)
    assert all(sha(Path(p))==h for p,h in proc["inputSHA256"].items())
    hashes={"frameSHA256":"restart_00001.csv","meshSHA256":"mesh.su2","seedSHA256":"seed_00000.csv",
            "configSHA256":"case.cfg","historySHA256":"history_00001.csv","wallSHA256":"surface_00001.vtk"}
    assert all(sha(step/file)==r[k] for k,file in hashes.items())
    assert manifest["fineSeedSHA256"]==r["seedSHA256"] and manifest["fineMeshSHA256"]==r["meshSHA256"]
    assert r["configSHA256"]==CFG_HASH and r["elapsedUs"]==6.25
    h=numbers(step/"history_00001.csv")
    w=residual_window(np.column_stack([h[k] for k in RES]))
    assert w["passed"] and all(r[k]==v for k,v in w.items())
    assert sha(COARSE/"restart_00001.csv")==COARSE_FRAME and sha(COARSE/"case.cfg")==CFG_HASH
    return r,h


def reference_flux(a,c):
    """用细网本身的远界状态重算特征单侧能流，供原生边界状态对齐。"""
    end=slice(-109,None)
    scale=np.sqrt(c["freePressurePa"]/c["densityKgM3"])
    xy=np.column_stack((a["x"][end],a["y"][end]))
    fields=np.column_stack((a["Density"][end]*c["densityKgM3"],a["Velocity_x"][end]*scale,
                           a["Velocity_y"][end]*scale,a["Pressure"][end]*c["freePressurePa"],a["Temperature"][end]*c["freeTemperatureK"]))
    return calculate(xy,fields,c)["characteristicNodal"]


def outer_layers(coarse_file,fine_file,href):
    """按父网线性基函数把新增中点的带符号项各半归入两端；保全域积分。"""
    def rows(path):
        a=np.genfromtxt(path,delimiter=",",names=True)
        return np.column_stack((a["geometryMinusUniform_kg_s"],a["geometryMinusUniform_W"]))
    before=rows(coarse_file);fine=rows(fine_file)
    assert before.shape==(89,2) and fine.shape==(91,2)
    restricted=np.zeros_like(before)
    low=LAYERS.astype(int);frac=(LAYERS-low)[:,None]
    np.add.at(restricted,low,(1-frac)*fine)
    np.add.at(restricted,np.ceil(LAYERS).astype(int),frac*fine)
    np.testing.assert_allclose(restricted.sum(axis=0),fine.sum(axis=0),atol=1e-8,rtol=0)
    def values(a):
        m,e=a[-2:].sum(axis=0)
        return dict(massKgS=float(m),energyW=float(e),excessW=float(e-href*m))
    return dict(coarse=values(before),fine=values(restricted))


def export(grid,step,replay,out):
    prep=read(grid/"prepared.json")
    record,h=check_receipts(step,prep)
    archive=read(ARCHIVE/"surface.json")
    model,_,area,matrix,meridian=mapping()
    assert model["fileSHA256"]==archive["surfaceSHA256"] and len(area)==9216
    c=read(OLD/"condition.json")
    a=numbers(step/"restart_00001.csv")
    nodes,gap,flux_error=wall_field(a,h,c,meridian)
    weights,group=np.unique(matrix.toarray(),axis=0,return_inverse=True)
    np.testing.assert_array_equal(group+1,read(GROUP_SOURCE)["faceGroup"])
    grouped=weights@nodes
    fine=matrix@nodes
    np.testing.assert_allclose(grouped[group],fine,rtol=1e-11,atol=1e-9)
    original=np.array(archive["groupValues"])[group]
    coarse_frame=numbers(COARSE/"restart_00001.csv")
    coarse_nodes,_,_=wall_field(coarse_frame,numbers(COARSE/"history_00001.csv"),c,meridian)
    np.testing.assert_allclose(matrix@coarse_nodes,original,rtol=1e-11,atol=1e-9)
    def reversed_faces(wall,frame):
        inward=(wall[:,4]<0)&(frame["Velocity_x"][109:218]<0)
        return np.isclose(matrix@inward.astype(float),1,rtol=0,atol=1e-12)
    scopes={k["name"]:slice(k["firstFace"]-1,k["lastFace"]) for k in model["regions"]}
    scopes.update(all=slice(None),contractionEnd=slice(2016,2112),noseEnd=slice(6816,6912))
    dyn=.5*c["densityKgM3"]*c["speedMps"]**2
    comparison=compare_pair((original[:,[0,1,4]],reversed_faces(coarse_nodes,coarse_frame)),
                            (fine[:,[0,1,4]],reversed_faces(nodes,a)),area,scopes,dyn)
    regions={}
    for key,s in scopes.items():
        v=comparison[key]
        w=area[s]
        p,q=fine[s,:2].T
        peak=np.flatnonzero(np.isclose(q,q.min(),rtol=1e-12,atol=1e-9))
        regions[key]=dict(faces=len(q),areaM2=float(sum(w)),researchIntoWallW=float(w@q),
                           coarseIntoWallW=archive["regions"][key]["researchIntoWallW"],
                           qRangeWm2=[float(q.min()),float(q.max())],pressureRangePa=[float(p.min()),float(p.max())],
                           peakOutwardFirstFace=int(np.arange(1,1+len(area))[s][peak[0]]),peakOutwardFaces=len(peak),
                           **v)
    surface=dict(frameSHA256=record["frameSHA256"],surfaceSHA256=model["fileSHA256"],faces=len(area),mach=.30,
                 altitudeM=5000,wallGivenK=c["givenWallK"],elapsedFromNumericalSeedUs=6.25,regions=regions)
    # 原初场以 17 位有效数写出；SU2 原生 ASCII 终场以小数点后 15 位的科学计数法写出。
    # 逐值核对序列化结果，储能积分的两个时间层共用实际加载的同一张静止网格。
    seed=numbers(step/"seed_00000.csv")
    xy=mesh_points(step/"mesh.su2",SHAPE[0]).reshape(-1,2)
    native_xy=np.fromfile(replay/"state.bin",dtype="<f8").reshape(-1,11)[:,:2]
    np.testing.assert_array_equal(xy,native_xy)
    for k,label in enumerate(("x","y")):
        np.testing.assert_array_equal(seed[label],xy[:,k])
        np.testing.assert_array_equal(a[label],[float(f"{v:.15e}") for v in xy[:,k]])
    saved=budget(seed,dict(a,x=xy[:,0],y=xy[:,1]),6.25e-6,dict(c,shape=list(SHAPE)))
    out.mkdir(parents=True,exist_ok=False)
    ledger=out/"native"
    audit_state(replay,ledger,step,surface,saved,reference_flux(a,c))
    old=read(BASE_BALANCE/"balance.json");new=read(ledger/"balance.json")
    def physical(s):
        mass,energy=s["nativeBoundaryPlusStorage"]
        return dict(massKgS=mass,energyW=energy,excessW=s["nativeExcessImbalanceW"])
    native=dict(coarse=physical(old),fine=physical(new))
    surface.update(coarseFrameSHA256=COARSE_FRAME,groupSourceSHA256=sha(GROUP_SOURCE),columns=FIELDS,
                   fineGroupValues=grouped.tolist(),areaM2=float(area.sum()),dynamicPressurePa=dyn,heatReadbackErrorWm2=flux_error,
                   nativeBalance=native,coarseNativeSHA256=sha(BASE_BALANCE/"balance.json"),fineNativeSHA256=sha(ledger/"balance.json"),
                   originalOuterBand=outer_layers(BASE_BALANCE/"layers.csv",ledger/"layers.csv",new["referenceEnthalpyJkg"]),
                   thermalSpatialComparisonPassed=all(v["heatRmsFraction"]<.02 and v["powerDifferenceFraction"]<.01 and
                       v["pressureRmsOverDynamic"]<.01 for v in comparison.values()),wholeBodyPowerW=None,physicalQualification=False)
    save(out/"comparison.json",surface);save(out/"step.json",record);save(out/"grid.json",prep)
    save_csv(out/"wall.csv",["xM","rM","normalGapM",*FIELDS],np.column_stack((meridian,gap,nodes)))
    with ZipFile(out/"history.zip","x") as z:
        z.writestr(ZipInfo("history_00001.csv"),(step/"history_00001.csv").read_bytes(),compress_type=ZIP_DEFLATED)
    print(native,regions["all"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for k in ("grid","step","replay","out"):p.add_argument(k,type=Path)
    a=p.parse_args()
    for path in (a.grid,a.step,a.replay):path.resolve().relative_to(BASE)
    export(a.grid.resolve(),a.step.resolve(),a.replay.resolve(),a.out.resolve())
