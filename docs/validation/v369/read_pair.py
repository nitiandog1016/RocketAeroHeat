"""从原生步末读出可用的单步研究场，记录全表面覆盖与未取得的同刻误差。"""
import argparse
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
import numpy as np
from short_step import BASE, DOC, OLD, RES, inspect, numbers, read, save, sha
from read_wall import mapping

FIELDS = ["pressurePa","qIntoWallWm2","wallFluidK","firstFluidK","shearDownstreamPa"]
GROUP_SOURCE = DOC.parent/"v365/research/surface.json"


def archive(pair, out):
    """此归档绑定本次实得的一步收敛及整步拒绝；其他终态须按其自己的证据处理。"""
    pair = pair.resolve()
    pair.relative_to(BASE)
    records = {key:inspect(pair/key,key) for key in ("full","half1")}
    for key,r in records.items():
        assert r == read(pair/key/"step.json")
    assert not records["full"]["passed"] and records["half1"]["passed"]
    assert not (pair/"half2").exists()
    model, _, area, matrix, meridian = mapping()
    c = read(OLD/"condition.json")
    a = numbers(pair/"half1/restart_00001.csv")
    h = numbers(pair/"half1/history_00001.csv")
    n = len(meridian)
    xy = np.column_stack((c["lengthM"]-a["x"], a["y"]))
    np.testing.assert_allclose(xy[:n],meridian,rtol=0,atol=5e-15)
    # 利用真实近壁两个节点和给定温度独立复核壁面热流的量纲与符号。
    tangent = np.diff(xy[:n],axis=0)
    normal = np.column_stack((-tangent[:,1],tangent[:,0]))
    dual = .5*np.vstack((normal[0],normal[:-1]+normal[1:],normal[-1]))
    length = np.linalg.norm(dual,axis=1)
    gap = np.sum((xy[n:2*n]-xy[:n])*dual,axis=1)/length
    assert min(gap)>0 and (length>0).all()
    tw = a["Temperature"][:n]*c["freeTemperatureK"]
    near = a["Temperature"][n:2*n]*c["freeTemperatureK"]
    k = 1.716e-5*(tw/273.15)**1.5*(273.15+110.4)/(tw+110.4)*1.4*287.05287/.4/.72
    q = a["Heat_Flux"][:n]
    qcheck = k*(near-c["givenWallK"])/gap
    np.testing.assert_allclose(q,qcheck,atol=2e-5,rtol=2e-8)
    np.testing.assert_allclose(length@q,h["HF"][-1],rtol=2e-8,atol=1e-5)
    dynamic = .5*c["densityKgM3"]*c["speedMps"]**2
    nodes = np.column_stack((a["Pressure"][:n]*c["freePressurePa"],q,tw,near,
                             dynamic*a["Skin_Friction_Coefficient_x"][:n]))
    weights, groups = np.unique(matrix.toarray(),axis=0,return_inverse=True)
    np.testing.assert_array_equal(groups+1,read(GROUP_SOURCE)["faceGroup"])
    grouped = weights@nodes
    faces = matrix@nodes
    np.testing.assert_allclose(grouped[groups],faces,rtol=5e-13,atol=1e-10)
    # 每个三角面按实际面积积分，底面使用跨全部半径的面积权重。
    scopes = {r["name"]:slice(r["firstFace"]-1,r["lastFace"]) for r in model["regions"]}
    scopes.update(all=slice(None),contractionEnd=slice(2016,2112),noseEnd=slice(6816,6912))
    regions = {}
    for name,s in scopes.items():
        p,qf = faces[s,:2].T
        regions[name] = dict(faces=len(qf),areaM2=float(sum(area[s])),
                             researchIntoWallW=float(area[s]@qf),pressureRangePa=[float(min(p)),float(max(p))],
                             qIntoWallRangeWm2=[float(min(qf)),float(max(qf))])
    signed_sum = sum(regions[r["name"]]["researchIntoWallW"] for r in model["regions"])
    np.testing.assert_allclose(signed_sum,regions["all"]["researchIntoWallW"],atol=1e-9,rtol=0)
    s = dict(mach=.3,altitudeM=5000,wallGivenK=c["givenWallK"],surfaceSHA256=model["fileSHA256"],
             faces=len(area),areaM2=float(sum(area)),elapsedFromNumericalSeedUs=records["half1"]["nominalEndUs"],
             frameSHA256=records["half1"]["frameSHA256"],groupSourceSHA256=sha(GROUP_SOURCE),
             columns=FIELDS,groupValues=grouped.tolist(),regions=regions,fullStepUs=records["full"]["deltaUs"],
             halfStepWorst=max(records["half1"]["last20Max"].values()),fullStepWorst=max(records["full"]["last20Max"].values()),
             firstHalfStepConverged=True,sameTimePairAvailable=False,sameTimeHeatError=None,
             maxHeatReadbackDifferenceWm2=float(max(abs(q-qcheck))),nativeFluxClosureProven=False,
             wholeBodyPowerW=None,physicalQualification=False)
    out.mkdir(parents=True,exist_ok=False)
    save(out/"surface.json",s)
    save(out/"runs.json",records)
    with (out/"wall.csv").open("x",encoding="ascii",newline="\n") as f:
        np.savetxt(f,np.column_stack((meridian,gap,nodes)),delimiter=",",fmt="%.17g",comments="",
                   header=",".join(["xM","rM","normalGapM",*FIELDS]))
    with ZipFile(out/"histories.zip","x") as z:
        for key in records:
            z.writestr(ZipInfo(key+".csv"),(pair/key/"history_00001.csv").read_bytes(),compress_type=ZIP_DEFLATED)
    print(regions["all"],s["maxHeatReadbackDifferenceWm2"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("pair",type=Path)
    p.add_argument("out",type=Path)
    args=p.parse_args()
    archive(args.pair,args.out)
