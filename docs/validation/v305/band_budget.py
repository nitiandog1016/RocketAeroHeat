"""用原生逐边热通量和同阶段节点残差求热敏珠各带的独立收支。正号朝带外。"""
import argparse
import json
from pathlib import Path

import numpy as np

from common_layer import BASE, ENDS, GRIDS, HERE, LABELS, REGIONS, SCALE, read, save, sha
from thermistorBead import D

TOL = 1e-11
N = len(LABELS)


def check_watts(left, right):
    np.testing.assert_allclose(left, right, atol=TOL, rtol=0)


def stage(sample, name, theta):
    a = np.loadtxt(sample/(name+".csv"), delimiter=",", skiprows=1)
    row = a[np.argsort(a[:,0])]
    np.testing.assert_array_equal(row[:,0], np.arange(len(theta)))
    np.testing.assert_array_equal(row[:,1], theta)
    return SCALE*row[:,2]


def transfers(i,j,groups,flux):
    # 行是流出的带，列是流入的带；最后一行和一列代表近壁层外侧。
    region = np.where(groups>=0, groups, N)
    left, right = region[i], region[j]
    sent = np.bincount(left*(N+1)+right, weights=flux, minlength=(N+1)**2).reshape(N+1,N+1)
    net = sent - sent.T
    check_watts(net + net.T, 0)
    return net


def account(level, sample, geom, info):
    native=GRIDS[level-2];key=info["grids"][level-2]
    rec=read(sample/"sample.json");step=read(sample/"step_process.json")
    assert rec["grid"]==level and rec["maxOfficialStepDifference"]<2e-7 and step["exitCode"]==0 and not (step["stopReason"] or step["livePids"])
    assert all(sha(sample/k)==v for k,v in rec["outputs"].items()) and sha(sample/"initial_flow.csv")==key["frameSHA256"]
    assert all(sha(Path(p))==h for p,h in rec["inputs"].items())
    assert sha(geom/f"grid{level}.npz")==key["mapSHA256"]
    with np.load(geom/f"grid{level}.npz") as z:
        canonical, groups, expected, faces=(z[k] for k in ("canonical","group","expectedEdges","faceBand"))
    with np.load(native/"nodes.npz") as z:
        wall=z["wallMap"];area=z["area"];wall_nodes=np.unique(wall)
    theta=np.loadtxt(sample/"initial_flow.csv",delimiter=",",skiprows=1)[:,-1]
    stages={k:stage(sample,k,theta) for k in ("start","internal","beforeFar","afterFar","beforeWall","afterWall","end")}
    start,interior,far0,far1,wall0,wall1,end=stages.values()
    chosen=groups>=0
    edge=np.loadtxt(sample/"edges.csv",delimiter=",",skiprows=1)
    assert np.isfinite(edge).all()
    i,j=edge[:,0].astype(int),edge[:,1].astype(int)
    assert np.all(i!=j) and min(i.min(),j.min())>=0 and max(i.max(),j.max())<len(theta)
    assert len(np.unique(np.sort(np.column_stack((i,j)),axis=1),axis=0))==len(i)
    contrib=edge[:,5:]*SCALE
    for col in (0,2,4):check_watts(contrib[:,col]+contrib[:,col+1],0)
    check_watts(contrib[:,:2]+contrib[:,2:4],contrib[:,4:6])
    assembled=np.bincount(i,contrib[:,4],len(theta))+np.bincount(j,contrib[:,5],len(theta))
    check_watts(assembled,interior-start)
    check_watts(start,0);check_watts(interior,far0);check_watts(far1,wall0);check_watts(wall1,end)
    check_watts(far1[chosen],far0[chosen])
    # 共同边的两端与有向对偶面积逐边配对，包含径向、沿带及带内的边。
    pick=chosen[i]|chosen[j]
    ends=np.column_stack((canonical[i[pick]],canonical[j[pick]]))
    assert ends.min()>=0
    direction=np.where(ends[:,0]<ends[:,1],1.,-1.)
    ends.sort(axis=1);order=np.lexsort((ends[:,1],ends[:,0]))
    np.testing.assert_array_equal(ends[order],expected)
    normals=(direction[:,None]*edge[pick,2:5])[order]
    norms=np.linalg.norm(normals,axis=1);assert np.all(norms>0)
    # 独立弱化强壁产生的反力、原生强壁装配和原面功率分别核对。
    native_wall=wall0-wall1
    reaction=np.genfromtxt(native/"reaction0/restart.csv",names=True,delimiter=",")
    np.testing.assert_array_equal(reaction["PointID"],np.arange(len(theta)))
    iswall=np.zeros(len(theta),bool);iswall[wall_nodes]=True
    check_watts(native_wall[~iswall],0)
    check_watts(native_wall[wall_nodes],-reaction["Residual_Species_0"][wall_nodes]*SCALE)
    summary=read(native/"surface/summary.json");w=np.genfromtxt(native/"surface/wall.csv",names=True,delimiter=",")
    assert sha(native/"surface/wall.csv")==key["faceTableSHA256"]
    np.testing.assert_array_equal(w["faceID"],np.arange(len(faces)))
    np.testing.assert_array_equal(w["regionID"],REGIONS[faces])
    np.testing.assert_array_equal(w["areaM2"],area)
    check_watts(native_wall.sum(),summary["totalOutwardW"]);check_watts((far1-far0).sum(),summary["energy"]["farOutW"])
    assert abs((native_wall.sum()-(far1-far0).sum())/native_wall.sum())<.005
    # 原面公用的壁节点按三角形面积/3分配。矩阵记录从节点组到原面带的全部转移。
    lumps=np.bincount(wall.ravel(),np.repeat(area/3,3),len(theta))
    part=native_wall[wall]*area[:,None]/(3*lumps[wall])
    face_power=-w["heatFluxWm2"]*area
    check_watts(part.sum(axis=1),face_power)
    transfer=np.bincount((groups[wall]*N+faces[:,None]).ravel(),part.ravel(),minlength=N*N).reshape(N,N)
    cvwall=np.bincount(groups[wall_nodes],native_wall[wall_nodes],minlength=N)
    surface=np.bincount(faces,face_power,minlength=N)
    check_watts(transfer.sum(axis=0),surface);check_watts(transfer.sum(axis=1),cvwall)
    for region,r in enumerate(summary["regions"],1):check_watts(surface[REGIONS==region].sum(),r["outwardW"])
    adv=transfers(i,j,groups,contrib[:,0]);diff=transfers(i,j,groups,contrib[:,2])
    residual=np.bincount(groups[chosen],end[chosen],minlength=N)
    outer_adv,outer_diff=adv[:N,N],diff[:N,N]
    side_adv,side_diff=adv[:N,:N].sum(axis=1),diff[:N,:N].sum(axis=1)
    check_watts(side_adv.sum(),0);check_watts(side_diff.sum(),0)
    # 非相邻母线组之间不共享控制体界面。
    disconnected=np.abs(np.arange(N)[:,None]-np.arange(N)[None,:])>1
    check_watts(adv[:N,:N][disconnected],0);check_watts(diff[:N,:N][disconnected],0)
    imbalance=outer_adv+outer_diff+side_adv+side_diff-cvwall-residual
    check_watts(imbalance,0)
    assert np.all(abs(residual)<.005*abs(cvwall))
    correction=surface-cvwall
    check_watts(outer_adv+outer_diff+side_adv+side_diff+correction-residual,surface)
    names=("outerConvW","outerDiffW","sideConvW","sideDiffW","residualW","controlWallW","mapCorrectionW","faceOutW")
    columns=(outer_adv,outer_diff,side_adv,side_diff,residual,cvwall,correction,surface)
    bands=[dict(zip(names,(float(col[k]) for col in columns))) for k in range(N)]
    notes=dict(frameSHA256=key["frameSHA256"],meshSHA256=key["meshSHA256"],wallSHA256=key["faceTableSHA256"],
               edgeFileSHA256=rec["outputs"]["edges.csv"],sampleSHA256=sha(sample/"sample.json"),edges=len(i),
               normalAreaRangeD2=[float(norms.min()),float(norms.max())],maxNodeAssemblyErrorW=float(abs(assembled-interior+start).max()),
               maxEdgeSplitErrorW=float(abs(contrib[:,:2]+contrib[:,2:4]-contrib[:,4:6]).max()),
               maxBandClosureErrorW=float(abs(imbalance).max()),maxBandResidualFraction=float(np.max(abs(residual/cvwall))),
               wallW=float(cvwall.sum()),farW=float((far1-far0).sum()),wallMappingW=transfer.tolist(),
               convectiveExchangeW=adv.tolist(),diffusiveExchangeW=diff.tolist(),bands=bands)
    return normals,faces,area,notes


def compare(geometry,samples,out):
    info=read(geometry/"geometry.json")
    data=[account(i,s,geometry,info) for i,s in enumerate(samples,2)]
    a,b=data
    norm_error=np.linalg.norm(a[0]-b[0],axis=1)
    rel_error=norm_error/np.linalg.norm(a[0],axis=1)
    assert rel_error.max()<1e-11
    np.testing.assert_array_equal(a[1],b[1]);np.testing.assert_array_equal(a[2],b[2])
    case=read(GRIDS[0]/"surface/summary.json")
    bands=[]
    for k,(medium,fine) in enumerate(zip(a[-1]["bands"],b[-1]["bands"])):
        select=a[1]==k
        delta={n:fine[n]-value for n,value in medium.items()}
        check_watts(delta["faceOutW"],sum(delta[n] for n in ("outerConvW","outerDiffW","sideConvW","sideDiffW","mapCorrectionW"))-delta["residualW"])
        bands.append(dict(id=k+1,name=LABELS[k],regionID=int(REGIONS[k]),faces=int(select.sum()),areaM2=float(a[2][select].sum()),
                          medium=medium,fine=fine,delta=delta))
    grids=[t[-1] for t in data]
    package=dict(schema="thermistor-b-band-heat-v1",mode="givenWall",wallK=305.15,airK=293.15,speedMps=3,
                 bodySHA256=case["bodySHA256"],caseSHA256=case["caseSHA256"],sourceGrids=[
                     {k:g[k] for k in ("meshSHA256","frameSHA256","wallSHA256","edgeFileSHA256","sampleSHA256")} for g in grids],
                 bands=bands,faceBand=(a[1]+1).tolist(),geometryAuditSHA256=sha(geometry/"geometry.json"))
    out.mkdir(exist_ok=False)
    # 用固定 UTF-8 字节输出工作台输入，跨平台换行规则不影响数据包指纹。
    with (out/"layerBands.json").open("xb") as stream:
        stream.write(json.dumps(package,ensure_ascii=False,separators=(",",":"),allow_nan=False).encode("utf8"))
    save(out/"budget.json",dict(layerBandsSHA256=sha(out/"layerBands.json"),grids=grids,bands=bands,
         matchedEdges=len(rel_error),maxNormalRelativeDifference=float(rel_error.max()),maxNormalAreaDifferenceM2=float(norm_error.max()*D**2)))
    print([(band["name"],round(band["delta"]["faceOutW"]*1e6,3)) for band in bands])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("geometry","grid2","grid3","output"):p.add_argument(name,type=Path)
    a=p.parse_args();paths=[v.resolve() for v in (a.geometry,a.grid2,a.grid3,a.output)]
    for path in paths:path.relative_to(BASE)
    compare(paths[0],paths[1:3],paths[3])
