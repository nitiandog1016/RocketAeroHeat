"""按共同三角面和历史光滑原面分别比较三次冻结工况，并发布带来源的工作台对照。"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import zipfile

import numpy as np
from scipy.spatial import cKDTree

from check_mesh import HERE,ROOT,shared_faces
from beadSurfaceResult import projection,state_and_reaction
from thermistorBead import A,B,D,AREA,DT,P0,RHO,U,read,sha
from thermistorLeadMesh import meridian
from readCM011TunnelSector import vtk_tokens


def read_csv(raw):
    return np.genfromtxt(io.BytesIO(raw),names=True,delimiter=",")


def integral(w,mask):
    return float(-np.sum(w["areaM2"][mask]*w["heatFluxWm2"][mask]))


def metrics(delta,area):
    return dict(weightedMean=float(np.average(delta,weights=area)),weightedRMS=float(np.sqrt(np.average(delta**2,weights=area))),
                minimum=float(delta.min()),maximum=float(delta.max()))


def project_to_original(folder,probe,old_body,old_area,smooth):
    """两种光滑解按同一历史面积和子三角积分规则读回，分别使用各自的真实母线角。"""
    spec=read(folder/"case.json");m=spec["meridianIntervals"]+1
    with np.load(folder/"nodes.npz") as z:
        xyz,ids=z["xyz"],z["ids"]
    profile=meridian(True)[0] if smooth else None
    angles=np.arctan2(profile[:,1]/B,-profile[:,0]/A) if smooth else None
    rays=np.linspace(0,np.pi/6,5)
    mat,nodes=projection(np.array(old_body["V"]),np.array(old_body["F"])-1,ids[:,:m],rays,angles)
    native,reaction,res,receipt=state_and_reaction(folder,probe,nodes,xyz)
    dual=mat.T@old_area
    assert np.all(dual>0)
    flux=mat@(reaction[nodes]/dual)
    assert abs(old_area@flux-reaction[nodes].sum())<1e-13
    pressure=P0+(RHO*U**2)*(mat@native[nodes,4])
    # 用同一原面的面积积分权重分别比较切向壁剪与首层速度；两分量均反向才计入。
    wall=ids[:,:m];near=ids[:,m:2*m]
    xyz_wall=xyz[wall[0]]
    theta=np.arctan2(xyz_wall[:,1]/B,-xyz_wall[:,0]/A)
    tangent=np.column_stack((-A*np.sin(theta),-B*np.cos(theta)))
    tangent/=np.linalg.norm(tangent,axis=1)[:,None]
    surf,fields=vtk_tokens(folder/"surface.vtk")
    assert sha(folder/"surface.vtk")==receipt["outputSHA256"]["surface.vtk"]
    dist,idx=cKDTree(surf).query(xyz[nodes]);assert dist.max()<2e-6 and len(set(idx))==len(nodes)==len(surf)
    cf=np.zeros(len(xyz));u=np.zeros(len(xyz))
    for j,phi in enumerate(rays):
        directions=np.column_stack((tangent[:,0],tangent[:,1]*np.cos(phi),tangent[:,1]*np.sin(phi)))
        cf[wall[j]]=np.sum(fields["Skin_Friction_Coefficient"][idx[np.searchsorted(nodes,wall[j])]]*directions,axis=1)
        u[wall[j]]=np.sum(native[near[j],5:8]*directions,axis=1)
    return pressure,flux,(mat@cf[nodes]<0)&(mat@u[nodes]<0),res


def compare(out):
    proof=HERE/"evidence"
    spec=read(proof/"inputs.json")
    folder=ROOT/spec["path"]
    assert all(sha(folder/name)==value for name,value in spec["inputs"].items())
    cases=ROOT/"cases/thermistor_b"
    with zipfile.ZipFile(cases/"leads3.zip") as z:
        wb=json.loads(z.read("body.json"));wired=read_csv(z.read("wall.csv"));wr=json.loads(z.read("summary.json"))
        assert hashlib.sha256(z.read("wall.csv")).hexdigest()==wr["wallSHA256"]
    with zipfile.ZipFile(cases/"surface.zip") as z:old_body=json.loads(z.read("body.json"))
    with zipfile.ZipFile(cases/"surface_7p5.zip") as z:
        history=read_csv(z.read("wall.csv"));old=json.loads(z.read("summary.json"))
        assert hashlib.sha256(z.read("wall.csv")).hexdigest()==old["wallTableSHA256"]
    new=read(folder/"surface/summary.json");smooth=read_csv((folder/"surface/wall.csv").read_bytes())
    sb=read(folder/"body.json")
    assert wr==read(HERE.parent/"v297/evidence/fine_surface.json")
    assert old==read(HERE.parent/"v283/evidence/refined.json")
    assert new["wallSHA256"]==sha(folder/"surface/wall.csv") and new["bodySHA256"]==sha(folder/"body.json")
    assert new["sourceFrameSHA256"]==sha(folder/"restart.csv")
    wp,sp=shared_faces(wb,sb).T
    assert np.array_equal(wired["areaM2"][wp],smooth["areaM2"][sp])
    wf=read_csv(gzip.decompress((HERE.parent/"v297/evidence/nearwall_faces.csv.gz").read_bytes()))
    sf=read_csv((folder/"surface/flow.csv").read_bytes())
    assert sha(folder/"surface/flow.csv")==read(folder/"surface/flow.json")["tableSHA256"]
    assert np.array_equal(wf["face_id"],wired["faceID"]+1) and np.array_equal(sf["face_id"],smooth["faceID"]+1)
    old_run=ROOT/"outputs/v283_validation/strong01"
    old_probe=old_run.parent/"probes01/reaction0"
    prior_energy=read(HERE.parent/"v283/evidence/energy.json")
    assert old["sourceFrameSHA256"]==sha(old_run/"restart.csv")==prior_energy["inputFrameSHA256"]
    assert sha(old_probe/"restart.csv")==old["sourceReactionSHA256"] and sha(old_run/"body.json")==old["bodySHA256"]
    a=history["areaM2"];body_v=np.array(old_body["V"])[np.array(old_body["F"])-1]
    np.testing.assert_allclose(np.linalg.norm(np.cross(body_v[:,1]-body_v[:,0],body_v[:,2]-body_v[:,0]),axis=1)/2,a,rtol=2e-13,atol=0)
    p_old,q_old,r_old,_=project_to_original(old_run,old_probe,old_body,a,False)
    p_new,q_new,r_new,_=project_to_original(folder,folder/"reaction0",old_body,a,True)
    np.testing.assert_allclose(p_old,history["pressurePa"],rtol=0,atol=1e-9)
    np.testing.assert_allclose(q_old,history["heatFluxWm2"],rtol=0,atol=1e-9)
    baseline=integral(history,slice(None));control=integral(smooth,slice(None));leads=wr["beadAndJointsW"]
    assert abs(-a@q_new-control)<1e-13 and abs(baseline-prior_energy["wallReleaseW"])<1e-13
    assert abs(control-new["energy"]["wallReleaseW"])<1e-13 and abs(new["energy"]["relativeDifference"])<.005
    assert abs(integral(wired,slice(None))-wr["energy"]["wallReleaseW"])<1e-13
    common_w=integral(wired,wp);common_s=integral(smooth,sp);caps=integral(smooth,smooth["regionID"]!=2)
    joints=integral(wired,np.isin(wired["regionID"],[3,5]))
    assert abs((common_w-common_s)+(joints-caps)-(leads-control))<1e-13
    qw=wired["heatFluxWm2"][wp];qs=smooth["heatFluxWm2"][sp];weights=smooth["areaM2"][sp]
    rev_w=wf["reverse_both"][wp].astype(bool);rev_s=sf["reverse_both"][sp].astype(bool)
    exp=read(HERE.parent/"v282/evidence/comparison.json")["experiment"]
    heads=dict(oldSmooth=baseline,pairedSmooth=control,wiredBeadAndJoints=leads)
    flows=dict(commonReverseWiredM2=float(weights[rev_w].sum()),commonReverseSmoothM2=float(weights[rev_s].sum()),
               commonReverseBothM2=float(weights[rev_w&rev_s].sum()),commonReversalChangedFaces=int(np.count_nonzero(rev_w!=rev_s)),
               commonReversalChangedAreaM2=float(weights[rev_w!=rev_s].sum()),oldGridReverseOnOriginalM2=float(a[r_old].sum()),
               newGridReverseOnOriginalM2=float(a[r_new].sum()),historicalReverseChangedAreaM2=float(a[r_old!=r_new].sum()),
               smoothNative=read(folder/"surface/flow.json"))
    result=dict(powerW=heads,meanHWm2K={k:v/(AREA*DT) for k,v in heads.items()},experimentHWm2K=exp["hWm2K"],
                experimentRelativeDifference={k:v/(AREA*DT*exp["hWm2K"])-1 for k,v in heads.items()},
                methodDeltaW=control-baseline,geometryDeltaW=leads-control,combinedDeltaW=leads-baseline,
                completeWiredReleaseW=wr["totalOutwardW"],commonAreaM2=float(weights.sum()),commonFaces=len(wp),
                commonWiredW=common_w,commonSmoothW=common_s,smoothCapsW=caps,wiredJointsW=joints,
                smoothCapsAreaM2=new["areaM2"]-float(weights.sum()),wiredJointsAreaM2=float(wired["areaM2"][np.isin(wired["regionID"],[3,5])].sum()),
                commonHeatDifferenceWm2=metrics(qw-qs,weights),commonPressureDifferencePa=metrics(wired["pressurePa"][wp]-smooth["pressurePa"][sp],weights),
                smoothMethodHeatDifferenceWm2=metrics(q_new-q_old,a),smoothMethodPressureDifferencePa=metrics(p_new-p_old,a),
                flow=flows,sourceFrames=dict(wired=wr["sourceFrameSHA256"],smooth=new["sourceFrameSHA256"],historical=old["sourceFrameSHA256"]),
                bodies=dict(wired=wr["bodySHA256"],smooth=new["bodySHA256"],historical=old["bodySHA256"]),
                primaryTables=dict(wired=wr["wallSHA256"],smooth=new["wallSHA256"],historical=old["wallTableSHA256"]))
    out.mkdir(exist_ok=False)
    np.savetxt(out/"shared.csv",np.column_stack((wp+1,sp+1)),fmt="%d",delimiter=",",comments="",header="wiredID,smoothID")
    np.savetxt(out/"historical.csv",np.column_stack((history["faceID"]+1,a,p_new,q_new,r_old,r_new)),
               fmt=["%d","%.17g","%.17g","%.17g","%d","%d"],delimiter=",",comments="",
               header="faceID,areaM2,pairedPressurePa,pairedHeatFluxWm2,originalReverse,pairedReverse")
    result["historicalTransferSHA256"]=sha(out/"historical.csv");result["sharedMapSHA256"]=sha(out/"shared.csv")
    with (out/"comparison.json").open("x",encoding="utf8") as f:json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
    with zipfile.ZipFile(out/"paired3.zip","x") as z:
        for p in (folder/"body.json",folder/"case.json",folder/"surface/wall.csv",folder/"surface/summary.json",
                  folder/"surface/flow.csv",out/"shared.csv",out/"historical.csv",out/"comparison.json"):
            z.writestr(zipfile.ZipInfo(p.name,date_time=(2026,1,1,0,0,0)),p.read_bytes(),compress_type=zipfile.ZIP_DEFLATED)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument("out",type=Path)
    out=cli.parse_args().out.resolve();out.relative_to(ROOT/"outputs/v298_validation")
    compare(out)
