"""从端点归档、独立角点面积、实际两物理层及释放边界行审核同刻结果。"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from build_time import ROOT,HERE,BASE,BASE_FRAME,sha,read
sys.path.insert(0,str(HERE.parent/"v331"))
from verify_direction import normals_and_weights,raw_check,unit
from thermistorWireSurface import residual
from shortSphere65 import conductivity

DATA=ROOT/"cases/short_sphere"


def table(z,name):return np.genfromtxt(io.BytesIO(z.read(name)),names=True,delimiter=",")


def verify(data,work=None,ui=None):
    pair=read(data/"supportedTime65.json");s=pair["condition"]
    assert pair["steadyFrame"]==BASE_FRAME and sha(DATA/"supportedAngularTransfer65.npz")==pair["transferSHA256"]
    with np.load(DATA/"supportedAngularTransfer65.npz") as m:grid=dict(m)
    area=grid["targetArea"];xyz=grid["nativeWallXYZD"];c=grid["nativeCorners"];sa=grid["nativeArea"]
    nsrc=len(sa);N=len(area);a,b=grid["target"],grid["source"]
    _,weights=normals_and_weights(xyz,c)
    B=csr_matrix(((weights*sa[:,None]).ravel(),(np.repeat(np.arange(nsrc),4),c.ravel())),shape=(nsrc,len(xyz)))
    lump=B.sum(axis=0).A1
    Q=csr_matrix((grid["sourceFraction"],(a,b)),shape=(N,nsrc))
    I=csr_matrix((grid["targetFraction"],(a,b)),shape=(N,nsrc))
    with ZipFile(DATA/"supported65.zip") as z:body=json.loads(z.read("geometry.json"))
    f=np.array(body["F"])-1;grid["f"]=f;p=np.array(body["V"])[f]
    grid["axis"]=np.any(np.linalg.norm(p[:,:,1:],axis=2)<1e-12,axis=1)&(grid["region"]==1)
    n=unit(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]));basis=np.cross(n,np.cross([-1.,0.,0.],n))
    lengths=np.linalg.norm(basis,axis=1);ok=lengths>1e-10;basis[ok]/=lengths[ok,None]
    def intensive(v):return I@((B@v)/sa)
    def close(x,y):np.testing.assert_allclose(x,y,atol=1e-10,rtol=1e-12)
    def rms(v):return np.sqrt(np.average(v*v,weights=area))
    proof=dict(faces=N,areaM2=area.sum(),elapsedS=pair["elapsedS"]);all_rows={}
    for branch,step in (("full",1),("half",2)):
        archive=data/("supportedTime65"+branch.title()+".zip")
        assert sha(archive)==pair[branch+"Archive"]
        with ZipFile(archive) as z:
            r=json.loads(z.read("summary.json"));w=table(z,"wall.csv");station=table(z,"stations.csv");d=table(z,"direction.csv")
            assert all(hashlib.sha256(z.read(k)).hexdigest()==v for k,v in r["files"].items())
        all_rows[branch]=w;e=r["energy"]
        assert r["frameSHA256"]==pair[branch+"Frame"]==e["frame"] and r["seedFrame"]==pair["seed"]
        assert e["old"]==(pair["seed"] if step==1 else pair["firstHalfFrame"]) and e["elapsedDU"]==.25
        close(e["elapsedS"],e["dt"]*step*s["diameterM"]/s["speedMps"]);close(e["elapsedS"],pair["elapsedS"])
        assert r["jointExperimentalAreaM2"]==0 and r["faces"]==32480
        np.testing.assert_array_equal(w["faceID"],np.arange(N)+1);np.testing.assert_array_equal(w["regionID"],grid["region"])
        np.testing.assert_array_equal(d["faceID"],w["faceID"]);close(w["areaM2"],area)
        close(-w["fluxInWm2"]@area,e["wallW"]);close(e["wallW"]-e["farW"]-e["storageW"],e["relativeBalance"]*e["wallW"])
        assert e["heatChange"]<.002 and abs(e["relativeBalance"])<.005 and max(e["window"])<-8 and e["rawMax"]<1e-8
        assert int(w["fluxInWm2"].argmin())+1==r["peakOutFace"]
        for region in r["regions"]:
            active=grid["region"]==region["id"];close(-w["fluxInWm2"][active]@area[active],region["outW"])
        close(station["computedNu"],pair[branch+"Nu"])
        close(area[d["status"]==-1].sum(),pair[branch+"Reverse"]["areaM2"])
        if work:
            stage=work/branch/f"step{step}";fname=f"{step:05d}"
            state=np.loadtxt(stage/("restart_"+fname+".csv"),delimiter=",",skiprows=1,usecols=range(10))
            prev=np.loadtxt(stage/f"seed_{step-1:05d}.csv",skiprows=1,delimiter=",",usecols=9)
            tree=cKDTree(state[:,1:4]);dist,wet=tree.query(xyz);df,first=tree.query(grid["firstXYZD"])
            assert max(dist)<1e-12 and max(df)<1e-12 and len(set(wet))==len(xyz)
            edges=np.fromfile(stage/"edges.bin","<f8").reshape(-1,6);ei,ej=edges[:,:2].astype(int).T
            delta=edges[:,4]*(1-edges[:,5])
            load=(residual(stage/"release"/("restart_"+fname+".csv"))+np.bincount(ei,delta,minlength=len(state))-np.bincount(ej,delta,minlength=len(state)))[wet]*e["scale"]
            gauge=intensive(s["rho"]*s["speedMps"]**2*state[wet,4])
            q=(Q@(B@(load/lump)))/area
            expect=dict(fluxInWm2=q,gaugePa=gauge,absolutePa=s["airPa"]+gauge,
                        wallInputK=intensive(s["airK"]+s["scalarUnitK"]*state[wet,9]),firstK=intensive(s["airK"]+s["scalarUnitK"]*state[first,9]))
            errors={k:max(abs(x-w[k])) for k,x in expect.items()};assert max(errors.values())<1e-8
            node=np.fromfile(stage/"nodes.bin","<f8").reshape(-1,18);node=node[np.argsort(node[:,0])]
            storage=node[:,5]@(state[:,9]-prev)*e["scale"]/e["dt"];close(storage,e["storageW"])
            for st in station:
                ids=np.searchsorted(wet,[int(st["leftNode"]),int(st["rightNode"])]);mix=st["rightWeight"]
                heat=(1-mix)*load[ids[0]]/lump[ids[0]]+mix*load[ids[1]]/lump[ids[1]]
                close(-heat*s["diameterM"]/((st["wallK"]-s["airK"])*conductivity(st["wallK"])),st["computedNu"])
            mapped={"angular"+k[0].upper()+k[1:]:d[k] for k in d.dtype.names}
            direction=raw_check("angular",dict(r["direction"],frame=e["frame"],vtk=e["vtk"]),grid,n,basis,ok,mapped,stage,128,fname)
            steady=BASE/"run01/steady05/restart_000199.csv";assert sha(steady)==BASE_FRAME
            base_u=np.loadtxt(steady,delimiter=",",skiprows=1,usecols=(5,6,7))
            xx,yy,zz=state[:,1:4].T
            wake=((xx+1.15)/.55)**2+(yy*yy+zz*zz)/.45**2<1
            difference=np.linalg.norm(state[wake,5:8]-base_u[wake],axis=1)
            proof[branch]=dict(faceErrors=errors,directionErrors=direction,storageW=storage,
                               maxVelocityChangeOverUInSeedSupport=max(difference),
                               volumeRMSVelocityChangeOverUInSeedSupport=np.sqrt(np.average(difference**2,weights=node[wake,5])))
        if ui and (ui/branch).exists():
            folder=ui/branch;page=json.loads((folder/"workbench/case-data.js").read_text("utf8").split("=",1)[1].strip().rstrip(";"))
            assert page["evidenceIdentity"]["sourceFrame"]["sha256"]==e["frame"] and page["evidenceIdentity"]["sourceFrame"]["previousSHA256"]==e["old"]
            close(page["resultCondition"]["elapsedPhysicalTimeS"],e["elapsedS"])
            fields={v["key"]:v for v in page["fields"]}
            for k,x in (("heatFlux",w["fluxInWm2"]/1000),("pressure",w["absolutePa"]/1000),("pressureGauge",w["gaugePa"]),
                        ("nearWallCellTemperature",w["firstK"]),("wallInput",w["wallInputK"])):
                close(fields[k]["values"],x)
            for col in ("shearPa","speedMps","shearMagnitudePa","speedMagnitudeMps","status"):
                x=d[col].copy()
                if col=="status":x[x==0]=np.nan
                close(np.array(fields["supportTime"+branch+col[0].upper()+col[1:]]["values"],float),x)
            assert all(all(v is None for v in fields[k]["values"]) for k in ("wallSolved","edgeTemperature","recoveryTemperature"))
            with (folder/"section.csv").open(encoding="utf-8-sig",newline="") as handle:cut=list(csv.DictReader(handle))
            assert cut
            for row in cut:
                assert row["field_key"]=="heatFlux" and row["unit"]=="kW/m²"
                close(float(row["value"]),w["fluxInWm2"][int(row["face_id"])-1]/1000)
            proof[branch+"SectionRows"]=len(cut)
    x,y=all_rows["full"],all_rows["half"];a=pair["stepComparison"]
    close(rms(x["gaugePa"]-y["gaugePa"]),a["pressureRMSPa"]);close(rms(x["fluxInWm2"]-y["fluxInWm2"]),a["heatRMSWm2"])
    assert a["pressureRMSOverDynamic"]<.01 and a["heatRMSOverPeak"]<.01 and a["heatMaxOverPeak"]<.05 and a["wallPowerRelativeDifference"]<.01
    return proof


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--data",type=Path,default=DATA);p.add_argument("--run",type=Path);p.add_argument("--ui",type=Path);p.add_argument("--out",type=Path)
    a=p.parse_args();proof=verify(a.data,a.run,a.ui)
    if a.out:
        with a.out.open("x",encoding="utf8",newline="\n") as out:json.dump(proof,out,indent=2);out.write("\n")
    print(proof)
