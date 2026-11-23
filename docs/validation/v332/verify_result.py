"""用独立角点面积、原生释放反力与逐边修正核查两端点；另可只检查发布小包。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from sector_grid import ROOT, HERE, read, sha

sys.path[:0]=[str(HERE.parent/v) for v in ("v331","v330","v329")]
from verify_direction import normals_and_weights, raw_check, table as zipped
from thermistorWireSurface import residual
from shortSphere65 import conductivity

DATA=ROOT/"cases/short_sphere"


def check_receipts(path,thermal):
    with ZipFile(path) as z:
        m=json.loads(z.read("manifest.json"))
        for name,expected in m["included"].items():
            assert hashlib.sha256(z.read(name)).hexdigest()==expected
        conv=json.loads(z.read("run01/converged.json"))
        assert conv["frame"]==thermal["early"]["sourceFrame"] and len(conv["trend"])*200==thermal["updates"]
        for k in ("early","late"):
            wall=np.loadtxt(z.open(k+"_wall.csv"),delimiter=",");far=np.loadtxt(z.open(k+"_far.csv"),delimiter=",")
            assert len(np.unique(wall[:,0]))==len(wall)==11906 and len(np.unique(far[:,0]))==len(far)==2690
            np.testing.assert_array_equal(far[:,2],far[:,3])
            watts=thermal["scale"]*np.sum(-wall[:,1]-wall[:,2]+wall[:,3])
            transported=thermal["scale"]*(np.maximum(far[:,2],0)@far[:,1])
            assert abs(watts-thermal[k]["wallW"])<1e-10 and abs(transported-thermal[k]["scalarFarW"])<1e-10


def verify(data,work=None,ui=None):
    comparison=read(data/"supportedAngular65.json");pair=read(DATA/"supportedNearWall65.json")
    transfer=data/"supportedAngularTransfer65.npz";surface=data/"supportedAngularSurface65.zip";directions=data/"supportedAngularDirection65.zip"
    assert sha(surface)==comparison["fineArchiveSHA256"] and sha(transfer)==comparison["transferSHA256"]
    assert sha(directions)==comparison["directionArchiveSHA256"] and pair["fineFrame"]==comparison["coarseFrame"]
    with ZipFile(surface) as z:
        s=json.loads(z.read("summary.json"));t=s["thermal"];cfg=s["condition"]
        assert s["frameSHA256"]==comparison["fineFrame"]==t["early"]["sourceFrame"] and t["passed"]
    w=zipped(surface,"wall.csv");stations=zipped(surface,"stations.csv");d=zipped(directions,"direction.csv")
    base=zipped(DATA/"supportedFine65.zip","wall.csv");prior=zipped(DATA/"supportedDirection65.zip","direction.csv")
    assert s["caseSHA256"]==pair["caseSHA256"] and s["jointExperimentalAreaM2"]==0
    with np.load(transfer) as data_file:m=dict(data_file)
    with ZipFile(DATA/"supported65.zip") as z:body=json.loads(z.read("geometry.json"))
    area=m["targetArea"];np.testing.assert_array_equal(area,w["areaM2"])
    ids=np.arange(1,len(area)+1)
    for x in (w,base,prior,d):np.testing.assert_array_equal(x["faceID"],ids);np.testing.assert_array_equal(x["regionID"],m["region"])
    q=w["fluxInWm2"];oldq=base["fluxInWm2"];status=d["status"];old=prior["fineStatus"]
    close=lambda x,y: np.testing.assert_allclose(x,y,rtol=1e-12,atol=1e-12)
    close(-q@area,comparison["fineW"]);close(-oldq@area,comparison["coarseW"])
    close(comparison["fineStationRMSE"],np.mean((stations["reportedNu"]-stations["computedNu"])**2)**.5)
    np.testing.assert_array_equal(np.where(q<=.99*q.min())[0]+1,comparison["finePeak"]["nearFaces"])
    for zone in comparison["regions"]:
        mask=m["region"]==zone["id"]
        close([-q[mask]@area[mask],-oldq[mask]@area[mask],-q[mask].min()],[zone["fineW"],zone["coarseW"],zone["finePeakWm2"]])
    close(area[(status==-1)&(old==-1)].sum(),comparison["reverseOverlap"]["sharedM2"])
    close(area[(status==-1)|(old==-1)].sum(),comparison["reverseOverlap"]["unionM2"])
    for a,label in [(-1,"reverse"),(1,"forward"),(0,"unknown")]:
        mask=status==a;rec=comparison["direction"]["total"][label]
        close([mask.sum(),area[mask].sum(),-q[mask]@area[mask]],[rec["faces"],rec["areaM2"],rec["outW"]])
    out=dict(faces=len(q),areaM2=area.sum(),frame=s["frameSHA256"],publicLedger=True)
    if work is None:
        check_receipts(HERE/"evidence/receipts.zip",t)
        out["archivedEndpointBudgetsVerified"]=True
    if work:
        cold=comparison["geometry"]["coldSHA256"];end=read(work/"converged.json")
        for i in range(1,len(end["trend"])+1):
            stage=work/f"steady{i:02d}";r=read(stage/"record.json");proc=read(stage/"debugger_process.json")
            assert r["start"]==cold and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
            assert not(proc["stopReason"] or proc["supervisorError"] or proc["survivingProcessIDs"])
            assert "Exit Success (SU2_CFD)" in (stage/"debugger.log").read_text(errors="replace")
            cold=sha(stage/"restart_000199.csv");assert cold==r["end"]
        assert cold==s["frameSHA256"] and max(r["window"])<-8 and r["rawMax20"]<1e-8
        ref=m["nativeWallXYZD"];c=m["nativeCorners"];srcarea=m["nativeArea"]
        _,portion=normals_and_weights(ref,c);a,b=m["target"],m["source"]
        B=csr_matrix(((srcarea[:,None]*portion).ravel(),(np.repeat(np.arange(len(c)),4),c.ravel())),shape=(len(c),len(ref)))
        lump=B.sum(axis=0).A1
        Q=csr_matrix((m["sourceFraction"],(a,b)),shape=(len(q),len(c)))
        I=csr_matrix((m["targetFraction"],(a,b)),shape=Q.shape)
        def field(v):return I@((B@v)/srcarea)
        for label,rel in (("early","releaseEarly"),("late","releaseLate")):
            state=np.loadtxt(work/label/"initial_flow.csv",skiprows=1,delimiter=",")
            tree=cKDTree(state[:,1:4]);dist,wall=tree.query(ref);df,first=tree.query(m["firstXYZD"])
            assert max(dist)<1e-12 and max(df)<1e-12 and len(np.unique(wall))==len(ref)
            edge=np.fromfile(work/label/"edges.bin","<f8").reshape(-1,6)
            qty=edge[:,4]*(edge[:,5]-1);ei,ej=edge[:,:2].astype(int).T
            corr=np.bincount(ei,qty,minlength=len(state))-np.bincount(ej,qty,minlength=len(state))
            load=(residual(work/rel/"restart.csv")[wall]-corr[wall])*t["scale"]
            actual=Q@(B@(load/lump))/area
            close(-actual@area,t[label]["wallW"])
            if label=="late":
                out["maxHeatFluxDriftWm2"]=max(abs(actual-q));continue
            gauge=field(cfg["rho"]*cfg["speedMps"]**2*state[wall,4])
            expected={"fluxInWm2":actual,"absolutePa":cfg["airPa"]+gauge,"gaugePa":gauge,
                      "wallInputK":field(cfg["airK"]+cfg["scalarUnitK"]*state[wall,9]),
                      "firstK":field(cfg["airK"]+cfg["scalarUnitK"]*state[first,9])}
            errors={name:max(abs(v-w[name])) for name,v in expected.items()};assert max(errors.values())<1e-8
            out["maxFaceErrors"]=errors
            lut={k:i for i,k in enumerate(wall)}
            for row in stations:
                nodes=np.array([lut[int(row[k])] for k in ("leftNode","rightNode")]);r=row["rightWeight"]
                qst=np.array([1-r,r])@(load[nodes]/lump[nodes])
                nu=-qst*cfg["diameterM"]/((row["wallK"]-cfg["airK"])*conductivity(row["wallK"]))
                close(nu,row["computedNu"])
        f=np.array(body["F"])-1;m["f"]=f;p=np.array(body["V"])[f]
        n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        unit=np.cross(n,np.cross([-1.,0.,0.],n));length=np.linalg.norm(unit,axis=1);defined=length>1e-10
        unit[defined]/=length[defined,None];m["axis"]=np.any(np.linalg.norm(p[:,:,1:],axis=2)<1e-12,axis=1)&(m["region"]==1)
        mapped={"angular"+k[0].upper()+k[1:]:d[k] for k in d.dtype.names}
        item=dict(comparison["direction"],frame=s["frameSHA256"],vtk=comparison["sourceVTKSHA256"])
        out["directionMaxAbsError"]=raw_check("angular",item,m,n,unit,defined,mapped,stage,128)
    if ui:
        page=json.loads((ui/"workbench/case-data.js").read_text("utf8").split("=",1)[1].strip().rstrip(";"))
        fields={v["key"]:v for v in page["fields"]}
        np.testing.assert_array_equal(np.array(page["faces"])+1,body["F"])
        for k,ref in (("pressure",w["absolutePa"]/1000),("pressureGauge",w["gaugePa"]),("heatFlux",q/1000),
                      ("wallInput",w["wallInputK"]),("nearWallCellTemperature",w["firstK"])):
            close(fields[k]["values"],ref)
        for col in ("shearPa","shearMagnitudePa","speedMps","speedMagnitudeMps","status"):
            key="supportAngular"+col[0].upper()+col[1:];ref=d[col].copy()
            if col=="status":ref[ref==0]=np.nan
            np.testing.assert_allclose(np.array(fields[key]["values"],float),ref,rtol=1e-12,atol=1e-12,equal_nan=True)
        assert all(all(v is None for v in fields[k]["values"]) for k in ("wallSolved","edgeTemperature","recoveryTemperature"))
        with (ui/"section.csv").open(encoding="utf-8-sig",newline="") as handle:cut=list(csv.DictReader(handle))
        assert cut
        for row in cut:
            assert row["field_key"]=="heatFlux" and row["unit"]=="kW/m²"
            close(float(row["section_x_m"]),.00137);close(float(row["value"]),q[int(row["face_id"])-1]/1000)
        out.update(offlinePageVerified=True,sectionRows=len(cut))
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--data",type=Path,default=DATA)
    parser.add_argument("--run",type=Path);parser.add_argument("--ui",type=Path);parser.add_argument("--out",type=Path)
    args=parser.parse_args();result=verify(args.data,args.run,args.ui)
    if args.out:
        with args.out.open("x",encoding="utf8",newline="\n") as f:json.dump(result,f,indent=2);f.write("\n")
    print(result)
