"""以原生 CSV、独立交点枚举和 STL 真值复核五条路径、首层位置及截线。"""
import csv
import hashlib
import io
import json
import sys
from pathlib import Path
from zipfile import ZipFile

import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/"v314"))
from audit_angular import independent_overlap


def read(path):return json.loads(path.read_text("utf-8-sig"))


def raw(stage):
    frame=stage/"restart.csv";reaction=stage/"reaction0/restart.csv"
    spec=stage/"case.json"
    if not spec.is_file():spec=stage.parent/"prepared/case.json"
    with reaction.open() as stream:
        col=next(csv.reader(stream)).index("Residual_Species_0")
        rhs=np.loadtxt(stream,delimiter=",",usecols=col)
    return read(spec),np.loadtxt(frame,delimiter=",",skiprows=1,usecols=(4,9)),rhs


def archive(stem,name="wall.csv"):
    with ZipFile(ROOT/"cases/thermistor_b"/(stem+".zip")) as z:
        return np.genfromtxt(io.BytesIO(z.read(name)),delimiter=",",names=True)


def validate(ui):
    proof=read(ui/"ui.json");log=read(ui/"debugger_process.json")
    assert log["status"]=="terminal" and log["debuggerExitCode"]==0 and all(log["inputsStable"].values())
    assert not(log["survivingProcessIDs"] or log["stopReason"] or log["supervisorError"])
    fixture=ui.parent/"fixture";prep=Path(proof["prepared"]);body=read(prep/"body.json")
    truth=np.loadtxt(fixture/"face_truth.csv",delimiter=",",dtype=int);order=truth[:,0]-1
    split=np.where(truth[:,1]>0)[0];np.testing.assert_array_equal(split+1,proof["children"])
    stl=(fixture/"probe-partition-mm.stl").read_bytes()
    assert hashlib.sha256(stl).hexdigest()==proof["stlSHA256"]
    dt=np.dtype([("normal","<f4",3),("points","<f4",(3,3)),("attribute","<u2")])
    tri=np.frombuffer(stl,dtype=dt,offset=84)["points"].astype(float)*.001
    coarse=np.asarray(body["V"])[np.array(body["F"])-1]
    def areas(t):return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
    area=areas(tri);parent_area=areas(coarse)
    # 用距离与正交投影独立重建文件角点在父平面的位置。顶点与细壁的交点由另一算法枚举。
    target=tri[split].copy()
    for k,i in enumerate(split):
        p=coarse[order[i]];normal=np.cross(p[1]-p[0],p[2]-p[0]);normal/=np.linalg.norm(normal)
        for j,x in enumerate(tri[i]):
            dist=np.linalg.norm(p-x,axis=1)
            target[k,j]=p[dist.argmin()] if dist.min()<2e-10 else x-(x-p[0])@normal*normal
    expected=areas(target)
    with np.load(prep/"nodes.npz") as z:
        wall,near,xyz=z["wallMap"],z["firstMap"],z["xyz"]*read(prep/"case.json")["diameterM"]
    # 对两个物理邻区独立扩大候选集，不能只依赖生产映射的非零交叠清单。
    native=xyz[wall];fine_area=areas(native)
    by_parent={}
    for p in set(order[split]):
        t=coarse[p];normal=np.cross(t[1]-t[0],t[2]-t[0]);normal/=np.linalg.norm(normal)
        intersects=np.all(native.min(axis=1)<=t.max(axis=0)+1e-14,axis=1)&np.all(native.max(axis=1)>=t.min(axis=0)-1e-14,axis=1)
        coplanar=np.max(abs((native-t[0])@normal),axis=1)<1e-14
        by_parent[p]=np.flatnonzero(intersects & coplanar)
    pairs=np.array([(i,j) for i,parent in enumerate(order[split]) for j in by_parent[parent]])
    dest,src=pairs.T;w=independent_overlap(target,native,dest,src)*fine_area[src,None]
    mapped=np.bincount(dest,w.sum(1));np.testing.assert_allclose(mapped,expected,rtol=1e-10,atol=1e-24)
    def integral(columns):return np.bincount(dest,np.einsum("ij,ij->i",w,columns[src]))
    xyz1=np.column_stack([integral(xyz[near,d])/mapped for d in range(3)])
    wall_center=np.mean(target,axis=1)
    normal=np.cross(target[:,1]-target[:,0],target[:,2]-target[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    distance=np.einsum("ij,ij->i",xyz1-wall_center,normal);assert np.all(distance>0)
    source_root=ROOT/"outputs"
    nativeJob=Path(read(source_root/"v319_validation/native01/ui.json")["native"]["directory"])/"run"
    stages=dict(wire=nativeJob/read(nativeJob/"completed.json")["stage"],corrected=source_root/"v318_validation/check01/steady01",
                uniform=source_root/"v316_validation/run01/steady04",six=source_root/"v315_validation/six01",twelve=source_root/"v315_validation/native01")
    lump=np.bincount(wall.ravel(),np.repeat(fine_area/3,3),minlength=len(xyz));wet=np.unique(wall)
    data={}
    for key,path in stages.items():
        c,state,rhs=raw(path);delta=c.get("scalarReferenceDeltaK",c["wallK"]-c["airK"])
        flux=np.zeros(len(state))
        flux[wet]=rhs[wet]*c["rho"]*c["speedMps"]*c["cp"]*delta*c["diameterM"]**2/lump[wet]
        data[key]=dict(flux=integral(flux[wall])/area[split],gauge=integral(state[wall,0]*c["rho"]*9)/mapped,
                       wall=integral(state[wall,1]*delta)/mapped,first=integral(state[near,1]*delta)/mapped)
    factor=read(ROOT/"cases/thermistor_b/resistance12.json")["factor"]
    cases=[("longWire_archive","wire",1,"longWire12","wall.csv"),("resistance_solve","wire",factor,"longWire12","wall.csv"),
           ("resistance_archive","corrected",1,"resistance12","wall.csv"),("uniform_solve","six",2,"windPoint12","six_wall.csv"),
           ("uniform_archive","twelve",1,"windPoint12","wall.csv")]
    result={}
    for action,key,scale,stem,member in cases:
        with (ui/action/"wall.csv").open(encoding="utf-8-sig",newline="") as stream:rows=list(csv.DictReader(stream))
        assert len(rows)==10276
        def column(name):return np.asarray([float(r[name]) for r in rows])
        np.testing.assert_array_equal(column("face_id"),np.arange(1,len(rows)+1))
        np.testing.assert_array_equal(column("reference_face_id"),order+1)
        np.testing.assert_allclose(column("area_m2"),area,rtol=1e-12,atol=0)
        assert all(r["joint_verified"] in ("0","false") for r in rows)
        fields=data[key];q=column("heat_flux")*1000;pressure=column("pressure")*1000
        np.testing.assert_allclose(q[split],scale*fields["flux"],rtol=1e-11,atol=1e-8)
        np.testing.assert_allclose(pressure[split],101325+fields["gauge"],rtol=0,atol=1e-8)
        np.testing.assert_allclose(column("wall_input_K")[split],293.15+scale*fields["wall"],rtol=0,atol=1e-9)
        np.testing.assert_allclose(column("first_fluid_K")[split],293.15+scale*fields["first"],rtol=0,atol=1e-9)
        actualXYZ=np.column_stack([column("subface_first_"+k+"_m")[split] for k in "xyz"])
        np.testing.assert_allclose(actualXYZ,xyz1,rtol=0,atol=1e-15)
        if "wire" in key or key=="corrected":
            control=factor if action.startswith("resistance") else 1.
            np.testing.assert_allclose(column("change_from_uniform_kW_m2")[split]*1000,q[split]-control*data["uniform"]["flux"],rtol=0,atol=1e-8)
        if action.startswith("resistance"):
            np.testing.assert_allclose(column("thermometry_change_kW_m2")[split]*1000,q[split]-data["wire"]["flux"],rtol=0,atol=1e-8)
        ref=archive(stem,member);nx=np.bincount(order,q*area,minlength=len(ref))
        whole=np.flatnonzero(truth[:,1]==0);src_id=order[whole]
        np.testing.assert_allclose(pressure[whole],ref["pressurePa"][src_id],rtol=0,atol=1e-8)
        np.testing.assert_allclose(column("first_fluid_K")[whole],293.15+scale*(ref["nearWallCellTemperatureK"][src_id]-293.15),rtol=0,atol=1e-9)
        wall_input=293.15+scale*(ref["wallInputK"][src_id]-293.15) if "wallInputK" in ref.dtype.names else 305.15
        np.testing.assert_allclose(column("wall_input_K")[whole],wall_input,rtol=0,atol=1e-9)
        groups=np.asarray(body["regionID"])[order]
        assert len(proof[action]["regions"])==7
        for j,region in enumerate(proof[action]["regions"],1):
            use=groups==j
            assert region["faces"]==np.sum(use) and abs(region["areaM2"]-area[use].sum())<1e-18
            assert abs(region["outwardW"]+q[use]@area[use])<1e-12
        diff=nx-scale*ref["areaM2"]*ref["heatFluxWm2"]
        assert abs(diff).max()<1e-12 and abs(nx.sum()-proof[action]["powerW"])<1e-12
        assert proof[action]["glPixels"]>0 and proof[action]["glError"]==0
        detail=[]
        for i in split:
            assert rows[i]["is_subface"] in ("1","true") and "交叠" in rows[i]["heat_flux_method"]
            detail.append(dict(id=int(i+1),parent=int(order[i]+1),areaM2=area[i],gaugePa=pressure[i]-101325,
                               fluxWm2=q[i],inputK=float(rows[i]["wall_input_K"]),firstK=float(rows[i]["first_fluid_K"])))
        result[action]=dict(faces=len(rows),wholeOutW=-nx.sum(),maxParentPowerErrorW=float(abs(diff).max()),children=detail,
                            wallCSV=hashlib.sha256((ui/action/"wall.csv").read_bytes()).hexdigest(),
                            shapeError=proof[action]["shapeError"])
        if action=="longWire_archive":
            for parent,x in [(183,.000574),(184,.000582)]:
                with (ui/action/("section%d.csv"%parent)).open(encoding="utf-8-sig",newline="") as f:cut=list(csv.DictReader(f))
                ids=np.array([int(r["face_id"])-1 for r in cut]);inter=np.flatnonzero((tri[:,:,0].min(axis=1)<x)&(tri[:,:,0].max(axis=1)>x))
                np.testing.assert_array_equal(np.sort(ids),inter)
                np.testing.assert_allclose([float(r["value"]) for r in cut],q[ids]/1000,rtol=0,atol=1e-12)
                assert (order[ids]==parent-1).sum()>=2
                for row,i in zip(cut,ids):
                    points=[]
                    for a,b in zip(tri[i],np.roll(tri[i],-1,axis=0)):
                        if (a[0]-x)*(b[0]-x)<0:points.append((a+(b-a)*(x-a[0])/(b[0]-a[0]))[1:])
                    ends=np.array([[float(row[tag+"_"+k+"_m"]) for k in "yz"] for tag in "ab"])
                    if np.linalg.norm(ends[0]-points[1])<np.linalg.norm(ends[0]-points[0]):points.reverse()
                    np.testing.assert_allclose(ends,points,rtol=0,atol=1e-15)
                result[action]["section%d"%parent]=dict(xM=x,segments=len(cut),childIDs=(ids[order[ids]==parent-1]+1).tolist())
    return dict(areaM2=float(area.sum()),repartitionedAreaM2=float(area[split].sum()),stl=proof["stlSHA256"],
                firstNormalDistanceM=distance.tolist(),uiSeconds=log["elapsedSeconds"],cases=result)


if __name__=="__main__":
    ui=Path(sys.argv[1]);out=Path(sys.argv[2]);data=validate(ui)
    with out.open("x",encoding="utf8") as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write("\n")
