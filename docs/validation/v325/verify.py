"""独立以原 XYZ 裁三角形，复核球面交叠，并用单壁片脉冲检查发布矩阵的五区守恒。"""
import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2];CASES=ROOT/"cases/short_sphere"
sys.path.insert(0,str(ROOT/"scripts"))
from runEmbeddedSU2 import sha


def cut(points,n):
    if not len(points):return points
    signed=points@n
    if np.min(signed)>=0:return points
    result=[]
    for j,p in enumerate(points):
        a,b=signed[j-1],signed[j]
        if (a<0)!=(b<0):result.append((a*p-b*points[j-1])/(a-b))
        if b>=0:result.append(p)
    return np.array(result)


def area(poly):
    if len(poly)<3:return 0.
    v=poly-poly[0]
    return np.linalg.norm(np.cross(v,np.roll(v,-1,axis=0)).sum(0))/2


def verify():
    report=json.loads((CASES/"supportedGrid65.json").read_bytes())
    assert sha(CASES/"supported65.zip")==report["bodySHA256"] and sha(CASES/"supportedTransfer65.npz")==report["publicTransferSHA256"]
    with ZipFile(CASES/"supported65.zip") as z:
        body=json.loads(z.read("geometry.json"));info=json.loads(z.read("summary.json"))
    D=info["diameterM"]
    points=np.array(body["V"])/D;tri=points[np.array(body["F"])-1]
    with np.load(CASES/"supportedTransfer65.npz") as z:
        fields={k:z[k] for k in z.files}
    i,j,ft,fs=(fields[k] for k in ("target","source","targetFraction","sourceFraction"))
    old,source=fields["targetArea"],fields["nativeArea"]
    p=fields["nativeWallXYZD"][fields["nativeCorners"]]
    sg=fields["nativeRegion"];tg=fields["region"]
    assert len(np.unique(np.c_[i,j],axis=0))==len(i) and np.all(tg[i]==sg[j])
    np.testing.assert_allclose(old/D**2,np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2,atol=2e-15,rtol=0)
    np.testing.assert_allclose(source/D**2,np.linalg.norm(np.cross(p[:,2]-p[:,0],p[:,3]-p[:,1]),axis=1)/2,atol=2e-15,rtol=0)
    normals=np.cross(p[:,2]-p[:,0],p[:,3]-p[:,1])
    old_normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    assert np.min(np.sum(old_normal[i]*normals[j],axis=1))>0
    # 每个原生球壁四边形的边和球心定义四个锥平面；接头的最外边由原实体折线代替。
    # 从目标真实三角形开始裁剪，与生产侧“把原面投到扇区平面”实现独立。
    err=0.;count=0
    for r in np.flatnonzero(tg[i]==1):
        f=p[j[r]];poly=tri[i[r]];row=j[r]%94
        for k in range(4):
            if (row==23 and k==0) or (row==70 and k==2):continue
            n=np.cross(f[k],f[(k+1)%4])
            if f.mean(0)@n<0:n=-n
            poly=cut(poly,n)
        value=area(poly)*D**2/old[i[r]]
        err=max(err,abs(value-ft[r]));count+=1
    assert err<1e-7,err
    np.testing.assert_allclose(np.bincount(i,ft),1,rtol=0,atol=1e-7)
    np.testing.assert_allclose(np.bincount(j,fs),1,rtol=0,atol=1e-7)
    # 给每个分区单独加一个非均匀正负载荷，检查跨区泄漏和原面积积分。
    pulse_errors=[]
    for group in range(1,6):
        load=np.where(sg==group,np.sin(1+np.arange(len(sg))*.17),0.)
        restored=np.bincount(i,fs*load[j]*source[j])/old
        assert np.all(restored[tg!=group]==0)
        pulse_errors.append(float(abs(restored@old-load@source)))
    assert max(pulse_errors)<1e-15
    resources={};base=ROOT/"outputs/v325_validation"
    for kind,stage,prep in (("hybrid","run01","prep02"),("structured","run02","prep03")):
        d=report[kind];path=base/stage/"steady01";proc=json.loads((path/"debugger_process.json").read_bytes())
        assert proc["status"]=="terminal" and proc["debuggerExitCode"]==0 and all(proc["inputsStable"].values())
        assert not(proc["stopReason"] or proc["survivingProcessIDs"] or proc["supervisorError"])
        assert sha(path/"mesh.su2")==d["meshSHA256"] and (base/stage/"failure.json").is_file() and not (base/stage/"endpoints.json").exists()
        with np.load(base/prep/"prepared/mapping.npz") as m:
            np.testing.assert_array_equal(p,m["xyz"][m["corners"]])
            np.testing.assert_array_equal(fields["firstXYZD"],m["xyz"][m["first"]])
        for frame in d["stages"][0]["frames"]:
            file=path/f"restart_{frame['index']:06d}.csv";assert sha(file)==frame["sha256"]
            scalar=np.loadtxt(file,skiprows=1,delimiter=",",usecols=9)
            assert float(scalar.min())==frame["minScalar"] and int((scalar< -1e-12).sum())==frame["below"]>0
        resources[kind]=[proc["elapsedSeconds"],proc["peakFamilyRSSBytes"]]
    ui=base/"ui01";receipt=json.loads((ui/"debugger_process.json").read_bytes());test=json.loads((ui/"ui.json").read_bytes())
    assert receipt["status"]=="terminal" and receipt["debuggerExitCode"]==0 and all(receipt["inputsStable"].values())
    assert not(receipt["stopReason"] or receipt["survivingProcessIDs"] or receipt["supervisorError"])
    rebuilt=Path(test["rebuiltJob"])
    assert sha(rebuilt/"mesh.su2")==report["hybrid"]["meshSHA256"] and test["transferFaces"]==len(old)
    with np.load(rebuilt/"transfer.npz") as z:
        for k in z.files:np.testing.assert_array_equal(z[k],fields[k])
    return dict(originalFaces=len(old),nativePatches=len(p),mappedPairs=len(i),sphereOverlapsChecked=count,
                maxIndependentSphereFractionError=err,pulseErrorsW=pulse_errors,programResources=resources,
                uiSeconds=receipt["elapsedSeconds"],uiMaxRSSBytes=receipt["peakFamilyRSSBytes"],ui=test)


if __name__=="__main__":
    result=verify()
    with Path(sys.argv[1]).open("x",encoding="utf8") as f:json.dump(result,f,indent=2);f.write("\n")
