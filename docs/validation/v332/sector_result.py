"""保存 Short Test 65 加密网格交叠，核对两个十层网的同面热流和流向。"""
import argparse
import io
import json
from pathlib import Path
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

import numpy as np
from sector_grid import HERE, ROOT, NEW_KEY, read, sha

sys.path[:0]=[str(HERE.parent/v) for v in ("v331","v330","v329")]
from direction import FLAGS, geometry, mapped_fields, native_vectors, report_grid
from compare import compare_surfaces, differences, load_surface
from map_surface import mapped

DATA=ROOT/"cases/short_sphere"


def store(path,value):
    with path.open("x",encoding="ascii",newline="\n") as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write("\n")


def publish_transfer(prepared,out):
    assert sha(prepared/"prepared.json")==NEW_KEY
    spec=read(prepared/"prepared.json")
    for name in ("mapping.npz","transfer.npz"):assert sha(prepared/name)==spec["files"][name]
    with np.load(prepared/"mapping.npz") as m,np.load(prepared/"transfer.npz") as t:
        wall=m["wet"];xyz=m["xyz"]
        with out.open("xb") as handle:
            np.savez_compressed(handle,**dict(t),nativeWallXYZD=xyz[wall],firstXYZD=xyz[m["first"]],
                                nativeCorners=np.searchsorted(wall,m["corners"]),nativeArea=m["nativeArea"],nativeRegion=m["group"])


def process(work,prep,grid,transfer,out):
    # 命令行入口已建立自己的任务目录，或由独立后处理在此创建；结果文件只写入新路径。
    out.mkdir(parents=True,exist_ok=True)
    assert grid["preparedSHA256"]==sha(prep/"prepared.json")==NEW_KEY
    pair=read(DATA/"supportedNearWall65.json");baseline=DATA/"supportedFine65.zip"
    assert sha(baseline)==pair["fineArchiveSHA256"]
    current=mapped(work,out/"surface",prepared=prep,expected_mesh=grid["meshSHA256"],
                   transfer_file=transfer,transfer_sha=grid["publicTransferSHA256"])
    comparison=compare_surfaces(baseline,out/"surface/surface.zip",grid,angular=True)
    base,bw,bs=load_surface(baseline)
    report,wall,stations=load_surface(out/"surface/surface.zip")
    c=read(work/"converged.json");assert c["frame"]==current["frameSHA256"]
    stage=work/c["stage"];vtk=sha(stage/"surface_000199.vtk")
    assert vtk==sha(stage/"surface.vtk")
    geom=geometry(current,transfer);keys=dict(frame=c["frame"],vtkSHA256=vtk)
    vectors,wet,scales,_=native_vectors(prep,stage,geom["m"],current["condition"],keys)
    fields=mapped_fields(geom,vectors,scales)
    d=report_grid("fine",geom,fields,wet,vectors,scales,wall,stations,comparison,sectors=128)
    prior=read(DATA/"supportedDirection65.json");assert base["frameSHA256"]==prior["fine"]["frame"]
    with ZipFile(DATA/"supportedDirection65.zip") as z:
        assert sha(DATA/"supportedDirection65.zip")==prior["archiveSHA256"]
        old=np.genfromtxt(io.BytesIO(z.read("direction.csv")),names=True,delimiter=",")
    np.testing.assert_array_equal(old["faceID"],wall["faceID"]);np.testing.assert_array_equal(old["regionID"],wall["regionID"])
    area=wall["areaM2"];a=old["fineStatus"]==-1;b=fields[:,4]==-1
    # 差值为同一物理量、同一原面。温界按同一给定极角函数重新取样，各面平均允许细小差异。
    comparison.update(schema="short-65-angular-64-128-v1",baseSectors=64,sectors=128,flagNames=FLAGS,
        geometry=grid,baseTransferSHA256=base["transferSHA256"],wallTemperatureDifferenceK=differences(abs(wall["wallInputK"]-bw["wallInputK"]),area),
        baseDirection=prior["fine"],direction=d,sourceVTKSHA256=vtk,
        reverseOverlap=dict(sharedM2=area[a&b].sum(),unionM2=area[a|b].sum(),baseOnlyM2=area[a&~b].sum(),refinedOnlyM2=area[b&~a].sum()),
        heatFluxAtBasePeakWm2=float(-wall["fluxInWm2"][comparison["coarsePeak"]["faceID"]-1]))
    assert comparison["coarsePeak"]["faceID"]==25055
    csv=io.StringIO();np.savetxt(csv,np.column_stack((wall["faceID"],wall["regionID"],fields)),delimiter=",",
        header="faceID,regionID,shearPa,shearMagnitudePa,speedMps,speedMagnitudeMps,status,flags",comments="",fmt=["%d","%d"]+["%.17g"]*4+["%d","%d"])
    direction=out/"direction.zip"
    with ZipFile(direction,"x") as z:z.writestr(ZipInfo("direction.csv"),csv.getvalue(),compress_type=ZIP_DEFLATED,compresslevel=9)
    comparison["directionArchiveSHA256"]=sha(direction)
    comparison["flowRangeMax"]={col:float(np.nanmax(abs(np.r_[old["fine"+col],fields[:,j]]))) for j,col in enumerate(
        ("ShearPa","ShearMagnitudePa","SpeedMps","SpeedMagnitudeMps","Status"))}
    store(out/"comparison.json",comparison)
    print(current["thermal"]["updates"],comparison["fineW"],comparison["finePowerError"],d["transitions"])
    return current


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("work",type=Path);parser.add_argument("prep",type=Path)
    parser.add_argument("grid",type=Path);parser.add_argument("transfer",type=Path);parser.add_argument("out",type=Path)
    a=parser.parse_args();process(a.work,a.prep,read(a.grid),a.transfer,a.out)
