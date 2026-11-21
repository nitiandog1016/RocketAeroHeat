"""为带线 B 型热敏珠生成 3 m/s 侧吹全周冷态输入；几何复制不复制任何流场。"""
import argparse
from pathlib import Path
import shutil

import numpy as np
from scipy.spatial import cKDTree

from runThermistorLeadCase import configuration
from thermistorBead import D, ENGINE, ENGINE_SHA, ROOT, read, sha, write_json
from thermistorLeadMesh import SECTORS, WEDGE, rotate
from thermistorOuterDomain import load_mesh, write_mesh


def expand(base,out):
    """沿已发布的真实扇区连接关系生成 48 个周向间隔，在轴与接缝处统一编号。"""
    contract=read(ROOT/"cases/thermistor_b/leadsMidbody3.json")
    sources={"mesh.su2":"meshSHA256","nodes.npz":"nodeMapSHA256","body.json":"bodySHA256","case.json":"caseSHA256",
             "case.cfg":"configSHA256","initial_flow.csv":"initialSHA256"}
    for name,key in sources.items():
        assert sha(base/name)==contract[key],name
    assert sha(ENGINE)==ENGINE_SHA
    cells,sector_xyz,markers=load_mesh(base/"mesh.su2")
    with np.load(base/"nodes.npz") as z:
        sector=z["ids"];areas=z["area"]
        np.testing.assert_array_equal(z["xyz"],sector_xyz)
    xy=sector_xyz[sector[0],:2]
    xyz,rays=rotate(xy,np.arange(SECTORS+1)*2*np.pi/SECTORS,closed=True)
    np.testing.assert_allclose(xyz[rays[:WEDGE+1]],sector_xyz[sector],atol=3e-14,rtol=0)
    np.testing.assert_array_equal(np.unique(sector),np.arange(len(sector_xyz)))
    full=[];boundary={"wall":[],"far":[]}
    wedges=SECTORS//WEDGE
    for start in range(0,SECTORS,WEDGE):
        ids=np.full(len(sector_xyz),-1,dtype=int)
        ids[sector]=rays[start:start+WEDGE+1]
        assert np.min(ids)>=0
        # 转换的是节点号。有序面与体单元的绕序沿旋转方向保持。
        local=ids.tolist()
        full.extend(tuple(local[v] for v in cell) for cell in cells)
        for name in boundary:
            boundary[name].extend(tuple(local[v] for v in face) for face in markers[name])
    assert len(full)==wedges*contract["cells"]
    out.mkdir(parents=True,exist_ok=False)
    write_mesh(out/"mesh.su2",full,xyz,boundary)
    body=read(base/"body.json");spec=read(base/"case.json");f=np.asarray(body["F"])-1
    width=spec["meridianIntervals"]+1;surface=rays[:-1,:width]
    ids=np.unique(surface)
    distance,index=cKDTree(xyz[ids]).query(np.asarray(body["V"])/D)
    assert max(distance)<1e-12 and len(np.unique(index))==len(ids)==len(index)
    wall=ids[index][f];next_node=np.full(len(xyz),-1,dtype=int)
    next_node[surface]=rays[:-1,width:2*width]
    first=next_node[wall];assert np.all(first>=0)
    np.savez_compressed(out/"nodes.npz",xyz=xyz,ids=rays,wallMap=wall,firstMap=first,area=areas)
    shutil.copyfile(base/"body.json",out/"body.json")
    spec.update(caseKind="thermistor-B-sidewind-two-leads",wallK=299.15,velocityMps=[0,-3,0],sectorDeg=360,angularCells=SECTORS)
    write_json(out/"case.json",spec)

    # 有量纲温差尺度为 6 K。强标量取环境 0、壁 1；弱焓独立作诊断。
    reference=configuration()
    assert (base/"case.cfg").read_text("ascii")==reference
    cfg=reference.replace("MARKER_SYM= (sym0, sym1)\n","").replace("INC_VELOCITY_INIT= (-1, 0, 0)","INC_VELOCITY_INIT= (0, -1, 0)")
    assert "MARKER_SYM" not in cfg and cfg.count("INC_VELOCITY_INIT= (0, -1, 0)")==1
    (out/"case.cfg").write_text(cfg.replace("OUTPUT_WRT_FREQ= (500, 500)","OUTPUT_WRT_FREQ= (100, 100)"),encoding="ascii")
    state=np.zeros((len(xyz),10));state[:,0]=np.arange(len(xyz));state[:,1:4]=xyz
    state[:,6]=-1;state[ids,6]=0
    state[:,8]=1-298.15;state[ids,8]=2-298.15;state[ids,9]=1
    with (base/"initial_flow.csv").open(encoding="ascii") as stream:
        header=stream.readline().strip()
    np.savetxt(out/"initial_flow.csv",state,delimiter=",",comments="",header=header,fmt=["%d"]+["%.17g"]*9)
    names=("mesh.su2","nodes.npz","body.json","case.json","case.cfg","initial_flow.csv")
    write_json(out/"inputs.json",dict(base=str(base),source={n:contract[k] for n,k in sources.items()},
        files={n:sha(out/n) for n in names},nodes=len(xyz),cells=len(full),axisNodes=int(np.sum(xy[:,1]==0)),
        wedgeCopies=wedges,sectionFaces=len(markers["sym0"]),engineSHA256=ENGINE_SHA))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base",type=Path);parser.add_argument("out",type=Path)
    args=parser.parse_args();dest=args.out.resolve();dest.relative_to(ROOT/"outputs")
    expand(args.base.resolve(),dest)
