"""比较单射线旋转与保存全部射线的同速基流初猜，读取同一内迭代的节点方程右端。"""
import argparse
import os
from pathlib import Path
import sys

import numpy as np
from matched_step import (CONTRACT, ENGINE, ENGINE_SHA, FULL, FULL_CONFIG_SHA, HERE, RMS, ROOT,
                          base_state, frame, observe, read, save, sha, terminal)

LIMITS = HERE/"TRANSFER_CHECK.md"


def reflect(base_dir, mapped, out):
    """绕 X 轴镜像延拓标量及矢量，端点共用一份值。此操作仅构造初始时刻。"""
    terminal(mapped)
    base,_,report=base_state(base_dir)
    assert report["acceptedForInitialization"]
    with np.load(FULL/"nodes.npz") as grid:
        xyz,ids=grid["xyz"],grid["ids"]
    with np.load(base_dir/"nodes.npz") as sector:
        old_ids=sector["ids"]
    nodes,first=np.unique(ids[:-1].ravel(),return_index=True)
    np.testing.assert_array_equal(nodes,np.arange(len(xyz)))
    ray,station=np.divmod(first,ids.shape[1])
    folded=ray%8
    k=np.minimum(folded,8-folded)
    handedness=np.where(folded>4,-1,1)
    native=base[old_ids[k,station]]
    alpha=k*np.pi/24
    phi=ray*np.pi/24
    ur=native[:,6]*np.cos(alpha)+native[:,7]*np.sin(alpha)
    ut=handedness*(-native[:,6]*np.sin(alpha)+native[:,7]*np.cos(alpha))
    proj=frame(mapped/"seed_00000.csv",xyz)
    impulse=proj[:,5]-base[old_ids[0,station],5]
    ux=native[:,5]+impulse
    full=np.column_stack((nodes,xyz,native[:,4],ux,ur*np.cos(phi)-ut*np.sin(phi),
                          ur*np.sin(phi)+ut*np.cos(phi),native[:,8:10]))
    # 复原原扇区的五条射线。每转 60° 标量、轴向速度和旋回局部坐标的横向分量重合。
    calm=full[:,4:].copy()
    calm[:,1]-=impulse
    np.testing.assert_allclose(calm[ids[:5]],base[old_ids,4:],rtol=0,atol=1e-13)
    ux0=calm[:,1][ids[:-1]]
    scalar=full[:,9][ids[:-1]]
    np.testing.assert_allclose(ux0,np.roll(ux0,8,axis=0),rtol=0,atol=1e-14)
    np.testing.assert_array_equal(scalar,np.roll(scalar,8,axis=0))
    vp=full[:,6:8][ids[:-1]]
    rotations=np.arange(48)[:,None]*np.pi/24
    for values in (vp[:,:,0]*np.cos(rotations)+vp[:,:,1]*np.sin(rotations),
                   -vp[:,:,0]*np.sin(rotations)+vp[:,:,1]*np.cos(rotations)):
        np.testing.assert_allclose(values,np.roll(values,8,axis=0),rtol=0,atol=1e-13)
    w=ids[:-1,:65]
    assert np.all(full[w,5:8]==0) and np.all(full[w,9]==1)
    mode=float((2*np.abs(np.fft.rfft(ux[ids[:-1]],axis=0)[1])/48).max())
    geom=read(mapped/"geometry.json")
    assert abs(mode-geom["peakAxialPerturbationOverU"])<1e-12
    with (mapped/"seed_00000.csv").open(encoding="ascii") as f:
        header=f.readline().strip()
    np.savetxt(out/"seed_00000.csv",full,fmt=["%d"]+["%.17g"]*9,delimiter=",",header=header,comments="")
    geom["initialGuessTransfer"]="mirrored five native rays and reflected tangential velocity"
    save(out/"geometry.json",geom)
    save(out/"mirror_check.json",dict(fullNodes=len(nodes),originalSectorNodes=len(base),
         matchedBaseSHA256=report["frameSHA256"],axialM1=mode,seedSHA256=sha(out/"seed_00000.csv"),
         peakAxialChangeFromProjection=float(np.max(abs(ux-proj[:,5])))))


def observed(stage,source,out,base):
    """仅在新目录运行，一次内迭代只作残差观察。完整复核需先满足已登记的收益条件。"""
    out.mkdir(parents=True,exist_ok=False)
    files=[CONTRACT,LIMITS,Path(__file__).resolve(),HERE/"matched_step.py"]
    if stage=="prepare":
        files.extend([base/"restart.csv",base/"nodes.npz",source/"seed_00000.csv",source/"geometry.json"])
        command=[sys.executable,str(Path(__file__).resolve()),"reflect",str(source),str(out),"--base",str(base)]
        seconds,space=90,80
    else:
        terminal(source)
        assert sha(FULL/"case.cfg")==FULL_CONFIG_SHA and sha(ENGINE)==ENGINE_SHA
        for path in (FULL/"mesh.su2",FULL/"nodes.npz",FULL/"body.json",source/"seed_00000.csv",source/"geometry.json"):
            os.link(path,out/path.name)
            files.append(out/path.name)
        cfg=(FULL/"case.cfg").read_text("ascii")
        if stage=="one":
            assert cfg.count("INNER_ITER= 180\n")==1 and cfg.count("VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)")==1
            cfg=cfg.replace("INNER_ITER= 180\n","INNER_ITER= 1\n").replace(
                "VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)","VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE, RESIDUAL)")
        elif stage=="full":
            decision=read(base)
            assert decision["runFull"] and sha(out/"seed_00000.csv")==decision["mirror"]["seedSHA256"]
            files.append(base)
        (out/"case.cfg").write_text(cfg,encoding="ascii")
        if stage=="full":
            assert sha(out/"case.cfg")==FULL_CONFIG_SHA
        files.extend([out/"case.cfg",ENGINE,FULL/"case.cfg",source/"debugger_process.json"])
        command=[str(ENGINE),"case.cfg","--threads","1"]
        seconds,space=(60,180) if stage=="one" else (360,180)
    rec=observe(command,out,{str(p):sha(p) for p in files},seconds,space*1024**2)
    return rec["debuggerExitCode"]==0 and rec["stopReason"] is None


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=("prepare","reflect","one","full"))
    parser.add_argument("source",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--base",type=Path)
    args=parser.parse_args()
    src,out=args.source.resolve(),args.output.resolve()
    for folder in (src,out):
        folder.relative_to(ROOT/"outputs/v288_validation")
    base=args.base.resolve() if args.base else None
    if args.stage=="reflect":
        reflect(base,src,out)
    else:
        raise SystemExit(0 if observed(args.stage,src,out,base) else 1)
