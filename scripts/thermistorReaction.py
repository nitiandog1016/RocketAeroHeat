"""从固定珠体的强等温终场读取壁面反力。两次读数使用各自完整的输入态。"""
import argparse
import os
from pathlib import Path
import re
import uuid

from thermistorBead import ENGINE_SHA, ROOT, read, sha, solve, write_json


def terminal(folder):
    record=read(folder/"engine_process.json")
    assert record["status"]=="terminal" and record["exitCode"]==0 and all(record["inputsStable"].values())
    assert record["inputSHA256"][record["command"][0]]==ENGINE_SHA
    assert not (record["timedOut"] or record["cancelled"] or record.get("supervisorError"))
    assert not read(folder/"resources.json")["resourceExceeded"]
    for name in ("case.cfg","mesh.su2","initial_flow.csv"):
        assert sha(folder/name)==record["inputSHA256"][str(folder/name)]
    for name in ("restart.csv","history.csv"):
        assert sha(folder/name)==record["outputSHA256"][name]
    return record


def observe(parent, out, *, release, wall_value=1.0, maximum_seconds=60,
            memory_bytes=2*1024**3, output_bytes=180*1024**2):
    terminal(parent)
    cfg=(parent/"case.cfg").read_text("ascii")
    strong_wall=f"MARKER_WALL_SPECIES= (wall, VALUE, {float(wall_value)})"
    assert strong_wall in cfg
    cfg,count=re.subn(r"^ITER= \d+$","ITER= "+("1" if release else "20"),cfg,flags=re.M)
    assert count==1
    if release:
        # 同一输入态取消该行的强约束，令外加壁通量为零。未约束 RHS 即原壁反力的相反数。
        cfg=cfg.replace(strong_wall,"MARKER_WALL_SPECIES= (wall, FLUX, 0.0)")
        cfg=cfg.replace("VOLUME_OUTPUT= (COORDINATES, SOLUTION, PRIMITIVE)",
                        "VOLUME_OUTPUT= (COORDINATES, SOLUTION, RESIDUAL)")
        cfg+="WRT_RESTART_COMPACT= NO\n"
    out.mkdir(parents=True,exist_ok=False)
    os.link(parent/"mesh.su2",out/"mesh.su2")
    os.link(parent/"restart.csv",out/"initial_flow.csv")
    (out/"case.cfg").write_text(cfg,encoding="ascii")
    sources=[out/"initial_flow.csv",Path(__file__).resolve(),parent/"engine_process.json",ROOT/"scripts/thermistorBead.py"]
    write_json(out/"engine_request.json",dict(workDirectory=str(out),configFile="case.cfg",meshFile="mesh.su2",
              maximumSeconds=maximum_seconds,threads=1,jobToken=str(uuid.uuid4()),sourceFiles=[p.relative_to(ROOT).as_posix() for p in sources]))
    solve(out,memory_bytes=memory_bytes,output_bytes=output_bytes)
    terminal(out)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent",type=Path)
    parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    out=args.new_output.resolve()
    out.relative_to(ROOT/"outputs")
    observe(args.parent.resolve(),out/"reaction0",release=True)
    observe(args.parent.resolve(),out/"plus20",release=False)
    observe(out/"plus20",out/"reaction20",release=True)
