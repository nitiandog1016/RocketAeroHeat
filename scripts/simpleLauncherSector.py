"""固定低速风洞试件的对称扇区及强等温空气温差。"""
import argparse
from pathlib import Path
import uuid

import numpy as np

from cm011Tunnel import laminar_config, write_json
from cm011TunnelSector import mesh, write_volume
from simpleLauncherTunnel import CASE, CP, GRIDS, K, MU, P_REF, RHO, ROOT, SPEED, T_AIR, T_WALL, geometry
from thermistorBead import solve


def prepare_and_solve(folder, grid):
    folder = folder.resolve()
    folder.relative_to(ROOT / "outputs")
    definition, _, diameter, wall = geometry()
    level, layers, extent = GRIDS[grid]
    refined = np.concatenate([a+(b-a)*np.arange(level)[:, None]/level
                              for a, b in zip(wall[:-1], wall[1:])] + [wall[-1:]])
    # 单个 7.5° 面片正好对应原 STL 周向相邻两列。径向及母线细化保持原壁。
    rays = np.deg2rad([0., 7.5])
    _, xyz, _, _, ids, cells, markers, stats = mesh(rays, level, layers, refined, extent, diameter)
    if abs(stats["fullChordWallAreaM2"]/definition["areaM2"]-1) > 1e-6:
        raise ValueError("扇区尚未覆盖全部原壁")
    folder.mkdir(parents=True, exist_ok=False)
    write_volume(folder, xyz, cells, markers)
    wet = np.unique(ids[:, :len(refined)])
    theta = np.zeros(len(xyz))
    theta[wet] = 1.
    velocity = np.zeros_like(xyz)
    velocity[:, 0] = -1.
    velocity[wet] = 0.
    initial = np.column_stack((np.arange(len(xyz)), xyz, np.zeros(len(xyz)), velocity, 1.-298.15+theta, theta))
    np.savetxt(folder/"initial_flow.csv", initial, delimiter=",", fmt=["%d"]+9*["%.17g"], comments="",
               header='"PointID","x","y","z","Pressure","Velocity_x","Velocity_y","Velocity_z","Enthalpy","Species_0"')
    re, pe = RHO*SPEED*diameter/MU, RHO*CP*SPEED*diameter/K
    cfg = laminar_config(re, pe, 16000).replace("AXISYMMETRIC= YES","AXISYMMETRIC= NO")
    cfg = cfg.replace("MARKER_SYM= (axis)", "MARKER_SYM= (sym0, sym1)")
    cfg = cfg.replace("RMS_VELOCITY-Y, RMS_ENTHALPY", "RMS_VELOCITY-Y, RMS_VELOCITY-Z, RMS_ENTHALPY, RMS_SPECIES_0")
    cfg = cfg.replace("OUTPUT_WRT_FREQ= (1500, 500, 18000)", "OUTPUT_WRT_FREQ= (500, 500, 16000)")
    cfg += f"""RESTART_SOL= YES
READ_BINARY_RESTART= NO
SOLUTION_FILENAME= initial_flow
RESTART_FILENAME= restart
KIND_SCALAR_MODEL= SPECIES_TRANSPORT
SPECIES_INIT= 0.0
DIFFUSIVITY_MODEL= CONSTANT_DIFFUSIVITY
DIFFUSIVITY_CONSTANT= {1/pe:.17g}
MARKER_WALL_SPECIES= (wall, VALUE, 1.0)
CONV_NUM_METHOD_SPECIES= SCALAR_UPWIND
MUSCL_SPECIES= NO
SPECIES_CLIPPING= NO
"""
    (folder/"case.cfg").write_text(cfg, encoding="ascii")
    np.savez_compressed(folder/"mapping.npz", xyz=xyz, ids=ids, wall=refined, coarseWall=wall)
    write_json(folder/"condition.json",dict(grid=grid, layers=layers, extentD=extent, diameterM=diameter,
               speedMps=SPEED, airK=T_AIR, wallK=T_WALL, pRefPa=P_REF, reynoldsD=re, pecletD=pe,
               copies=48, surfaceSHA256=definition["fileSHA256"], mesh=stats))
    sources=[CASE/"geometry.json", CASE/definition["file"], Path(__file__), ROOT/"scripts/simpleLauncherTunnel.py",
             ROOT/"scripts/cm011TunnelSector.py", ROOT/"scripts/cm011Tunnel.py", folder/"condition.json", folder/"mapping.npz"]
    write_json(folder/"engine_request.json",dict(workDirectory=str(folder),configFile="case.cfg",meshFile="mesh.su2",
               threads=1,maximumSeconds=720,jobToken=str(uuid.uuid4()),
               sourceFiles=[p.resolve().relative_to(ROOT).as_posix() for p in sources]))
    solve(folder, memory_bytes=4*1024**3, output_bytes=2*1024**3)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder",type=Path)
    parser.add_argument("grid",choices=GRIDS)
    args=parser.parse_args()
    prepare_and_solve(args.folder,args.grid)
