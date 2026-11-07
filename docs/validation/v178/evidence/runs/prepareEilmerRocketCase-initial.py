"""生成固定 Eilmer v5 完整域输入；仅准备，不安装、编译或运行求解器。

模块一：检查原 SU2 网格身份、边界和单元顺序，显式逆转项目楔体约定。
模块二：生成原生热完美空气/热壁/Wilcox k-omega 输入，不迁移 SST 初场。
模块三：独立目录只保存小输入与清单，原网格不复制，原生可运行性保持未知。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PureWindowsPath
import time

import numpy as np

COMMIT = "f53f4609a0331d48efee69a4e4f3c3598378cc03"
TAG = "v5.0.0"
WIDTH = {10: 4, 12: 8, 13: 6, 14: 5}
WEDGE_PERMUTATION = [1, 0, 2, 4, 3, 5]
FACES = {
    10: [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
    12: [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
         [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]],
    13: [[0, 2, 1], [3, 4, 5], [0, 1, 4, 3], [1, 2, 5, 4], [2, 0, 3, 5]],
    14: [[0, 3, 2, 1], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value, name, positive=True):
    require(not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value),
            f"{name} 必须是有限实数")
    require(not positive or value > 0, f"{name} 必须为正")
    return float(value)


def validate_condition(condition):
    require(isinstance(condition, dict), "condition 必须是对象")
    c = dict(condition)
    mode = c.setdefault("flowMode", "explicitVelocity")
    common = {"flowMode", "pressurePa", "temperatureK", "wallTemperatureK", "maxNewtonSteps",
              "stopRelativeResidual", "boundaryMarkers"}
    modes = {
        "explicitVelocity": {"velocityMps", "tkeM2s2", "omegaPerS", "targetMach"},
        "nativeMach": {"Mach", "alphaDeg", "betaDeg", "turbulenceIntensity", "eddyViscosityRatio"},
    }
    require(mode in modes, "flowMode 仅支持 explicitVelocity 或 nativeMach")
    require(not set(c) - common - modes[mode], "工况含未知参数")
    for key in ["pressurePa", "temperatureK", "wallTemperatureK", "maxNewtonSteps", "stopRelativeResidual"]:
        require(key in c, f"缺少 {key}")
        c[key] = number(c[key], key)
    require(c["maxNewtonSteps"].is_integer(), "maxNewtonSteps 必须是正整数")
    c["maxNewtonSteps"] = int(c["maxNewtonSteps"])
    require(c["stopRelativeResidual"] < 1, "相对残差目标必须小于1；不等价于物理收敛")
    if mode == "nativeMach":
        for key in modes[mode]:
            require(key in c, f"缺少 {key}")
            c[key] = number(c[key], key, key not in {"alphaDeg", "betaDeg"})
        require(c["Mach"] > 1, "当前准备合同仅支持超声速外流")
        require(abs(c["alphaDeg"]) < 90 and abs(c["betaDeg"]) < 90, "攻角/侧滑须在(-90,90)度")
        require(c["turbulenceIntensity"] < 1, "湍流强度须在(0,1)")
    else:
        require(isinstance(c.get("velocityMps"), (list, tuple)) and len(c["velocityMps"]) == 3,
                "velocityMps 必须为3个实际速度分量")
        c["velocityMps"] = [number(v, "velocityMps", False) for v in c["velocityMps"]]
        require(math.hypot(*c["velocityMps"]) > 0, "速度不能全为零")
        for key in ["tkeM2s2", "omegaPerS"]:
            require(key in c, f"缺少 {key}")
            c[key] = number(c[key], key)
        if "targetMach" in c:
            c["targetMach"] = number(c["targetMach"], "targetMach")
    markers = c.setdefault("boundaryMarkers", {"wall": "wall", "farfield": "farfield"})
    require(isinstance(markers, dict) and set(markers) == {"wall", "farfield"}, "边界必须显式映射wall/farfield")
    require(all(isinstance(v, str) and v.strip() == v and v and not any(x in v for x in "\r\n\x00")
                for v in markers.values()) and len(set(markers.values())) == 2, "边界名称无效或重复")
    return c


def read_mesh(path):
    """读取原文件，不改节点/单元/边界序；只接受本项目输出的三维SU2子集。"""
    with Path(path).open("r", encoding="utf-8-sig") as stream:
        def line():
            for value in stream:
                value = value.split("%", 1)[0].strip()
                if value:
                    return value
            raise ValueError("SU2文件意外结束")

        def header(name):
            text = line().split("=", 1)
            require(len(text) == 2 and text[0].strip() == name, f"预期SU2头 {name}")
            return text[1].strip()

        require(header("NDIME") == "3", "仅接受三维完整域")
        count = int(header("NELEM"));require(count > 0, "体单元为空")
        groups = {key: [] for key in WIDTH}
        for i in range(count):
            tokens = [int(v) for v in line().split()];kind = tokens[0]
            require(kind in WIDTH, "不支持的体单元类型")
            width = WIDTH[kind]
            require(len(tokens) in (width + 1, width + 2), "单元节点数错误")
            require(len(tokens) == width + 1 or tokens[-1] == i, "单元编号未保持原序")
            ids = tokens[1:width + 1]
            require(len(set(ids)) == width, "单元含重复节点")
            groups[kind].append(ids)
        n = int(header("NPOIN"));require(n > 0, "节点为空")
        points = np.empty((n, 3), dtype=np.float64)
        for i in range(n):
            tokens = line().split();require(len(tokens) in (3, 4), "节点格式错误")
            require(len(tokens) == 3 or int(tokens[-1]) == i, "节点编号未保持原序")
            points[i] = [float(v) for v in tokens[:3]]
        require(np.isfinite(points).all(), "节点含非有限值")
        groups = {key: np.asarray(value, dtype=np.int64).reshape(-1, WIDTH[key]) for key, value in groups.items()}
        for cells in groups.values():
            require(not cells.size or (cells.min() >= 0 and cells.max() < n), "体单元节点越界")
        marker_count = int(header("NMARK"));markers = {}
        for _ in range(marker_count):
            tag = header("MARKER_TAG");require(tag not in markers, "边界名称重复")
            m = int(header("MARKER_ELEMS"));require(m > 0, "空边界标记")
            digest = hashlib.sha256();area = 0.0
            for _ in range(m):
                tokens = [int(v) for v in line().split()];kind = tokens[0]
                require(kind in (5, 9), "边界须为三角或四边形")
                ids = tokens[1:];require(len(ids) == (3 if kind == 5 else 4), "边界节点数错误")
                require(len(set(ids)) == len(ids) and min(ids) >= 0 and max(ids) < n, "边界节点无效")
                digest.update(np.asarray(tokens, dtype="<i8").tobytes())
                p = points[ids];center = p.mean(axis=0)
                piece = sum(np.linalg.norm(np.cross(p[(j + 1) % len(p)] - p[j], center - p[j])) / 2
                            for j in range(len(p)))
                require(piece > 0 and math.isfinite(piece), "边界面积非正")
                area += piece
            markers[tag] = {"faces": m, "areaM2": float(area), "orderedConnectivitySHA256": digest.hexdigest()}
        for value in stream:
            require(not value.split("%", 1)[0].strip(), "原网格尾部含未支持内容")
    return points, groups, markers


def audit_order(points, groups):
    """按Eilmer原生转换得到VTK节点序；中位面离散体积不是全参考Jacobian证明。"""
    total = 0.0;minimum = math.inf
    for kind, original in groups.items():
        for start in range(0, len(original), 40000):
            cells = original[start:start + 40000]
            if kind == 13:
                cells = cells[:, WEDGE_PERMUTATION]
            p = points[cells];center = p.mean(axis=1);vol = np.zeros(len(cells))
            for face in FACES[kind]:
                fc = p[:, face].mean(axis=1) - center
                for j in range(len(face)):
                    x = p[:, face[j]] - center;y = p[:, face[(j + 1) % len(face)]] - center
                    sub = np.einsum("ij,ij->i", x, np.cross(y - x, fc - x)) / 6
                    require(np.isfinite(sub).all() and np.all(sub > 0),
                            f"类型{kind}存在非正中位子体积；拒绝错误节点序/反向/非凸单元")
                    minimum = min(minimum, float(sub.min()));vol += sub
            total += float(np.sum(vol))
    require(total > 0 and math.isfinite(total), "总体积无效")
    return {"discreteVolumeM3": total, "minimumMedianSubVolumeM3": minimum,
            "scope": "每个转换后单元的定向中位面离散体积；未证明全参考Jacobian、全域拓扑、近壁质量或原生读取"}


def native_path(path):
    """Windows本地盘仅生成标准WSL挂载引用；不声称挂载已经存在。"""
    path = Path(path).resolve()
    if os.name == "nt":
        p = PureWindowsPath(str(path))
        require(len(p.drive) == 2 and p.drive[1] == ":", "当前合同不支持UNC网络网格")
        return "/mnt/" + p.drive[0].lower() + "/" + "/".join(p.parts[1:])
    return path.as_posix()


def lua_string(value):
    return json.dumps(value, ensure_ascii=False)


def make_job(mesh_native, c):
    """Lua名称逐项依据固定v5的bc/flowstate/nkconfig及官方SU2案例，不执行原生。"""
    f = lambda value: format(value, ".17g")
    if c["flowMode"] == "nativeMach":
        flow = f'''-- 按同一原生热完美气体求声速/输运，不能使用定gamma替代。
local speed = {f(c['Mach'])} * gs.a
local alpha, beta = math.rad({f(c['alphaDeg'])}), math.rad({f(c['betaDeg'])})
local ux = -speed*math.cos(alpha)*math.cos(beta)
local uy = -speed*math.sin(beta)
local uz = speed*math.sin(alpha)*math.cos(beta)
local tke = 1.5*({f(c['turbulenceIntensity'])}*speed)^2
local omega = gs.rho*tke/(gs.mu*{f(c['eddyViscosityRatio'])})'''
    else:
        flow = f'''-- 显式速度模式；未核对它与任何目标Mach相等。
local ux, uy, uz = {', '.join(f(v) for v in c['velocityMps'])}
local tke, omega = {f(c['tkeM2s2'])}, {f(c['omegaPerS'])}'''
    return f'''-- RocketAeroHeat：仅准备的Eilmer {TAG}输入，原生尚未验证。
-- 固定源码 {COMMIT}；Wilcox k-omega，非SST、非旧重启动。
config.dimensions = 3
config.axisymmetric = false
config.solver_mode = "steady"
config.field_format = "rawbinary"
config.grid_format = "rawbinary"
config.viscous = true
config.reacting = false
config.turbulence_model = "k_omega"
config.turbulence_prandtl_number = 0.89
config.adjust_invalid_cell_data = false
config.max_invalid_cells = 0
local nsp, nmodes, gm = setGasModel("rah-air.gas")
local gs = GasState:new{{gm}}
gs.p, gs.T, gs.massf = {f(c['pressurePa'])}, {f(c['temperatureK'])}, {{air=1.0}}
gm:updateThermoFromPT(gs)
gm:updateSoundSpeed(gs)
gm:updateTransCoeffs(gs)
{flow}
local inflow = FlowState:new{{p=gs.p, T=gs.T, massf=gs.massf,
    velx=ux, vely=uy, velz=uz, tke=tke, omega=omega}}
-- 原文件不复制；此开关是固定v5公开的Lua参数，显式逆转项目SU2楔体序。
registerFluidGrid{{grid=UnstructuredGrid:new{{filename={lua_string(mesh_native)},
    fmt="su2text", scale=1.0, expect_gmsh_order_for_wedges=true}}, fsTag="initial"}}
local bcDict = {{
    [{lua_string(c['boundaryMarkers']['wall'])}] = WallBC_NoSlip_FixedT:new{{
        Twall={f(c['wallTemperatureK'])}, wall_function=false, group="wall"}},
    [{lua_string(c['boundaryMarkers']['farfield'])}] = InOutFlowBC_Ambient:new{{flowState=inflow}}
}}
makeFluidBlocks(bcDict, {{initial=inflow}})
config.boundary_groups_for_loads = "wall"
config.flux_calculator = "adaptive_hanel_ausmdv"
config.interpolation_order = 2
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
-- 步数只是未来任务上限；本准备器无运行命令，无自动接续。
NewtonKrylovGlobalConfig{{number_of_phases=1,
    max_newton_steps={c['maxNewtonSteps']}, stop_on_relative_residual={f(c['stopRelativeResidual'])},
    max_consecutive_bad_steps=2, use_preconditioner=true, preconditioner="ilu", ilu_fill=0,
    max_linear_solver_iterations=50, max_linear_solver_restarts=0,
    use_line_search=true, use_physicality_check=true,
    write_loads=true, steps_between_loads_update={max(1,c['maxNewtonSteps'])},
    total_snapshots=2, steps_between_snapshots={max(1,c['maxNewtonSteps'])},
    write_snapshot_on_last_step=true, write_loads_on_last_step=true,
    steps_between_status=1, steps_between_diagnostics=1}}
NewtonKrylovPhase:new{{use_local_timestep=true, use_auto_cfl=false,
    start_cfl=0.1, max_cfl=0.1, ignore_stopping_criteria=false,
    residual_interpolation_order=2, jacobian_interpolation_order=1,
    linear_solve_tolerance=0.01}}
'''


def prepare_case(mesh_path, output_dir, condition, *, source_order="rah-su2-v1"):
    started = time.monotonic();source = Path(mesh_path).resolve();target = Path(output_dir).resolve()
    require(source.is_file(), "源网格不存在")
    require(not target.exists(), "目标目录已存在，禁止覆盖")
    require(source_order == "rah-su2-v1", "必须显式声明项目SU2楔体/边界序；其他来源未实现")
    c = validate_condition(condition);before = sha(source)
    points, groups, markers = read_mesh(source)
    require(set(markers) == set(c["boundaryMarkers"].values()), "源边界与显式wall/farfield映射不完全一致")
    geometry = audit_order(points, groups)
    require(sha(source) == before, "审核期间源网格发生变化")
    mesh_native = native_path(source)
    wall, farfield = (markers[c["boundaryMarkers"][name]] for name in ("wall", "farfield"))
    mesh = dict(nodes=len(points), cells=sum(len(v) for v in groups.values()),
                cellTypeCounts={str(k): len(v) for k, v in groups.items()},
                wallFaces=wall["faces"], farfieldFaces=farfield["faces"], wallAreaM2=wall["areaM2"],
                sourceSHA256=before, sourcePath=str(source), nativeSourcePath=mesh_native,
                sourceOrder=source_order, expectGmshOrderForWedges=True,
                wedgePermutationZeroBased=WEDGE_PERMUTATION, markers=markers, **geometry)
    outputs = {
        "gas.inp": '-- 使用固定官方数据库；不改物性核心或私设60K系数。\nmodel = "ThermallyPerfectGas"\nspecies = {"air"}\n',
        "job.lua": make_job(mesh_native, c),
    }
    notes = """# Eilmer 完整域输入：仅准备，原生未验证

固定 Eilmer v5.0.0 / f53f4609a0331d48efee69a4e4f3c3598378cc03。
本目录不包含求解结果；不表示原生准备/收敛/热流/实验验证已通过。
引用原 SU2 文件，禁止改动；没有复制大网格。实际 WSL 路径尚须验证。
原生读取必须保留 expect_gmsh_order_for_wedges=true。边界标签和全部原壁面保留。
使用 Wilcox k-omega，不能称为 SST，不能接续旧 SST 保存场。
官方伪空气热力表 200–6000 K；表外采用端点定cp延拓h/s，不裁温，表外物理资格未知。
不能把有限值或数值收敛当成全面积物理求解资格。RANS近壁/壁热和分离区须独立验证。
loads 的 q_total=q_cond+q_diff 以内部面法向定义，必须连同outsign核对入壁正号；
原生非结构y+使用L_min/2估距，不能直接充当项目严格法向距离验收。

安装/编译与系统变更另行处理，本准备器没有启动它们。固定运行时就绪且输入SHA复核后，
官方命令顺序是 prep-gas gas.inp rah-air.gas、lmr prep-grid、lmr prep-sim --job=job.lua，
然后才可在已登记预算与监督下执行 lmr run。本轮这些命令全部未执行。
原始面顺序的结果映射、q符号/面积积分与稀疏输出读取均待独立实现，不能直接装载研究结果。
"""
    outputs["PREPARATION.md"] = notes
    manifest = dict(schema="RAH.EilmerPreparedCase.v1", status="prepared-only", nativeVerified=False,
        runtimeReady=False, productionQualified=False, newQualifiedAreaM2=0, nativeTasks=0,
        meshCopied=False, sourceStable=True, sourceCommit=COMMIT, sourceTag=TAG,
        turbulenceModel="Wilcox k-omega (2006); not SST", condition=c, mesh=mesh,
        thermalProperties={"model": "official ThermallyPerfectGas species air", "tableMinimumK": 200,
            "tableMaximumK": 6000, "outsideBehavior": "endpoint Cp with continuous enthalpy/entropy extension; temperature not clipped",
            "givenTemperaturesOutsideTable": any(not 200 <= c[k] <= 6000 for k in ("temperatureK", "wallTemperatureK")),
            "computedFieldTemperatureRangeKnown": False, "physicalQualification": False},
        pending=["WSL/runtime and pinned native preparation", "native geometric interpretation/topology",
                 "native gas/velocity/turbulence receipt", "signed heat and original-face mapping",
                 "near-wall, conservation, convergence and public physical validation"],
        sourceURLs={name: f"https://github.com/gdtk-uq/gdtk/blob/{COMMIT}/{name}" for name in
            ["src/geom/luawrap/luausgrid.d", "src/geom/grid/usgrid.d", "src/lmr/lua-modules/bc.lua",
             "src/lmr/lua-modules/nkconfig.lua", "src/lmr/turbulence.d", "src/lmr/loads.d",
             "src/gas/species-database/air.lua", "src/gas/thermo/cea_thermo_curves.d"]})
    target.mkdir(parents=True, exist_ok=False)
    for name, content in outputs.items():
        (target / name).write_text(content, encoding="utf-8", newline="\n")
    manifest["filesSHA256"] = {name: sha(target / name) for name in outputs}
    manifest["elapsedSeconds"] = time.monotonic() - started
    (target / "case_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", required=True);parser.add_argument("--output", required=True)
    parser.add_argument("--condition-json", required=True)
    parser.add_argument("--source-order", default="rah-su2-v1")
    args = parser.parse_args()
    condition = json.loads(Path(args.condition_json).read_text(encoding="utf-8-sig"))
    result = prepare_case(args.mesh, args.output, condition, source_order=args.source_order)
    print(json.dumps({"status": result["status"], "nativeVerified": False,
                      "manifest": str(Path(args.output).resolve() / "case_manifest.json")}))


if __name__ == "__main__":
    main()
