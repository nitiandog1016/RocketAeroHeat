-- RocketAeroHeat：仅准备的Eilmer v5.0.0输入，原生尚未验证。
-- 固定源码 f53f4609a0331d48efee69a4e4f3c3598378cc03；Wilcox k-omega，非SST、非旧重启动。
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
local gs = GasState:new{gm}
gs.p, gs.T, gs.massf = 1197.0279007383137, 226.50908361133003, {air=1.0}
gm:updateThermoFromPT(gs)
gm:updateSoundSpeed(gs)
gm:updateTransCoeffs(gs)
-- 显式速度模式；未核对它与任何目标Mach相等。
local ux, uy, uz = -1400, 500, 0
local tke, omega = 2, 120
local inflow = FlowState:new{p=gs.p, T=gs.T, massf=gs.massf,
    velx=ux, vely=uy, velz=uz, tke=tke, omega=omega}
-- 原文件不复制；此开关是固定v5公开的Lua参数，显式逆转项目SU2楔体序。
registerFluidGrid{grid=UnstructuredGrid:new{filename="/mnt/e/RocketAeroHeat/outputs/v178_validation/fixture01/mixed_original.su2",
    fmt="su2text", scale=1.0, expect_gmsh_order_for_wedges=true}, fsTag="initial"}
local bcDict = {
    ["wall"] = WallBC_NoSlip_FixedT:new{
        Twall=300, wall_function=false, group="wall"},
    ["farfield"] = InOutFlowBC_Ambient:new{flowState=inflow}
}
makeFluidBlocks(bcDict, {initial=inflow})
config.boundary_groups_for_loads = "wall"
config.flux_calculator = "adaptive_hanel_ausmdv"
config.interpolation_order = 2
config.spatial_deriv_calc = "least_squares"
config.spatial_deriv_locn = "cells"
-- 步数只是未来任务上限；本准备器无运行命令，无自动接续。
NewtonKrylovGlobalConfig{number_of_phases=1,
    max_newton_steps=20, stop_on_relative_residual=1e-08,
    max_consecutive_bad_steps=2, use_preconditioner=true, preconditioner="ilu", ilu_fill=0,
    max_linear_solver_iterations=50, max_linear_solver_restarts=0,
    use_line_search=true, use_physicality_check=true,
    write_loads=true, steps_between_loads_update=20,
    total_snapshots=2, steps_between_snapshots=20,
    write_snapshot_on_last_step=true, write_loads_on_last_step=true,
    steps_between_status=1, steps_between_diagnostics=1}
NewtonKrylovPhase:new{use_local_timestep=true, use_auto_cfl=false,
    start_cfl=0.1, max_cfl=0.1, ignore_stopping_criteria=false,
    residual_interpolation_order=2, jacobian_interpolation_order=1,
    linear_solve_tolerance=0.01}
