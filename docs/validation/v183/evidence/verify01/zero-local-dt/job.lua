-- RocketAeroHeat：仅准备的Eilmer v5.0.0输入，原生尚未验证。
-- 固定源码 f53f4609a0331d48efee69a4e4f3c3598378cc03；Wilcox k-omega，非SST、非旧重启动。
config.dimensions = 3
config.axisymmetric = false
config.solver_mode = "time-marching"
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
-- 按同一原生热完美气体求声速/输运，不能使用定gamma替代。
local speed = 5 * gs.a
local alpha, beta = math.rad(20), math.rad(0)
local ux = -speed*math.cos(alpha)*math.cos(beta)
local uy = -speed*math.sin(beta)
local uz = speed*math.sin(alpha)*math.cos(beta)
local tke = 1.5*(0.01*speed)^2
local omega = gs.rho*tke/(gs.mu*10)
-- 保存实际原生物性及来流，供独立审核；给定壁温与流体温度分开。
local receipt = assert(io.open("native-input-state.json", "w"))
receipt:write(string.format('{"p":%.17g,"T":%.17g,"rho":%.17g,"a":%.17g,"u":%.17g,"cp":%.17g,"mu":%.17g,"k":%.17g,"velx":%.17g,"vely":%.17g,"velz":%.17g,"tke":%.17g,"omega":%.17g,"Tw_given":%.17g}',
    gs.p, gs.T, gs.rho, gs.a, gs.u, gm:Cp(gs), gs.mu, gs.k, ux, uy, uz, tke, omega, 300))
receipt:close()
local inflow = FlowState:new{p=gs.p, T=gs.T, massf=gs.massf,
    velx=ux, vely=uy, velz=uz, tke=tke, omega=omega}
-- 网格已由grid.lua独立准备，流场准备时不重复读取/构造大SU2对象。
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
-- 显式TVD RK3与局部伪时间；此时间不代表真实轨迹时间。
config.solver_mode = "time-marching"
config.gasdynamic_update_scheme = "tvd_rk3"
config.with_local_time_stepping = true
config.local_time_stepping_limit_factor = 1000
config.with_super_time_stepping = false
config.residual_smoothing = false
config.stringent_cfl = true
config.cfl_value = 0.1
config.cfl_count = 1
config.viscous_signal_factor = 1.0
config.turbulent_signal_factor = 1.0
config.fixed_time_step = false
config.dt_init = 1e-10
config.dt_max = 0.001
config.max_time = 1.0
config.max_step = 6
config.max_attempts_for_step = 1
config.print_count = 1
config.write_transient_residuals = true
config.write_loads = true
config.dt_plot_value = 1.0
config.dt_loads = 1.0
config.write_loads_at_step = 3
config.write_flow_solution_at_step = 3
