-- 模块一：原生CEA空气和真实20度来流；壁温仅为后续热目标，不施加给Euler。
config.dimensions=3
config.solver_mode="steady"
config.field_format="rawbinary"
config.grid_format="rawbinary"
config.viscous=false
config.reacting=false
config.turbulence_model="none"
config.adjust_invalid_cell_data=false
config.max_invalid_cells=0
local nsp,nmodes,gm=setGasModel("rah-air.gas")
local gs=GasState:new{gm}
gs.p=1197.0279007383137; gs.T=226.50908361133003; gs.massf={air=1.0}
gm:updateThermoFromPT(gs); gm:updateSoundSpeed(gs); gm:updateTransCoeffs(gs)
local speed=5.0*gs.a
local alpha,beta=math.rad(20),math.rad(0)
local ux=-speed*math.cos(alpha)*math.cos(beta)
local uy=-speed*math.sin(beta)
local uz=speed*math.sin(alpha)*math.cos(beta)
local inflow=FlowState:new{p=gs.p,T=gs.T,massf=gs.massf,velx=ux,vely=uy,velz=uz}
local receipt=assert(io.open("native-input-state.json","w"))
receipt:write(string.format('{"p":%.17g,"T":%.17g,"rho":%.17g,"e":%.17g,"a":%.17g,"cp":%.17g,"velx":%.17g,"vely":%.17g,"velz":%.17g}',
 gs.p,gs.T,gs.rho,gs.u,gs.a,gm:Cp(gs),ux,uy,uz)); receipt:close()
-- 模块二：成熟反射滑移壁及环境外边界；无热壁/湍流/局部截断。
makeFluidBlocks({wall=WallBC_WithSlip:new{group="wall"},
 farfield=InOutFlowBC_Ambient:new{flowState=inflow}},{initial=inflow})
-- 模块三：保持已核对的Euler算子；本轮max_step=0，准备不授予推进权限。
config.flux_calculator="ausmdv"
config.interpolation_order=2
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
-- 模块三：原生全表面谱半径、物理性检查和残差回溯；复数导数。
config.extrema_clipping=false
config.thermo_interpolator="rhop"
config.unstructured_limiter="hvenkat"
config.smooth_limiter_coeff=1.0
config.apply_unstructured_limiter_stagnation_point_filter=true
NewtonKrylovGlobalConfig{
 number_of_phases=1, use_preconditioner=true, preconditioner_perturbation=1.0e-50,
 preconditioner="ilu", ilu_fill=0, max_newton_steps=120, max_consecutive_bad_steps=10,
 stop_on_relative_residual=1.0e-8, frechet_derivative_perturbation=1.0e-50,
 use_scaling=true, max_linear_solver_iterations=50, max_linear_solver_restarts=0,
 inviscid_cfl_only=true, use_line_search=true, line_search_order=3,
 use_scaling_in_line_search=true, use_physicality_check=true,
 allowable_relative_mass_change=0.5, min_relaxation_factor_for_update=0.1,
 min_relaxation_factor_for_cfl_growth=0.5,
 number_of_steps_for_setting_reference_residuals=0, steps_between_status=1,
 write_loads=false, total_snapshots=8, steps_between_snapshots=20, steps_between_diagnostics=1
}
NewtonKrylovPhase:new{
 frozen_preconditioner=true, use_adaptive_preconditioner=true, steps_between_preconditioner_update=5,
 linear_solve_tolerance=0.01, residual_interpolation_order=2, jacobian_interpolation_order=1,
 use_residual_smoothing=true, use_auto_cfl=true, use_local_timestep=true,
 threshold_relative_residual_for_cfl_growth=0.99, start_cfl=0.1, max_cfl=100.0,
 auto_cfl_exponent=1.0, limit_on_cfl_decrease_ratio=1.0
}
