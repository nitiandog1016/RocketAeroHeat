-- 模块三：仅原生算子控制，禁止牛顿/时间发展；人工层顶不是物理终端资格。
config.dimensions=3
config.solver_mode="steady"
config.field_format="rawbinary"
config.grid_format="rawbinary"
config.viscous=true
config.reacting=false
config.turbulence_model="spalart_allmaras_edwards"
config.adjust_invalid_cell_data=false
config.max_invalid_cells=0
setGasModel("phase-air.lua")
local initial=FlowState:new{p=2500.0,T=300.0,velx=-1000.0,vely=0.0,velz=200.0,nuhat=1.0e-5}
makeFluidBlocks({wall=WallBC_NoSlip_FixedT:new{Twall=300.0,group="wall"},
 outer=InFlowBC_Supersonic:new{flowCondition=initial},
 exchange=InFlowBC_Supersonic:new{flowCondition=initial,label="exchange"}}, {initial=initial})
config.flux_calculator="ausmdv"
config.interpolation_order=2
config.thermo_interpolator="rhou"
config.extrema_clipping=false
config.unstructured_limiter="hvenkat"
config.smooth_limiter_coeff=1.0
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
config.viscous_least_squares_type="weighted_qr"
NewtonKrylovGlobalConfig{number_of_phases=1,max_newton_steps=0,use_preconditioner=true,
 preconditioner="ilu",ilu_fill=0,max_linear_solver_iterations=50,max_linear_solver_restarts=0,
 use_physicality_check=true,use_line_search=true,stop_on_relative_residual=1.0e-8}
NewtonKrylovPhase:new{residual_interpolation_order=2,jacobian_interpolation_order=1,
 use_residual_smoothing=false,start_cfl=0.1}
