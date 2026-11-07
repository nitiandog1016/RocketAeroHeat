-- 仅生成，未原生验证。局部伪时间不代表飞行时间。
config.dimensions = 3
config.solver_mode = "time-marching"
config.field_format = "rawbinary"
config.grid_format = "rawbinary"
config.viscous = false
config.reacting = false
config.turbulence_model = "none"
config.adjust_invalid_cell_data = false
config.max_invalid_cells = 0
local nsp,nmodes,gm=setGasModel("rah-air.gas")
local gs = GasState:new{gm}
gs.p=5000.0; gs.T=300.0; gs.massf={air=1.0}
gm:updateThermoFromPT(gs); gm:updateSoundSpeed(gs); gm:updateTransCoeffs(gs)
local inflow=FlowState:new{p=gs.p,T=gs.T,massf=gs.massf,velx=-3.0*gs.a,vely=0.0,velz=0.0}
local receipt=assert(io.open("native-input-state.json","w"))
receipt:write(string.format('{"p":%.17g,"T":%.17g,"rho":%.17g,"e":%.17g,"a":%.17g,"cp":%.17g,"velx":%.17g}',
 gs.p,gs.T,gs.rho,gs.u,gs.a,gm:Cp(gs),-3.0*gs.a)); receipt:close()
makeFluidBlocks({wall=WallBC_WithSlip:new{group="wall"},span=WallBC_WithSlip:new{},
 outer=InFlowBC_Supersonic:new{flowCondition=inflow},
 inflow=InFlowBC_Supersonic:new{flowCondition=inflow},outflow=OutFlowBC_Simple:new{}},
 {initial=inflow})
config.flux_calculator="ausmdv"
config.interpolation_order=2
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
config.gasdynamic_update_scheme="tvd_rk3"
config.with_local_time_stepping=true
config.local_time_stepping_limit_factor=1000
config.with_super_time_stepping=false
config.residual_smoothing=false
config.stringent_cfl=true
config.cfl_value=0.3
config.cfl_count=1
config.fixed_time_step=false
config.dt_init=1.0e-10
config.dt_max=1.0e-3
config.max_time=1.0
config.max_step=400
config.max_attempts_for_step=1
config.print_count=1
config.write_transient_residuals=true
config.dt_plot_value=1.0
config.write_flow_solution_at_step=200
