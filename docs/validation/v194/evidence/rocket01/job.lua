-- 模块一：原生CEA空气和真实20度来流；壁温仅为后续热目标，不施加给Euler。
config.dimensions=3
config.solver_mode="time-marching"
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
config.max_step=0
config.max_attempts_for_step=1
config.print_count=1
config.write_transient_residuals=true
config.dt_plot_value=1.0
