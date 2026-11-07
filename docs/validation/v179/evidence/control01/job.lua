-- 具名理想N2层流控制，忽略振动是明示模型假设，不替换真实箭体变cp模型。
config.dimensions=2
config.axisymmetric=true
config.solver_mode="steady"
config.field_format="rawbinary"
config.grid_format="gziptext"
config.viscous=true
config.turbulence_model="none"
config.reacting=false
config.adjust_invalid_cell_data=false
config.max_invalid_cells=0
config.flux_calculator="adaptive_hanel_ausmdv"
config.interpolation_order=2
config.spatial_deriv_calc="least_squares"
config.spatial_deriv_locn="cells"
local _,_,gm=setGasModel("nitrogen.gas")
local gs=GasState:new{gm}
local f=assert(io.open("native-gas.csv","w"))
f:write("T,p,rho,a,e,cp,mu,k\n")
for _,T in ipairs{100,120.4,150,200,295.2,500,1000,2000} do
 gs.p=31.88;gs.T=T;gs.massf={N2=1}
 gm:updateThermoFromPT(gs);gm:updateSoundSpeed(gs);gm:updateTransCoeffs(gs)
 f:write(string.format("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",T,gs.p,gs.rho,gs.a,gs.u,gm:Cp(gs),gs.mu,gs.k))
end
f:close()
local inflow=FlowState:new{p=31.88,T=120.4,velx=2304,massf={N2=1}}
local initial=FlowState:new{p=31.88/5,T=120.4,massf={N2=1}}
makeFluidBlocks({inflow=InFlowBC_Supersonic:new{flowState=inflow},
 wall=WallBC_NoSlip_FixedT:new{Twall=295.2,group="wall"},
 outflow=OutFlowBC_FixedP:new{p_outside=31.88/5}}, {inflow=inflow,initial=initial})
config.boundary_groups_for_loads="wall"
-- 原生两阶段精度推进：50步一阶启动后转二阶；自动CFL按残差反馈，不扫描。
NewtonKrylovGlobalConfig{number_of_phases=2,max_steps_in_initial_phases={50},
 max_newton_steps=2000,stop_on_relative_residual=1e-8,
 min_cfl=0.01,max_cfl=1000,use_line_search=true,use_physicality_check=true,
 max_consecutive_bad_steps=2,use_preconditioner=true,preconditioner="ilu",ilu_fill=0,
 max_linear_solver_iterations=50,max_linear_solver_restarts=0,
 write_loads=true,steps_between_loads_update=50,steps_between_snapshots=250,total_snapshots=3,
 write_snapshot_on_last_step=true,write_loads_on_last_step=true,
 steps_between_status=10,steps_between_diagnostics=1}
NewtonKrylovPhase:new{use_auto_cfl=true,start_cfl=0.1,max_cfl=1000,
 threshold_relative_residual_for_cfl_growth=0.9,auto_cfl_exponent=0.75,
 residual_interpolation_order=1,jacobian_interpolation_order=1,linear_solve_tolerance=0.01}
NewtonKrylovPhase:new{start_cfl=1,residual_interpolation_order=2,jacobian_interpolation_order=1}
