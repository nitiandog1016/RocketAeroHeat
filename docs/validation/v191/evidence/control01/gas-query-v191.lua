-- 独立原生查询：输出实际CEA焓/熵参考，不人为平移。
local gm=GasModel:new{"rah-air.gas"}
local gs=GasState:new{gm};gs.massf={air=1.0}
local states={{5000,200},
{5000,250},
{5000,300},
{5000,400},
{5000,700},
{5000,799},
{5000,800},
{5000,900},
{5000,1000},
{5000,1100},
{5000,1200},
{5000,1201},
{5000,2500},
{5000,6000},
{5000,300},
{10268.789830080517,372.26527282381682},
{5012.280709824975,303.4472812396429},
{101325,298.14999999999998},
{25000,1000}}
local f=assert(io.open("gas-query-v191.csv","w"))
f:write("p,T,rho,e,h,s,a,cp,mu,k\n")
for _,v in ipairs(states) do
 gs.p=v[1];gs.T=v[2]
 gm:updateThermoFromPT(gs);gm:updateSoundSpeed(gs);gm:updateTransCoeffs(gs)
 f:write(string.format("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
 gs.p,gs.T,gs.rho,gs.u,gm:enthalpy(gs),gm:entropy(gs),gs.a,gm:Cp(gs),gs.mu,gs.k))
end
f:close()
