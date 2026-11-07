# -*- coding: utf-8 -*-
"""用准确原生Lua气体API导出锚点，供独立Python回调核对；无CFD推进。"""
from pathlib import Path
import argparse
import json
import shutil
from runEilmerProcess import run, sha, PREFIX


def query(case, reference, record):
    case,reference,record=map(Path,(case,reference,record))
    assert not record.exists()
    refs=json.loads(reference.read_text())['states']
    states=[(5000.,t) for t in [200,250,300,400,700,799,800,900,1000,1100,1200,1201,2500,6000]]
    states += [(s['pressurePa'],s['temperatureK']) for s in refs]
    states += [(101325.,298.15),(25000.,1000)]
    job=case/'gas-query-v191.lua';assert not job.exists()
    rows=',\n'.join('{%.17g,%.17g}'%row for row in states)
    text='''-- 独立原生查询：输出实际CEA焓/熵参考，不人为平移。
local gm=GasModel:new{"rah-air.gas"}
local gs=GasState:new{gm};gs.massf={air=1.0}
local states={STATES}
local f=assert(io.open("gas-query-v191.csv","w"))
f:write("p,T,rho,e,h,s,a,cp,mu,k\\n")
for _,v in ipairs(states) do
 gs.p=v[1];gs.T=v[2]
 gm:updateThermoFromPT(gs);gm:updateSoundSpeed(gs);gm:updateTransCoeffs(gs)
 f:write(string.format("%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\\n",
 gs.p,gs.T,gs.rho,gs.u,gm:enthalpy(gs),gm:entropy(gs),gs.a,gm:Cp(gs),gs.mu,gs.k))
end
f:close()
'''.replace('STATES',rows)
    job.write_text(text,encoding='utf8')
    locks=[job,case/'rah-air.gas',reference,Path(__file__)]+[PREFIX/'source'/p for p in
        ['src/gas/composite_gas.d','src/gas/thermo/therm_perf_gas_mix.d',
         'src/gas/thermo/cea_thermo_curves.d','src/gas/species-database/air.lua',
         'src/gas/physical_constants.d','src/gas/luagas_model.d']]
    result=run(['lmr','custom-script','--job='+job.name],case,record,20,128,locks,4096,4096)
    shutil.copyfile(__file__,record/'executed_query.py')
    assert result['passed'] and (case/'gas-query-v191.csv').is_file()
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--reference',required=True)
    p.add_argument('--record',required=True);a=p.parse_args()
    r=query(a.case,a.reference,a.record);print(json.dumps(dict(passed=r['passed'],elapsedSeconds=r['elapsedSeconds'])))
