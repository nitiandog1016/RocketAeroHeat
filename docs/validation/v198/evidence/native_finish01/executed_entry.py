"""闭合实测原生常数差：复用已编译查询器和原网格，只修正显式EOS参数。"""
from pathlib import Path
import json,shutil,time
import numpy as np
from prepareOpenFoamControl import ROOT,BASE,sha,save,header
from openFoamAirThermo import mixture_dictionary,contract,evaluate,R,HF,FOAM_R_UNIVERSAL,MIX_MOLAR_MASS_GMOL


def main():
    """模块一：旧失败留存；新配置通过独立原生查询后再检查原生三维网格。"""
    from runOpenFoamProcess import run
    start=time.monotonic();out=BASE/'native_finish01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py')
    old=BASE/'native_entry02';before=(old/'query/stdout.txt').read_text()
    observed=float(next(line.split()[1] for line in before.splitlines() if line.startswith('R ')))
    assert abs(observed*MIX_MOLAR_MASS_GMOL/1000/FOAM_R_UNIVERSAL-1)<2e-15
    for rel in ['src/OpenFOAM/global/constants/fundamental/fundamentalConstants.C',
        'src/OpenFOAM/global/constants/physicoChemical/physicoChemicalConstants.C']:
        p=Path('/opt/openfoam13')/rel;dst=out/'native_constants'/p.name;dst.parent.mkdir(exist_ok=True);shutil.copyfile(p,dst)
    (out/'mixture.dict').write_text(mixture_dictionary())
    binary=old/'probe/rahThermoProbe'
    prior=json.loads((old/'query/report.json').read_text());assert sha(binary)==prior['inputSHA256'][str(binary)]
    q=run([binary,out/'mixture.dict'],out,out/'query',20,[__file__,binary,out/'mixture.dict',ROOT/'scripts/openFoamAirThermo.py'],32)
    assert q['passed']
    text=(out/'query/stdout.txt').read_text();assert 'out of temperature range' not in text
    rows=np.array([[float(x) for x in line.split()[1:]] for line in text.splitlines() if line.startswith('ROW ')])
    meta=next(line.split() for line in text.splitlines() if line.startswith('R '))
    assert len(rows)==9 and abs(float(meta[1])/R-1)<2e-15 and abs(float(meta[3])-HF)<2e-9
    fields=['cpJkgK','cvJkgK','enthalpyJkg','sensibleEnthalpyJkg','internalEnergyJkg','densityKgm3','soundSpeedMps']
    g=evaluate(5000,rows[:,0]);errors={}
    for i,key in enumerate(fields):
        errors[key]=float(np.max(abs(rows[:,i+1]-g[key])/np.maximum(1,abs(g[key]))))
        assert errors[key]<2e-10,key
    inv=float(np.max(abs(rows[:,-1]-rows[:,0])));assert inv<2e-7
    case=BASE/'control02';shutil.copytree(BASE/'control01',case)
    (case/'constant/physicalProperties').write_text(header('physicalProperties')+
        'thermoType { type hePsiThermo; mixture pureMixture; transport const; thermo janaf; equationOfState perfectGas; specie specie; energy sensibleInternalEnergy; }\nmixture {\n'+mixture_dictionary()+'}\n')
    save(case/'thermo_contract.json',contract())
    for p in (case/'constant/polyMesh').iterdir():assert sha(p)==sha(BASE/'control01/constant/polyMesh'/p.name)
    mesh=run(['checkMesh','-allTopology','-allGeometry'],case,out/'mesh',60,
        [__file__,*[p for p in (case/'constant/polyMesh').iterdir()]],64)
    log=(out/'mesh/stdout.txt').read_text();assert mesh['passed'] and 'Mesh OK.' in log
    result=dict(passed=True,propertyMaxRelativeErrors=errors,inverseTemperatureMaxAbsK=inv,
        query=q,mesh=mesh,contract=contract(),nativeRows=rows.tolist(),
        originalMolarParameterRRelativeDifference=observed/R-1,engineBuilds=0,
        originalNativeFailureRetained=True,inputSHA256={str(p):sha(p) for p in case.rglob('*') if p.is_file()},
        elapsedSeconds=time.monotonic()-start,physicalQualification=False,CFDRuns=0)
    save(out/'report.json',result);print(json.dumps({k:result[k] for k in ['passed','elapsedSeconds','propertyMaxRelativeErrors','inverseTemperatureMaxAbsK']}))

if __name__=='__main__':main()
