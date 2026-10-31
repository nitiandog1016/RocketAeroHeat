"""真实初态来源核对和拒绝反例；不改变历史查询、初态或实际失败记录。"""
from pathlib import Path
import hashlib,json,sys,copy
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from su2ThermalRuntime import validate_coolprop_manifest,validate_coolprop_initialization,resolve
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    folder=ROOT/'outputs/v122_validation/actual04';out=ROOT/'outputs/v122_validation/initialization_contracts_v3.json';assert not out.exists()
    data=read(folder/'case_manifest.json');config=(folder/'case.cfg').read_text()
    _,rt,flavor,_,locks=resolve(folder,config);assert flavor=='rah-coolprop-recon-v3'
    assert str(ROOT/'outputs/v122_validation/linked_property_equivalence/report.json') in locks
    source=ROOT/'outputs/v120_validation/study01/control';old=np.loadtxt(str(source/'restart_flow.csv'),delimiter=',',skiprows=1)
    new=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1)
    props=np.genfromtxt(str(folder/'initial_thermodynamics/DE_output.csv'),names=True,delimiter=',')
    v=old[:,5:8]/old[:,4,None];internal=old[:,8]/old[:,4]-.5*np.sum(v*v,axis=1)
    p=.4*old[:,4]*internal;T=internal*.4/287.05287
    errors=dict(pressure=float(max(abs(props['p_Pa']/p-1))),temperature=float(max(abs(props['T_K']/T-1))),
        velocity=float(np.max(abs(new[:,5:8]/new[:,4,None]-v)/np.maximum(abs(v),1))))
    assert max(errors.values())<1e-8 and np.array_equal(old[:,:4],new[:,:4])
    assert data['initialization']['sourceCumulativeIteration']==799 and data['initialization']['newPhysicsLocalIterationStartsAt']==0
    cases={}
    for name,mutation in [('old_energy_label',lambda d:d['initialization'].update(kind='same-physics-continuation')),
                          ('wrong_eos_label',lambda d:d.update(thermodynamicModel='ideal-gas')),
                          ('stale_initial_certificate',lambda d:d.update(coolPropInitializationSHA256='0'*64))]:
        bad=copy.deepcopy(data);mutation(bad)
        try:
            validate_coolprop_manifest(bad,config)
            validate_coolprop_initialization(folder,bad,ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-recon-v3/runtime.json')
        except AssertionError:cases[name]=True
        else:raise AssertionError(name)
    for name,text in [('ordinary_farfield',config+'\nMARKER_FAR= (farfield)\n'),('unadapted_rans',config.replace('KIND_TURB_MODEL= NONE','KIND_TURB_MODEL= SST')),
                      ('binary_old_energy',config.replace('READ_BINARY_RESTART= NO','READ_BINARY_RESTART= YES'))]:
        try:validate_coolprop_manifest(data,text)
        except AssertionError:cases[name]=True
        else:raise AssertionError(name)
    report=dict(passed=True,rejections=cases,independentPrimitiveErrors=errors,nodeCount=len(new),sourceCumulativeIteration=799,
        newPhysicsStartsAt=0,strictConservativeRemapClaimed=False,actualFieldProduced=True,propertyRebindingProofLocked=True,physicsQualified=False,
        sources={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/su2ThermalRuntime.py',folder/'case_manifest.json',folder/'initial_flow.csv',source/'restart_flow.csv']})
    out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
