"""同气体显式回调与旧默认回归；不运行任何CFD。"""
from pathlib import Path
import json,sys,time,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from nistIntegratedAir import evaluate,R,MODEL
from openFoamAirThermo import evaluate as old_evaluate
from openFoamFields import audit_state
from openFoamConservation import geometry,flux_and_residual,reference_scales
from hybridEulerReference import control_reference
from prepareOpenFoamControl import save


def main():
    start=time.monotonic();base=ROOT/'outputs/v202_validation';out=base/'callback_check01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py');rejections=[]
    def reject(name,action):
        try:action()
        except (ValueError,AssertionError):rejections.append(name);return
        raise AssertionError('应拒绝但通过 '+name)
    gas=evaluate(1000.,180.);fields={k:np.array([v]) for k,v in [('p',1000.),('T',180.),('rho',gas['densityKgm3']),('e',gas['internalEnergyJkg'])]}
    fields['U']=np.zeros((1,3));audit_state(fields,evaluate)
    reject('new_low_state_with_old_default',lambda:audit_state(fields))
    for name,key,value in [('pressure_limit','p',100001.),('temperature_limit','T',99.9),('wrong_energy','e',0.),('wrong_density','rho',1.)]:
        changed={k:v.copy() for k,v in fields.items()};changed[key][:]=value
        reject(name,lambda f=changed:audit_state(f,evaluate))
    case=ROOT/'outputs/v199_validation/small01';g=geometry(case)
    f,flux,res=flux_and_residual(case/'4e-07',g)
    f2,flux2,res2=flux_and_residual(case/'4e-07',g,old_evaluate,R)
    assert np.array_equal(flux,flux2) and np.array_equal(res,res2)
    reject('old_saved_energy_with_new_model',lambda:flux_and_residual(case/'4e-07',g,evaluate,R))
    # 模块二：默认路径与v200保存审核数值逐项核对，不重跑历史CFD。
    oldcase=ROOT/'outputs/v200_validation/restart01';gg=geometry(oldcase)
    _,_,rr=flux_and_residual(oldcase/'0.003',gg);scale,speed=reference_scales()
    norms=np.sum(abs(rr)*gg['volume'][:,None],axis=0)/gg['volume'].sum()/scale/speed
    prior=json.loads((ROOT/'outputs/v200_validation/restart_audit02/report.json').read_text(encoding='utf-8'))['frames'][-1]['residualVolumeL1Scaled']
    assert np.array_equal(norms,prior)
    save(out/'report.json',dict(passed=True,rejectionCount=len(rejections),rejections=rejections,oldResidualBitwiseEqual=True,
        explicitDefaultEqual=True,elapsedSeconds=time.monotonic()-start,CFDRuns=0))
    refout=base/'reference01';refout.mkdir(exist_ok=False)
    ref=control_reference(evaluate,R,MODEL,upstream_pressure=1000,upstream_temperature=180,temperature_min=100,entropy_pressure=100000)
    assert max(ref['residuals'][k] for k in ['massRelative','normalMomentumRelative','totalEnthalpyRelative'])<1e-10
    assert ref['residuals']['entropyExpansionAbsolute']<1e-7 and ref['residuals']['expansionAngleRad']<1e-10
    assert ref['states'][2]['temperatureK']<200 and ref['states'][1]['temperatureK']>200
    ref['passed']=True;save(refout/'report.json',ref)
    print(json.dumps(dict(passed=True,rejections=len(rejections),referenceTemperatures=[r['temperatureK'] for r in ref['states']],elapsedSeconds=time.monotonic()-start)))


if __name__=='__main__':main()
