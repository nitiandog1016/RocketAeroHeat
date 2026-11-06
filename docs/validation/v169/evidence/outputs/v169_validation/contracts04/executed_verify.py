"""模块用途：独立局部归一化、入口语义拒绝和旧合同回归，全部只读。"""
import argparse
from pathlib import Path
import sys
import subprocess
import types
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from naca3344Model import D,DATA,source_data
from naca3344Normalization import published_recovery,compare_local,require_experimental_qualification
from su2NACALaminarInterface import inlet_state_check
from su2ThermalRuntime import read,sha,resolve
from validateIdealAir import save


def main(out):
    out.mkdir(parents=True,exist_ok=False)
    (out/'executed_verify.py').write_bytes(Path(__file__).read_bytes())
    locks={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/naca3344Normalization.py',
        ROOT/'scripts/su2IdealAirInterface.py',ROOT/'scripts/su2IdealAirCheckpoints.py',
        ROOT/'scripts/su2NACALaminarInterface.py',DATA/'normalization/manifest.json']}
    s=np.array([.1,.3,.6,.9,1.04])*D
    coefficient=published_recovery(s)
    contract,_,_=source_data();T0=contract['freeStream']['totalTemperature_K']
    Te=np.array([315,300,270,200,190.]);Tw=np.full(5,335.)
    rho=np.array([1.,.9,.8,.6,.5]);u=np.array([60.,100.,200.,300.,400.]);mu=np.full(5,2e-5);k=np.full(5,.03)
    # 独立由已知目标参数反推合成q，核验换算回环；没有使用任何实验q。
    target=np.array([.65,.6,.55,.4,.35])
    Tr=coefficient*T0+(1-coefficient)*Te
    q=-target*k*(Tw-Tr)*np.sqrt(rho*u/(mu*s))
    edge=dict(temperatureK=Te,densityKgM3=rho,speedMS=u,viscosityPas=mu,conductivityWmK=k)
    compared=compare_local(s,q,Tw,edge,'synthetic-local-edge-independent')
    error=float(max(abs(compared['parameter']-target)));assert error<1e-14
    assert not compared['measurementOperatorQualified'] and not compared['physicalQualification']
    rejected=[]
    for name,call in [('missing-electrode-segments',lambda:require_experimental_qualification(None,True,True,True)),
        ('missing-convergence',lambda:require_experimental_qualification(np.tile([.01,.02],(20,1)),True,False,True)),
        ('freestream-as-edge',lambda:compare_local(s,q,Tw,edge,'freestream-substituted-for-local-edge')),
        ('recovery-extrapolation',lambda:published_recovery(1.3*D))]:
        try:call()
        except AssertionError:rejected.append(name)
        else:raise AssertionError('未拒绝：'+name)
    base=ROOT/'outputs/v169_validation'
    correct=read(base/'small3d03/case_manifest.json')
    inlet=inlet_state_check((base/'small3d03/case.cfg').read_text(),contract['freeStream'],correct['expectedVelocity'])
    try:inlet_state_check((base/'small3d02/case.cfg').read_text(),contract['freeStream'],correct['expectedVelocity'])
    except AssertionError:rejected.append('native-small3d02-wrong-Mach')
    else:raise AssertionError('原错误M1输入未被拒绝')
    # 旧小例的可变源码已跨版本变化，完整来源审核应拒绝。此处只比较新旧公共前端：
    # 在相同末端替身处停止，不把这项回归称为旧初场/完整来源重新通过。
    old=[]
    baseline_commit='ae471ae0b9220ace54c67777df2a4a73bf46a790'
    baseline_bytes=subprocess.check_output(['git','show',baseline_commit+':scripts/su2IdealAirInterface.py'],cwd=str(ROOT))
    (out/'baseline_interface.py').write_bytes(baseline_bytes)
    baseline=types.ModuleType('baseline_ideal_air_interface')
    baseline.__file__=str(ROOT/'scripts/su2IdealAirInterface.py')
    exec(compile(baseline_bytes,baseline.__file__,'exec'),baseline.__dict__)
    import su2IdealAirInterface
    for relative in ['outputs/v159_validation/small01','outputs/v156_validation/cases01/restart']:
        folder=ROOT/relative
        m=read(folder/'case_manifest.json');flavor=m['requiredRuntimeFlavor']
        runtime=read(ROOT/'tools/runtime/su2'/('8.6.0-'+flavor)/'runtime.json')
        text=(folder/'case.cfg').read_text()
        target='prepareSU2AxisymmetricAir.validate_source' if 'axisymmetric' in flavor else 'prepareSU2IdealAirRestart.validate_source'
        with patch(target,return_value={}) as endpoint:
            a=baseline.validate(folder,m,text,runtime)
            b=su2IdealAirInterface.validate(folder,m,text,runtime)
            assert endpoint.call_count==2 and a==b
            invalid=text.replace('CONV_RESIDUAL_MINVAL= -8','CONV_RESIDUAL_MINVAL= -7')
            assert invalid!=text
            for validator in (baseline.validate,su2IdealAirInterface.validate):
                try:validator(folder,m,invalid,runtime)
                except AssertionError:pass
                else:raise AssertionError('旧残差门槛被放宽')
            assert endpoint.call_count==2
        old.append(dict(case=relative,runtimeFlavor=flavor,frontEndDispatchAndRejectionMatchBaseline=True,
            sourceValidationReexecuted=False,fullHistoricalCaseRequalified=False))
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    save(out/'report.json',dict(passed=True,syntheticParameterMaximumError=error,inlet=inlet,
        rejected=rejected,oldContracts=old,baselineCommit=baseline_commit,inputSHA256=locks,inputsStable=stable,
        experimentalQualification=False,newQualifiedArea_m2=0))
    print(dict(passed=True,rejected=len(rejected),oldContracts=len(old)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True)
    args=p.parse_args();main(Path(args.output).resolve())
