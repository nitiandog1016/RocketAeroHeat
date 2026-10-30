"""模块用途：独立检查一阶显式薄层压力的轴线能量极限，不把此极限当新解。"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq
from auditPressureDrivenSurfaceFlow import rows
from auditShockSurfaceInlet import thermal,jump
from auditThermalPressureSurfaceFlow import sha


def audit(evidence,target):
    assert not target.exists();saved=loadmat(evidence,simplify_cells=True);fs=saved['fs'];cases=[]
    up={'pressure':fs['pressure'],'temperature':fs['temperature'],'density':fs['density'],'velocityVector':fs['flowVelocityVector']}
    frame=rows(saved['cases'])[0]['column']['shape']['frame'];U=np.linalg.norm(up['velocityVector']);A=frame@np.vstack([np.zeros(2),np.eye(2)])
    for case in rows(saved['cases']):
        gas=case['column']['shock']['gas'];jumpState,velocity,n=jump(up,gas,A,1);ps,Ts,rhos,_,S=jumpState;vs=velocity@n
        H0=thermal(up['temperature'],up['pressure'],gas)[1]+U*U/2
        T0=brentq(lambda T:thermal(T,ps,gas)[1]-H0,*gas['temperatureRangeK'],xtol=2e-11)
        p0=ps*np.exp((thermal(T0,ps,gas)[2]-S)/gas['R'])
        profile=[]
        for eta in [1.,.75,.5,.25,0.]:
            # q→0时：C→1，P1→0，P2→rho_inf*U_inf*vs/2。
            pressure=ps+.5*up['density']*U*vs*(1-eta*eta);v=vs*eta
            T=brentq(lambda T:thermal(T,pressure,gas)[2]-S,*gas['temperatureRangeK'],xtol=2e-11)
            squared=2*(H0-thermal(T,pressure,gas)[1])-v*v
            if eta!=1:assert squared>0
            profile.append({'eta':eta,'explicitPressurePa':float(pressure),'temperatureK':float(T),
                            'remainingTangentialSpeedSquaredM2S2':float(squared),'normalSpeedMS':float(v)})
        final=profile[-1];assert final['explicitPressurePa']<p0
        cases.append({'gasMode':gas['thermalPropertyMode'],'normalShockPressurePa':float(ps),
                      'normalShockTemperatureK':float(Ts),'normalShockNormalSpeedMS':float(vs),
                      'isentropicStagnationPressurePa':float(p0),'isentropicStagnationTemperatureK':float(T0),
                      'explicitAxisEndPressurePa':final['explicitPressurePa'],
                      'pressureUnderpredictionRelativeToExactAxisStagnation':float(1-final['explicitPressurePa']/p0),
                      'spuriousAxisTangentialSpeedMS':float(np.sqrt(final['remainingTangentialSpeedSquaredM2S2'])),
                      'profile':profile})
    report={'status':'confirmed-first-order-pressure-axis-energy-defect-before-global-matching',
            'evidenceSHA256':sha(evidence),'programSHA256':sha(Path(__file__)),'cases':cases,
            'scope':'显式一阶薄层压力的解析轴线极限诊断；不是实际驻点或新增热流解',
            'stagnationGlobalMatchingSolved':False,'productionQualified':False,'newProductionSolvedAreaM2':0}
    target.write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode('utf8'));print(json.dumps(report))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('evidence',type=Path);parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();audit(args.evidence,args.report)
