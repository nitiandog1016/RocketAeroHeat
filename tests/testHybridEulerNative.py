# -*- coding: utf-8 -*-
"""原生Euler准备/物性/映射检查；独立守恒参照，不给热资格。"""
from pathlib import Path
import argparse
import json
import sys
import shutil
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from readHybridEulerNative import read_wall_frame
from hybridEulerReference import control_reference
from eilmerAirThermo import evaluate, bind_model
from readEilmerNative import sha


def check(case, preparation, out, native=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    root=Path(__file__).resolve().parents[1]
    sources=[Path(__file__)]+[root/'scripts'/n for n in ['readHybridEulerNative.py','eilmerAirThermo.py',
                       'hybridEulerReference.py','readEilmerNative.py','hybridExternalState.py']]
    for p in sources:shutil.copyfile(p,out/p.name)
    hashes={str(p):sha(p) for p in sources}
    result=dict(passed=False,physicalQualification=False)
    try:
        if native is None:
            wall,audit=read_wall_frame(case,0,preparation);reference=control_reference()
            r=reference['residuals']
            assert max(r[k] for k in ['massRelative','normalMomentumRelative','totalEnthalpyRelative'])<1e-11
            assert r['entropyExpansionAbsolute']<1e-8 and r['expansionAngleRad']<1e-11
            assert r['cpIntegralEnthalpyAbsolute']<1e-6
            inp=json.loads((Path(case)/'native-input-state.json').read_text())
            g=evaluate(inp['p'],inp['T'])
            for target,native_key in [('densityKgm3','rho'),('internalEnergyJkg','e'),
                                      ('soundSpeedMps','a'),('cpJkgK','cp')]:
                assert abs(float(g[target])-inp[native_key])<2e-10*max(1,abs(inp[native_key]))
            assert np.allclose(wall['fields']['velocityMps'],[inp['velx'],0,0],rtol=0,atol=1e-12)
            assert wall['wallHeatFluxWm2']==[None]*240 and not wall['engineeringHeatReady']
            for name,obj in [('initial_wall.json',wall),('reference.json',reference),('prepared_audit.json',audit)]:
                (out/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
            result.update(prepared=audit,reference=reference,fullVolumeEOSErrors=wall['fullVolumeEOSErrors'])
        else:
            data=np.genfromtxt(native,delimiter=',',names=True)
            g=evaluate(data['p'],data['T']);errors={}
            for key,col in [('densityKgm3','rho'),('internalEnergyJkg','e'),('enthalpyJkg','h'),
                            ('entropyJkgK','s'),('soundSpeedMps','a'),('cpJkgK','cp'),
                            ('viscosityPas','mu'),('conductivityWmK','k')]:
                # 输运不以1作分母，防止微小数被绝对容差淹没。
                scale=np.maximum(1e-30,abs(data[col]))
                error=float(np.max(abs(g[key]-data[col])/scale));errors[key]=error
                assert error<2e-10,(key,error)
            rejected=[]
            for name,p,t in [('below',5000,199.99),('above',5000,6000.01),('nan',5000,float('nan')),
                             ('pressure',0,300),('infinity',float('inf'),300)]:
                try:evaluate(p,t)
                except ValueError:rejected.append(name)
                else:raise AssertionError('未拒绝: '+name)
            try:bind_model(Path(case)/'rah-air.gas','0'*64)
            except ValueError:rejected.append('gas-source')
            else:raise AssertionError('气体身份未拒绝')
            # 同源EOS吻合不代表热力学微分一致：混合区明确量化导数差。
            derivative=[]
            for t in [300.,700.,900.,1000.,1100.,1300.]:
                dt=.001;cp=float(evaluate(5000,t)['cpJkgK'])
                dh=float((evaluate(5000,t+dt)['enthalpyJkg']-evaluate(5000,t-dt)['enthalpyJkg'])/(2*dt))
                derivative.append(dict(T=t,cp=cp,dh_dT=dh,relativeDifference=(dh-cp)/cp))
            result.update(nativeStates=len(data),errors=errors,rejected=rejected,enthalpyDerivative=derivative,
                          nativeCSV_SHA256=sha(native))
        result['passed']=True
    except Exception as error:result['error']=repr(error)
    result['sourcesStable']=all(sha(p)==digest for p,digest in hashes.items())
    result['sourceSHA256']=hashes;result['elapsedSeconds']=time.monotonic()-start
    (out/'report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(passed=result['passed'],error=result.get('error'),elapsedSeconds=result['elapsedSeconds'])))
    return result['passed'] and result['sourcesStable']


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--preparation',required=True)
    p.add_argument('--out',required=True);p.add_argument('--native')
    a=p.parse_args();raise SystemExit(0 if check(a.case,a.preparation,a.out,a.native) else 1)
