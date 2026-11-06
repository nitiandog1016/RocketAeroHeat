"""模块用途：独立解析/积分真值核验重叠测量段，拒绝外插、错号和无来源。"""
import argparse
from pathlib import Path
import sys
import time
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from nacaMeasurementSegments import segments,integrate_profile,reduce_segment,transport_reduction_ratio
from naca3344Model import D,R,SHOULDER,END,heat_parameter
from su2ThermalRuntime import sha
from validateIdealAir import save


def verify(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    paths=[Path(__file__).resolve(),ROOT/'scripts/nacaMeasurementSegments.py',ROOT/'scripts/naca3344Model.py',
           ROOT/'scripts/naca3344Normalization.py',ROOT/'tests/data/naca_tn3344/segments_v170/summary.json']
    locks={str(p):sha(p) for p in paths};save(out/'input_manifest.json',locks)
    (out/'executed_verify.py').write_bytes(Path(__file__).read_bytes())
    (out/'executed_operator.py').write_bytes((ROOT/'scripts/nacaMeasurementSegments.py').read_bytes())
    source=dict(fieldSource='analytic-independent-fixture',circumferentialReduction='analytic-axisymmetric-field')
    s=np.unique(np.r_[np.linspace(0,END,89),SHOULDER]);q=-2000-30000*s
    answer=integrate_profile(s,q,source);errors=[]
    for row in answer['segments']:
        a,b=row['boundsM']
        weight=lambda x:2*np.pi*R*(np.sin(x/R) if x<SHOULDER else 1.)
        oracle=quad(lambda x:(-2000-30000*x)*weight(x),a,b,points=[SHOULDER] if a<SHOULDER<b else [],epsabs=1e-11)[0]
        errors.append(abs(row['powerIntoWallW']/oracle-1))
    assert max(errors)<2e-12 and answer['sumOverlappingAreasM2']>1.9*answer['nominalMeasuredUnionAreaM2']
    constant=integrate_profile(s,np.full(len(s),-3210.),source)
    assert max(abs(row['meanQIntoWallWm2']+3210) for row in constant['segments'])<1e-8
    mixed=integrate_profile(s,50000*(s-.03),source)
    assert min(r['meanQIntoWallWm2'] for r in mixed['segments'])<0<max(r['meanQIntoWallWm2'] for r in mixed['segments'])
    edge=dict(temperatureK=240.,densityKgM3=.8,speedMS=350.,viscosityPas=1.5e-5,conductivityWmK=.021)
    rows=[reduce_segment(row,'oppositeTemperatureStation',335.,edge,'synthetic-local-edge') for row in answer['segments']]
    assert all(not row['experimentalQualification'] for row in rows)
    assert sum(row['strictLocalAbscissaIdentityClosed'] for row in rows)==14
    second=dict(edge,densityKgM3=.9,speedMS=345.,viscosityPas=1.6e-5,conductivityWmK=.023)
    a=heat_parameter(.05,-3000.,335.,edge,300.,dict(edge='a',recovery='a'))['parameter']
    b=heat_parameter(.05,-3000.,335.,second,302.,dict(edge='b',recovery='b'))['parameter']
    ratio=transport_reduction_ratio(edge,second,335.,300.,302.);ratio_error=abs(ratio/(b/a)-1)
    assert ratio_error<1e-14
    rejected=[]
    for name,call in [('missing-source',lambda:integrate_profile(s,q,{})),
                      ('missing-start',lambda:integrate_profile(s[1:],q[1:],source)),
                      ('unordered',lambda:integrate_profile(s[::-1],q[::-1],source)),
                      ('nonfinite',lambda:integrate_profile(s,np.full(len(s),np.nan),source)),
                      ('wrong-heat-sign',lambda:reduce_segment(mixed['segments'][-1],'plot',335.,edge,'synthetic-local-edge')),
                      ('freestream-substitution',lambda:reduce_segment(answer['segments'][0],'plot',335.,edge,'freestream-substituted-for-local-edge'))]:
        try:call()
        except AssertionError:rejected.append(name)
        else:raise AssertionError('应拒绝 '+name)
    stable={p:sha(Path(p))==h for p,h in locks.items()};assert all(stable.values())
    report=dict(passed=True,independentLinearProfileIntegralMaxRelativeError=max(errors),
        transportRatioIndependentRelativeError=ratio_error,nominalSegments=20,strictCenterIdentityRows=14,
        unionAreaM2=answer['nominalMeasuredUnionAreaM2'],sumOverlappingAreasM2=answer['sumOverlappingAreasM2'],
        rejected=rejected,inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started,
        experimentalQualified=False,actualCFDSurfaceMapped=False)
    save(out/'report.json',report);save(out/'synthetic_segments.json',answer);save(out/'synthetic_reduction.json',rows)
    print({k:v for k,v in report.items() if k not in ('inputSHA256','inputsStable')})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True)
    verify(Path(parser.parse_args().output).resolve())
