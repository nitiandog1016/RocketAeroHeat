"""新物性独立积分/导数、公开表、源接缝与旧模型重叠域检查。"""
from pathlib import Path
import json,sys,time,shutil
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from nistIntegratedAir import AIR,evaluate,mixture_reference,species_properties,segment_properties,contract,R,RU,DATA,M,Y,SMIX
from openFoamAirThermo import evaluate as old_evaluate


def main():
    start=time.monotonic();out=ROOT/'outputs/v202_validation/model_check01';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_test.py');shutil.copyfile(ROOT/'scripts/nistIntegratedAir.py',out/'executed_model.py')
    temperatures=np.unique(np.r_[np.linspace(100,2500,24001),[499.999999,500,699.999999,700,1999.999999,2000,298.15]])
    g=evaluate(5000,temperatures);cp,h,s=mixture_reference(temperatures);_,rawh,raws=mixture_reference(temperatures,False)
    error=dict(cpRelative=float(abs(g['cpJkgK']/cp-1).max()),hIntegralJkg=float(abs(g['enthalpyJkg']-h).max()),
        sIntegralJkgK=float(abs(g['entropyJkgK']+R*np.log(5000/100000)-s).max()),
        rawNistHJkg=float(abs(h-rawh).max()),rawNistSJkgK=float(abs(s-raws).max()))
    seams=[]
    for name in ['N2','O2']:
        segments=DATA[name]['segments']
        for a,b in zip(segments,segments[1:]):
            t=b['rangeK'][0];left=segment_properties(t,a['coefficients']);right=segment_properties(t,b['coefficients'])
            seams.append(dict(species=name,temperatureK=t,CpJumpJmolK=float(right[0]-left[0]),HJumpJmol=float(right[1]-left[1]),SJumpJmolK=float(right[2]-left[2])))
    # 模块一：独立数值积分从原Shomate Cp出发，显式在分段处拆分。
    integral_errors=[]
    for t in [100,170,250,600,1000,1800,2500]:
        lo,hi=sorted([298.15,t]);sign=1 if t>=298.15 else -1;breaks=[b for b in [500,700,2000] if lo<b<hi]
        hv=sign*quad(lambda x:float(mixture_reference(x)[0]),lo,hi,points=breaks,epsabs=1e-7,epsrel=1e-12)[0]
        sv=sign*quad(lambda x:float(mixture_reference(x)[0])/x,lo,hi,points=breaks,epsabs=1e-9,epsrel=1e-12)[0]+AIR.sf
        _,he,se=mixture_reference(t);integral_errors.append([float(abs(he-hv)),float(abs(se-sv))])
    ts=np.array([120.17,226.50908361133,450.23,650.19,1000.27,1800.13,2400.23]);eps=.001
    plus=evaluate(5000,ts+eps);minus=evaluate(5000,ts-eps);middle=evaluate(5000,ts)
    derivative_h=float(abs((plus['enthalpyJkg']-minus['enthalpyJkg'])/(2*eps)/middle['cpJkgK']-1).max())
    derivative_e=float(abs((plus['internalEnergyJkg']-minus['internalEnergyJkg'])/(2*eps)/middle['cvJkgK']-1).max())
    # 模块二：不同出版表/拟合的差异保留，不把表值的三位舍入当机器精度。
    public={'N2':[[100,29.104,159.811,-5.768],[200,29.107,179.985,-2.857],[298.15,29.124,191.609,0],[1000,32.697,228.170,21.463],[2500,36.616,260.176,74.296]],
        'O2':[[100,29.106,173.307,-5.779],[200,29.126,193.485,-2.868],[298.15,29.376,205.147,0],[1000,34.870,243.578,22.703],[2500,38.856,277.290,78.328]],
        'Ar':DATA['Ar']['janafRows_T_Cp_S_HkJmol']}
    public_errors={}
    for name,rows in public.items():
        a=np.array(rows);cc,hh,ss=species_properties(name,a[:,0]);public_errors[name]=dict(cpRelative=float(abs(cc/a[:,1]-1).max()),
            hJmol=float(abs(hh-a[:,3]*1000).max()),sJmolK=float(abs(ss-a[:,2]).max()))
    # 模块三：重叠域先对齐参考量，再报告物性变化，不把旧气体数据当新真值。
    to=np.linspace(200,2500,2301);old=old_evaluate(5000,to);new=evaluate(5000,to);o_ref=old_evaluate(5000,298.15);n_ref=evaluate(5000,298.15)
    overlap=dict(cpMaxRelative=float(abs(new['cpJkgK']/old['cpJkgK']-1).max()),
        alignedHMaxJkg=float(abs((new['enthalpyJkg']-n_ref['enthalpyJkg'])-(old['enthalpyJkg']-o_ref['enthalpyJkg'])).max()),
        alignedSMaxJkgK=float(abs((new['entropyJkgK']-n_ref['entropyJkgK'])-(old['entropyJkgK']-o_ref['entropyJkgK'])).max()),
        nativeSpeedAtRocketM5Mps=5*float(evaluate(1197.0279007383137,226.50908361133003)['soundSpeedMps']))
    rejects=[]
    for p,t in [(0,300),(-1,300),(100001,300),(5000,99.99),(5000,2500.01),(5000,float('nan')),(float('inf'),300)]:
        try:evaluate(p,t)
        except ValueError:rejects.append([str(p),str(t)]);continue
        raise AssertionError('未拒绝超域输入')
    result=dict(passed=False,contract=contract(),representationErrors=error,sourceSeams=seams,independentIntegralErrors=integral_errors,
        derivativeHRelative=derivative_h,derivativeERelative=derivative_e,publicTableErrors=public_errors,oldNASA7Overlap=overlap,
        rejects=rejects,minimumCvJkgK=float(g['cvJkgK'].min()),elapsedSeconds=time.monotonic()-start,physicalQualification=False)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    assert error['cpRelative']<=2e-4 and error['hIntegralJkg']<=5 and error['sIntegralJkgK']<=.02
    assert error['rawNistHJkg']<=50 and error['rawNistSJkgK']<=.1
    assert np.max(integral_errors)<1e-6 and max(derivative_h,derivative_e)<2e-8 and g['cvJkgK'].min()>0
    assert all(v['cpRelative']<.001 and v['hJmol']<5 and v['sJmolK']<.05 for v in public_errors.values())
    result['passed']=True;(out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    (out/'mixture.dict').write_text(AIR.dictionary(),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['passed','representationErrors','publicTableErrors','oldNASA7Overlap','elapsedSeconds']}))


if __name__=='__main__':main()
