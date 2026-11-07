"""同基态原生导数的分变量独立比较；读CSV，不以程序exit0授数值或物理资格。"""
from pathlib import Path
import argparse,json,time
import numpy as np


def check(base,out):
    start=time.monotonic();base=Path(base);out=Path(out);out.mkdir(exist_ok=False,parents=True)
    records={};arrays={}
    def relative(a,b,ref):
        return [float(np.linalg.norm((a-b)[ref['variable']==v])/max(1.,np.linalg.norm(b[ref['variable']==v]))) for v in range(7)]
    for name in ['real-state0','complex-state0','real-state2','complex-state2']:
        a=np.genfromtxt(base/name/'operator.csv',delimiter=',',names=True);arrays[name]=a
        assert len(a)==28 and all(np.isfinite(a[k]).all() for k in a.dtype.names)
        d={}
        for key in ['J_native','J_sparse','J_mixed']:
            d[key+'_vs_central']=relative(a[key],a['J_central'],a)
        d['sparse_vs_native']=relative(a['J_sparse'],a['J_native'],a)
        for key in ['R_after_matrix','R_after_jv','R_restored']:
            d[key+'_vs_baseline']=relative(a[key],a['R'],a)
        for key in ['U_after_matrix','U_after_jv','U_after_R0','U_final']:
            d[key+'_vs_original']=relative(a[key],a['U'],a)
        d['high_R_after_R0']=relative(a['R_high_after_R0'],a['R_high_before'],a)
        d['all_derivatives_within_1e_4']=all(max(d[key])<=1e-4 for key in ['J_native_vs_central','J_sparse_vs_central','J_mixed_vs_central'])
        d['restoration_within_1e_12']=all(max(v)<=1e-12 for k,v in d.items() if k.endswith('_vs_baseline') or k.endswith('_vs_original'))
        records[name]=d
    cross={}
    for snapshot in [0,2]:
        a,b=arrays['real-state'+str(snapshot)],arrays['complex-state'+str(snapshot)]
        cross[str(snapshot)]={k:relative(a[k],b[k],a) for k in ['U','z','R','J_native','J_sparse','J_central','J_mixed']}
    r=dict(passed=all(x['all_derivatives_within_1e_4'] and x['restoration_within_1e_12'] for x in records.values()),
           variableOrder=['mass','x-momentum','y-momentum','z-momentum','total-energy','rho-k','rho-omega'],
           records=records,crossType=cross,nativeDevelopmentSteps=0,physicalQualification=False,newQualifiedAreaM2=0,
           note='Norms are per conservative-variable family across four cells; denominator max(1,norm(reference)). Central differences are independent numerical comparisons, not exact truth at nonsmooth states.',
           elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(dict(passed=r['passed'],records={n:{k:(max(v) if isinstance(v,list) else v) for k,v in d.items()} for n,d in records.items()},crossType=cross)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--out',required=True);a=p.parse_args();check(a.base,a.out)
