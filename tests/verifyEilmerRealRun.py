"""只读核对新实数小例与既有复数小例，保留原生字段名与第一次读取失败。"""
from pathlib import Path
import argparse,hashlib,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from readEilmerNative import read_field,export_loads,audit_numerics
from auditEilmerCompletion import audit_completion


def verify(out):
    start=time.monotonic();out=Path(out);out.mkdir(exist_ok=False,parents=True)
    base=ROOT/'outputs/v181_validation';case=base/'small01';ref=ROOT/'outputs/v180_validation/small01'
    r=json.loads((base/'small01_run/report.json').read_text());assert r['passed']
    assert r['commands'][0]['argv'][0]=='lmr-run' and r['commands'][0]['nativeExitCode']==0
    m=json.loads((case/'case_manifest.json').read_text());audit_numerics(json.loads((case/'lmrsim/nk-config').read_text()),m['condition'])
    h=np.genfromtxt(case/'lmrsim/diagnostics/nk-diagnostics',names=True,deletechars='')
    assert len(h)==8 and all(np.isfinite(h[k]).all() for k in h.dtype.names)
    new=read_field(case,2,4);old=read_field(ref,2,4);errors={}
    for k in ['p','T','rho','vel.x','vel.y','vel.z','tq-tke','tq-omega']:
        errors[k]=float(np.max(np.abs(new[k]-old[k]))/max(1.,np.max(np.abs(old[k]))))
    load=case/'lmrsim/loads/0003/blk-0000-bndry-0-wall.dat'
    wall=export_loads(case,load,out/'wall.csv',base/'small01_run/report.json')
    a=np.genfromtxt(load,names=True,deletechars='');b=np.genfromtxt(ref/'lmrsim/loads/0003/blk-0000-bndry-0-wall.dat',names=True,deletechars='')
    delta={k:float(np.sum(a['area']*np.abs(a[k]-b[k]))/np.sum(a['area']*np.abs(b[k]))) for k in ['p','q_total']}
    wrong=out/'wrong-number-type.txt';wrong.write_text((base/'small01_run/process00/stdout.txt').read_text().replace('Number-type: real','Number-type: complex'))
    try:audit_completion(case,wrong)
    except ValueError as e:rejection=str(e)
    else:raise AssertionError('错二进制没有拒绝')
    assert rejection=='wrong-native-number-type-for-real-profile'
    passed=max(errors.values())<=1e-3 and delta['p']<=1e-3 and delta['q_total']<=1e-2
    result=dict(passed=passed,fieldMaxScaledDifferences=errors,wallAreaL1Differences=delta,wallExport=wall,
                phase=h['phase'].tolist(),lastGlobalRelativeResidual=float(h['global-residual-rel'][-1]),
                maximumLinearResidual=float(np.max(h['linear-solve-residual-rel-achieved'])),linearIterations=h['n-iters'].tolist(),
                peakGroupRSSBytes=r['commands'][0]['peakObservedGroupRSSBytes'],
                binaryMismatchRejected=rejection,physicalQualification=False,newQualifiedAreaM2=0,newNativeRuns=0,
                sourceSHA256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
                              ['tests/verifyEilmerRealRun.py','scripts/readEilmerNative.py','scripts/auditEilmerCompletion.py']},
                firstReadFailure='run_small.py requested tke; accurate metadata uses tq-tke and tq-omega. Native run succeeded; no rerun.',
                elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));assert passed


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);verify(p.parse_args().out)
