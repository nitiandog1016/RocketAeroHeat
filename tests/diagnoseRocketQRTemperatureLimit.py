"""独立复算真实拒绝界面的温度和局部可行半径，不改温域或构造新物理解。"""
from pathlib import Path
import json, os
import numpy as np
from scipy.io import loadmat
from scipy.optimize import brentq
from riemannReservoirReference import ROOT, sha, exact, PUBLIC_BINDINGS
OUT=ROOT/'outputs/v85_validation'


def primitive(q,gas):
    rho=q[0];v=q[1:4]/rho
    return np.r_[rho,v,(q[4]-.5*rho*float(v@v))/(rho*gas['R']/(gas['gamma']-1))]


def main():
    rejected=OUT/'rejected01/rejected_interface_evidence.mat'
    state=OUT/'qr_real01/krylov01/wall_velocity_krylov_evidence.mat'
    bindings={str(p.relative_to(ROOT)):sha(p) for p in (rejected,state)}
    s=loadmat(state,simplify_cells=True,variable_names=['result','gas','boundary']);r=s['result'];gas=s['gas']
    d=loadmat(rejected,simplify_cells=True);report=d['report'];normal=d['normal'];i=int(report['ownerCV'])-1;j=int(report['neighbourCV'])-1
    trace=r['trace'] if isinstance(r['trace'],dict) else r['trace'][0]
    initial=r['initialQ'][[i,j]];direction=trace['linearUpdateU'][[i,j]]*r['residualScale'];accepted_alpha=float(trace['alpha']);rejected_alpha=float(report['alpha'])
    np.testing.assert_array_equal(initial+rejected_alpha*direction,np.row_stack([d['ql'],d['qr']]))
    rows=[]
    def sample(alpha,save=False):
        q=initial+alpha*direction;left=primitive(q[0],gas);right=primitive(q[1],gas)
        sample,flux,region=exact(left,right,normal,gas)
        if save:
            cp=gas['gamma']*gas['R']/(gas['gamma']-1)
            rows.append(dict(alpha=alpha,cellTemperatureK=[float(left[4]),float(right[4])],sampleTemperatureK=float(sample[4]),
                             cellStagnationTemperatureK=[float(x[4]+x[1:4]@x[1:4]/(2*cp)) for x in (left,right)],
                             sampleRegion=region,primitiveLeft=left.tolist(),primitiveRight=right.tolist()))
        return float(sample[4])
    for alpha in (0.,accepted_alpha,rejected_alpha):sample(alpha,True)
    assert rows[0]['sampleTemperatureK']<=2000 and rows[1]['sampleTemperatureK']<=2000 and rows[2]['sampleTemperatureK']>2000
    critical=brentq(lambda a:sample(a)-2000,accepted_alpha,rejected_alpha,xtol=1e-15,rtol=1e-14)
    expected=float(__import__('re').search(r'T=([0-9.]+)',report['message']).group(1))
    difference=abs(rows[2]['sampleTemperatureK']-expected)/expected;assert difference<1e-12
    out=dict(status='passed-independent-public-Riemann-reference-real-rejected-interface-and-local-upper-temperature-root',
             globalInteriorFaceID=int(report['globalInteriorFaceID']),ownerCV=i+1,neighbourCV=j+1,
             sourceFacet=int(report['ownerSourceFacet']),centroidM=report['centroidM'].tolist(),
             samples=rows,localFace2000KCrossingAlpha=critical,acceptedFractionOfLocalFaceRadius=accepted_alpha/critical,
             globalAdmissibilityRadiusClaimed=False,temperatureRootRelativeDifference=difference,
             referenceKind='public-pressure-curves-and-independent-Brent-sampling-arithmetic-not-experimental-validation',
             completeNSRecomputed=False,newActualNSEvaluations=0,newFlowAcceptances=0,newProductionSolvedAreaM2=0,
             physicalSolution=False,inputSHA256={**bindings,**PUBLIC_BINDINGS},scriptSHA256=sha(__file__),processID=os.getpid())
    assert all(sha(ROOT/path)==h for path,h in bindings.items())
    target=OUT/'qr_temperature_limit.json';assert not target.exists();target.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in out.items() if k!='inputSHA256'}))


if __name__=='__main__':main()
