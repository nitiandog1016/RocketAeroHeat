# -*- coding: utf-8 -*-
"""由已保存真实方向还原候选五量，定位密度/内能/温域护栏；不重新求解NS。"""
from pathlib import Path
import hashlib, json, os, sys
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1];VAL=ROOT/'outputs/v83_validation'
case=sys.argv[1];assert case.replace('_','').isalnum()
target=VAL/(case+'_state_guards.json');assert not target.exists()
evidence=VAL/case/'wall_velocity_krylov_evidence.mat'
domain_file=ROOT/'outputs/v80_validation/domain02/shell_domain_evidence.mat'
def sha(path):
    # 大MAT分块读取，避免在已加载完整状态时额外申请一个整文件缓冲区。
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()
bindings={str(p.relative_to(ROOT)):sha(p) for p in [evidence,domain_file]}
s=loadmat(evidence,simplify_cells=True,variable_names=['result','gas'])
d=loadmat(domain_file,simplify_cells=True,variable_names=['domain'])['domain']
r=s['result'];gas=s['gas'];cv=gas['R']/(gas['gamma']-1);scale=r['residualScale']
nf=int(d['quality']['sourceFaces']);nl=int(d['quality']['layerCount'])
assert len(d['cells']['volume'])==nf*nl and len(d['walls']['ownerCell'])==nf
np.testing.assert_array_equal(d['walls']['ownerCell'],np.arange(1,nf+1))
rows=lambda x:[x] if isinstance(x,dict) else list(x)
traces=rows(r['trace']);trials=rows(r['lineLog']);reports=[];unavailable=[]
def location(index):
    return dict(cellID=int(index+1),layer=int(index//nf+1),sourceFacet=int(index%nf+1),
                centroidM=d['cells']['centroid'][index].tolist())
for row in trials:
    trace=next((t for t in traces if t['iteration']==row['iteration'] and t['pseudoDt']==row['pseudoDt']),None)
    if trace is None:
        unavailable.append(dict(iteration=int(row['iteration']),retry=int(row['retry']),alpha=float(row['alpha']),
                                reason='此重试方向未在接受trace内保存，不能用另一方向重建'))
        continue
    q=trace['initialQ']+row['alpha']*trace['linearUpdateU']*scale
    assert np.isfinite(q).all()
    rho=q[:,0];velocity=q[:,1:4]/rho[:,None]
    internal=q[:,4]-.5*rho*np.sum(velocity*velocity,axis=1);temperature=internal/(rho*cv)
    density_ok=bool(np.all(rho>0));internal_ok=bool(np.all(internal>0))
    temperature_ok=bool(np.all((temperature>=gas['temperatureRange'][0])&(temperature<=gas['temperatureRange'][1])))
    error=row['errorIdentifier'] if isinstance(row['errorIdentifier'],str) else ''
    if error=='CutNS:State':assert not density_ok
    elif error=='CutNS:InternalEnergy':assert density_ok and not internal_ok
    elif error=='CutNS:Temperature':assert density_ok and internal_ok and not temperature_ok
    else:assert density_ok and internal_ok and temperature_ok
    minimum_rho=int(np.argmin(rho));minimum_internal=int(np.argmin(internal))
    minimum_temperature=int(np.argmin(temperature));maximum_temperature=int(np.argmax(temperature))
    reports.append(dict(iteration=int(row['iteration']),retry=int(row['retry']),alpha=float(row['alpha']),accepted=bool(row['accepted']),
                        originalErrorIdentifier=error,positiveDensity=density_ok,positiveInternalEnergy=internal_ok,
                        cellTemperatureInOriginalRange=temperature_ok,minimumDensityKgM3=float(rho[minimum_rho]),
                        minimumDensityLocation=location(minimum_rho),minimumInternalEnergyDensityJm3=float(internal[minimum_internal]),
                        minimumInternalEnergyLocation=location(minimum_internal),algebraicTemperatureRangeK=[float(temperature.min()),float(temperature.max())],
                        temperatureValuesArePhysicalCandidate=density_ok and internal_ok and temperature_ok,
                        minimumTemperatureLocation=location(minimum_temperature),maximumTemperatureLocation=location(maximum_temperature),
                        interfaceOrBoundaryRecomputed=False))
assert all(sha(ROOT/p)==h for p,h in bindings.items())
report=dict(status='reconstructed-real-saved-directions-and-cell-guard-locations',case=case,processID=os.getpid(),
            trials=reports,unavailableDirections=unavailable,fullNSRecomputed=False,interfaceStatesRecomputed=False,
            candidatePromotedToAcceptedSolution=False,physicalSolution=False,newProductionSolvedAreaM2=0,
            inputSHA256=bindings,scriptSHA256=sha(__file__))
target.write_bytes((json.dumps(report,ensure_ascii=False,indent=2)+'\n').encode('utf8'))
print(json.dumps(dict(status=report['status'],reconstructedTrials=len(reports),unavailableDirections=len(unavailable))))
