"""真实全壁初态面波速和已测成本核验；仅只读预算，不运行CFD。"""
from pathlib import Path
import json,re,time,shutil
import numpy as np
from prepareOpenFoamControl import ROOT,sha,save
from openFoamConservation import geometry,boundary_values
from openFoamFields import internal_field
from openFoamAirThermo import evaluate


def main():
    start=time.monotonic();base=ROOT/'outputs/v201_validation';out=base/'budget01';out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');source=ROOT/'outputs/v200_validation/rocket03'
    accepted=json.loads((ROOT/'outputs/v200_validation/rocket_check03/report.json').read_text());assert accepted['passed']
    prep=json.loads((ROOT/'outputs/v200_validation/rocket_relocation_audit01/report.json').read_text());assert prep['passed']
    for p,h in prep['inputSHA256'].items():assert sha(p)==h
    # 模块一：准确原生源码的amaxSf/correctCoNum；不是有量纲残差或最短边估算。
    src=ROOT/'outputs/v199_validation/source_review01/source/applications/modules/shockFluid'
    for name in ['fluxPredictor.C','shockFluid.C']:
        shutil.copyfile(src/name,out/name)
    text=(src/'shockFluid.C').read_text();assert '0.5*gMax(sumAmaxSf/mesh.V().primitiveField())*runTime.deltaTValue()' in text
    g=geometry(source);n=len(g['cells']);ni=len(g['neighbor']);owner=g['owner'];neighbor=g['neighbor']
    f={name:internal_field(source/'0'/name,n,name=='U') for name in ['p','T','U']}
    b={name:boundary_values(source/'0',name,f[name],g) for name in f}
    U=np.r_[f['U'][owner[:ni]],b['U']];T=np.r_[f['T'][owner[:ni]],b['T']];p=np.r_[f['p'][owner[:ni]],b['p']]
    ur=np.r_[f['U'][neighbor],b['U']];tr=np.r_[f['T'][neighbor],b['T']];pr=np.r_[f['p'][neighbor],b['p']]
    phi=np.einsum('ij,ij->i',U,g['sf']);phir=np.einsum('ij,ij->i',ur,g['sf'])
    c=evaluate(p,T)['soundSpeedMps']*g['area'];cr=evaluate(pr,tr)['soundSpeedMps']*g['area']
    ap=np.maximum.reduce([phi+c,phir+cr,np.zeros(len(phi))]);am=np.minimum.reduce([phi-c,phir-cr,np.zeros(len(phi))])
    weight=ap/(ap-am);asf=am*weight;alpha=weight*phi-asf;beta=(1-weight)*phir+asf
    amax=np.maximum(abs(alpha),abs(beta));summed=np.bincount(owner,weights=amax,minlength=n)+np.bincount(neighbor,weights=amax[:ni],minlength=n)
    rate=.5*summed/g['volume'];index=int(np.argmax(rate));target_dt=.2/float(rate.max())
    speed=float(np.linalg.norm(f['U'][0]));mapping=json.loads((source/'wall_mapping.json').read_text());centers=np.array(mapping['wallCentersM'])
    surface=json.loads((ROOT/'outputs/v194_validation/rocket_design03/surface.json').read_text());vertices=np.array(surface['vertices']);length=float(np.ptp(vertices[:,0]));transit=length/speed
    control=json.loads((ROOT/'outputs/v200_validation/restart_run01/report.json').read_text());steps=3894
    step_seconds=control['process']['elapsedSeconds']/steps*n/14892
    report=dict(passed=True,cellCount=n,initialMaxCoAtDt1e7=float(rate.max()*1e-7),targetMaxCo=.2,
        initialTargetDeltaTSeconds=target_dt,limitingCellZeroBased=index,limitingCenterM=g['centers'][index].tolist(),
        minimumCellVolumeM3=float(g['volume'].min()),lengthReferenceM=length,speedReferenceMps=speed,bodyTransitSeconds=transit,
        bodyTransitStepsAtInitialDt=transit/target_dt,linearCellStepCostSeconds=step_seconds,
        extrapolatedSecondsOneBodyTransit=step_seconds*transit/target_dt,
        timingScope='linear extrapolation from actual 14892-cell control; not measured rocket cost or convergence guarantee',
        originalWallFaces=19752,originalAreaM2=mapping['originalAreaM2'],
        bytesInitialCase=sum(p.stat().st_size for p in source.rglob('*') if p.is_file()),
        sourceSHA256={str(p):sha(p) for p in [source/'exact_mesh.npz',source/'0/p',source/'0/T',source/'0/U',src/'fluxPredictor.C',src/'shockFluid.C',ROOT/'outputs/v200_validation/rocket_check03/report.json']},
        CFDRuns=0,physicalQualification=False,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',report);print(json.dumps(report))


if __name__=='__main__':main()
