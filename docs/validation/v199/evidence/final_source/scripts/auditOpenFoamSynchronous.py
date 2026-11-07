"""四步原生控制的独立增量/面通量/全域收支审核，不重跑CFD。"""
from pathlib import Path
import json,time,shutil
import numpy as np
from prepareOpenFoamControl import ROOT,save,sha
from openFoamConservation import geometry,flux_and_residual,conservative,reference_scales,native_mass_flux
from openFoamFields import read_state

BASE=ROOT/'outputs/v199_validation'

def main():
    start=time.monotonic();case=BASE/'small01';out=BASE/'small_audit01';out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/'executed_entry.py');shutil.copyfile(ROOT/'scripts/openFoamConservation.py',out/'executed_oracle.py')
    run=json.loads((BASE/'small_run01/report.json').read_text());assert run['passed']
    for path,h in run['saved'].items():
        for name,digest in h['fieldsSHA256'].items():assert sha(case/path/name)==digest
    g=geometry(case);vol=g['volume'];ni=len(g['neighbor']);scale,speed=reference_scales();rows=[]
    for before,after in [('1e-07','2e-07'),('2e-07','3e-07'),('3e-07','4e-07')]:
        f,flux,residual=flux_and_residual(case/before,g);following=read_state(case/after,len(vol))
        Q=conservative(f);nextQ=conservative(following);dt=float(after)-float(before)
        delta=nextQ-Q;defect=delta+dt*residual
        errors=np.max(abs(defect)/scale,axis=0)
        nativephi=native_mass_flux(case/after,g);phi_error=float(np.max(abs(nativephi-flux[:,0]))/(scale[1]*g['area'].max()))
        total_change=np.sum(delta*vol[:,None],axis=0);boundary_flux=flux[ni:].sum(axis=0)
        balance=abs(total_change+dt*boundary_flux)/(scale*vol.sum())
        assert np.max(errors)<2e-10 and phi_error<2e-10 and max(balance)<2e-10
        rows.append(dict(before=before,after=after,incrementScaledMax=errors.tolist(),nativeMassFluxScaledMax=phi_error,
            integratedChange=total_change.tolist(),boundaryTransport=boundary_flux.tolist(),globalBalanceScaled=balance.tolist(),
            residualVolumeL1Scaled=(np.sum(abs(residual)*vol[:,None],axis=0)/vol.sum()/scale/speed).tolist()))
    result=dict(passed=True,pairs=rows,scaleDefinition='rho_inf, rho_inf*U_inf for each momentum, rho_inf*(abs(es_inf)+U_inf^2); L=1 m',
        scalarScales=scale.tolist(),velocityScaleMps=speed,source='actual stock solver consecutive saved fields and actual phi',
        physicalQualification=False,numericallyConverged=False,elapsedSeconds=time.monotonic()-start)
    save(out/'report.json',result);print(json.dumps(result))

if __name__=='__main__':main()
