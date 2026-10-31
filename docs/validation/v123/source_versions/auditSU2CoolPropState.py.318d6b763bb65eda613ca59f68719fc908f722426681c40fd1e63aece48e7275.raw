"""用同库独立程序审核真实原生守恒场；不套定比热EOS或强制内能零点。"""
from pathlib import Path
import argparse,csv,hashlib,json,os
import numpy as np
import buildSU2ThermalRuntime as supervisor
from su2ThermalRuntime import COOLPROP_FLAVORS
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def audit(folder,process_runner=None):
    folder=Path(folder).resolve();folder.relative_to(ROOT/'outputs')
    engine=read(folder/'engine_process.json');assert engine['status']=='terminal' and engine['exitCode']==0 and engine['runtimeFlavor'] in COOLPROP_FLAVORS
    assert all(engine['inputsStable'].values()) and not engine['timedOut'] and not engine['cancelled']
    restart=folder/'restart_flow.csv';assert sha(restart)==engine['outputSHA256'][restart.name]
    runtimeFile=ROOT/'tools/runtime/su2'/('8.6.0-'+engine['runtimeFlavor'])/'runtime.json';runtime=read(runtimeFile)
    assert sha(runtimeFile)==engine['inputSHA256'][str(runtimeFile)]
    helper=Path(runtime['propertyExecutable']);assert sha(helper)==next(q['sha256'] for q in runtime['files'] if q['path']=='bin/RAH_CoolProp_State.exe')
    state=np.loadtxt(str(restart),delimiter=',',skiprows=1);nd=(state.shape[1]-3)//2
    assert nd in [2,3] and state.shape[1]==2*nd+3 and np.array_equal(state[:,0],np.arange(len(state))) and np.isfinite(state).all()
    rho=state[:,nd+1];assert min(rho)>0;velocity=state[:,nd+2:2*nd+2]/rho[:,None]
    e=state[:,-1]/rho-.5*np.sum(velocity**2,axis=1)
    out=folder/'thermodynamics';assert not out.exists();out.mkdir();src=out/'input.csv';dst=out/'properties.csv'
    src.write_text('mode,a,b\n'+''.join('DE,%.17g,%.17g\n'%q for q in zip(rho,e)))
    locks={str(p):sha(p) for p in [Path(__file__).resolve(),restart,folder/'engine_process.json',runtimeFile,helper,src]}
    supervisor.OUT=out
    run=supervisor.run if process_runner is None else lambda *a:process_runner(out,*a)
    run('state',[str(helper),str(src),str(dst)],ROOT,dict(os.environ,OMP_NUM_THREADS='1'),120)
    assert all(sha(Path(p))==h for p,h in locks.items())
    data=np.atleast_1d(np.genfromtxt(str(dst),names=True,delimiter=','));assert np.array_equal(data['row'],np.arange(len(state)))
    assert np.max(abs(data['rho_kg_m3']/rho-1))<1e-8 and np.max(abs(data['e_J_kg']-e)/np.maximum(abs(e),100000))<1e-8
    assert min(data['p_Pa'])>0 and min(data['T_K'])>0
    report=dict(thermodynamicModel='CoolProp-7.2.0-HEOS-Air',sourceStateSHA256=sha(restart),propertyOutputSHA256=sha(dst),
       processReceiptSHA256=sha(out/'state_process.json'),sourceSHA256=locks,sourceStable=True,rows=len(state),spatialDimensions=nd,
       minimumSpecificInternalEnergy=float(min(e)),internalEnergyZeroPointAssumed=False,
       temperatureRangeK=[float(min(data['T_K'])),float(max(data['T_K']))],physicsQualified=False,newQualifiedArea=0)
    (out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--folder',type=Path,required=True);args=parser.parse_args()
    report=audit(args.folder);print(json.dumps({k:v for k,v in report.items() if k!='sourceSHA256'}))
if __name__=='__main__':main()
