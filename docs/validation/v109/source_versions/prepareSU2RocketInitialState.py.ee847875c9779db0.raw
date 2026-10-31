"""为已有完整原生三维网格建立满足壁面/远场的平滑初值；仅初始猜测，不是物理解。"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import numpy as np
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[1]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def prepare(source,target):
    source=Path(source).resolve();target=Path(target).resolve();target.relative_to(ROOT/'outputs')
    assert not target.exists();native=source/'native_input'
    receipt=json.loads((native/'engine_process.json').read_text(encoding='utf8'))
    assert receipt['status']=='terminal' and receipt['exitCode']==0 and all(receipt['inputsStable'].values())
    for name in ('mesh.su2','case.cfg'):
        assert sha(native/name)==receipt['inputSHA256'][str(native/name)]
    s=loadmat(str(source/'normal_volume.mat'),squeeze_me=True,struct_as_record=False);v=s['volume'];V=v.vertices
    assert v.meshConstruction.method=='piecewise-parallel-offset-normal-layers-with-miter-corners'
    manifest=json.loads((native/'case_manifest.json').read_text(encoding='utf8'));f=manifest['flight']
    assert len(V)==manifest['nodes']==59122 and len(v.cells)==manifest['cells']==58560
    n=int(v.domainOptions.radialLayers);k=int(v.meshConstruction.normalLayers);assert n==40 and k==24
    nw=len(V)//(n+1);assert nw*(n+1)==len(V)
    layers=V.reshape(n+1,nw,3);wall=layers[0];shell=layers[k]
    distance=np.linalg.norm(layers[:k+1]-wall,axis=2);span=np.linalg.norm(shell-wall,axis=1)
    eta=distance/span;assert np.max(eta)<=1+1e-13 and np.min(eta)>=0
    eta[0]=0;eta[-1]=1;blend=np.ones((n+1,nw));blend[:k+1]=eta**2*(3-2*eta)
    assert np.all(blend>=0) and np.all(blend<=1) and np.all(np.diff(blend,axis=0)>=0)
    b=blend.ravel();U=np.asarray(manifest['expectedVelocity']);T=f['wallTemperature']+(f['temperature']-f['wallTemperature'])*b
    velocity=b[:,None]*U;R=287.05287;gamma=1.4;rho=f['pressure']/(R*T)
    momentum=rho[:,None]*velocity;energy=rho*(R*T/(gamma-1)+np.sum(velocity**2,axis=1)/2)
    assert np.all(np.isfinite(energy)) and np.min(rho)>0 and np.min(T)>0
    assert np.all(velocity[:nw]==0) and np.all(T[:nw]==f['wallTemperature'])
    assert np.all(velocity[k*nw:]==U) and np.all(T[k*nw:]==f['temperature'])
    target.mkdir(parents=True);shutil.copyfile(str(native/'mesh.su2'),str(target/'mesh.su2'))
    config=(native/'case.cfg').read_text(encoding='utf8')
    assert config.count('RESTART_SOL= NO')==1
    config=config.replace('RESTART_SOL= NO','RESTART_SOL= YES')
    config=re.sub(r'^ITER= \d+$','ITER= 30000',config,flags=re.M)
    config=re.sub(r'^OUTPUT_WRT_FREQ= \d+$','OUTPUT_WRT_FREQ= 250',config,flags=re.M)
    config+='SOLUTION_FILENAME= initial_flow\nREAD_BINARY_RESTART= NO\n'
    (target/'case.cfg').write_bytes(config.encode('utf8'));shutil.copyfile(str(native/'case_manifest.json'),str(target/'source_case_manifest.json'))
    manifest['flight']['iterations']=30000;manifest['flight']['outputEvery']=250
    manifest['initialization']='smooth-wall-to-freestream-guess; not-a-physical-boundary-layer'
    (target/'case_manifest.json').write_bytes((json.dumps(manifest,indent=2)+'\n').encode('utf8'))
    data=np.column_stack((np.arange(len(V)),V,rho,momentum,energy))
    header='"PointID","x","y","z","Density","Momentum_x","Momentum_y","Momentum_z","Energy"'
    np.savetxt(str(target/'initial_flow.csv'),data,delimiter=',',fmt=['%d']+['%.17e']*8,header=header,comments='')
    files=[source/'normal_volume.mat',native/'mesh.su2',native/'case.cfg',native/'case_manifest.json',native/'engine_process.json',
           Path(__file__).resolve(),target/'case.cfg',target/'mesh.su2',target/'initial_flow.csv']
    report=dict(method='C1 cubic blend on existing normal shell; constant initial pressure; ideal-gas conservative state',
                initialGuessOnly=True,boundaryLayerPredicted=False,normalShellLayers=k,radialLayers=n,wallNodes=nw,nodes=len(V),
                wallNoSlip=True,wallTemperature=f['wallTemperature'],outerExactFreestream=True,minimumDensity=float(min(rho)),
                temperatureRange=[float(min(T)),float(max(T))],qualifiedArea=0,
                inputSHA256={str(p):sha(p) for p in files})
    (target/'initialization_manifest.json').write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('target');args=p.parse_args()
    r=prepare(args.source,args.target);print(json.dumps({k:r[k] for k in ('nodes','wallNodes','temperatureRange','minimumDensity')}))
