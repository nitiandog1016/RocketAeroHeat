"""保持完整原表面和四面体外域，仅重新分布法向棱柱层；不调用Gmsh。"""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time
import numpy as np
from scipy.io import loadmat,savemat
from scipy.optimize import brentq
from meshSU2HybridVolume import prepare_layers
ROOT=Path(__file__).resolve().parents[1]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))

def build(source,folder,first,layers):
    receipt=read(source/'mesh_process.json');assert receipt['status']=='terminal' and receipt['exitCode']==0
    assert receipt['inputStable'] and receipt['sourceStable']
    for name in ['hybrid_volume.mat','mesh_report.json']:assert sha(source/name)==receipt['outputSHA256'][name]
    assert sha(Path(receipt['input']))==receipt['inputSHA256']
    old=loadmat(str(source/'hybrid_volume.mat'));inp=loadmat(receipt['input'])
    opt=json.loads(str(inp['optionsJSON'][0]));oldH=np.asarray(opt['layerHeights']);n=len(old['sourceVertices']);oldN=len(oldH)-1
    V=old['sourceVertices'];F=old['sourceFaces'].astype(int)-1
    reproduced=prepare_layers(V,F,opt)[0]
    assert np.array_equal(reproduced,old['vertices'][:len(oldH)*n])
    assert 0<first<oldH[1] and 2<=layers<=64 and layers*first<oldH[-1]
    ratio=brentq(lambda g:first*np.sum(g**np.arange(layers))-oldH[-1],1.,2.)
    heights=np.r_[0,np.cumsum(first*ratio**np.arange(layers))];heights[-1]=oldH[-1]
    opt['layerHeights']=heights.tolist()
    P,W,areas,body,alignment,jac=prepare_layers(V,F,opt)
    assert np.array_equal(P[-n:],old['vertices'][oldN*n:(oldN+1)*n])
    P=np.concatenate([P,old['vertices'][(oldN+1)*n:]])
    shift=(layers-oldN)*n
    T=old['tetrahedra'].astype(int)-1;far=old['farfield'].astype(int)-1
    assert min(T.ravel())>=oldN*n and min(far.ravel())>=oldN*n
    Tnew=T+shift;farnew=far+shift
    assert np.array_equal(P[Tnew],old['vertices'][T]) and np.array_equal(P[farnew],old['vertices'][far])
    assert len(P)<=400000 and len(W)+len(Tnew)<=900000, '本次显式400k/900k准备预算超限'
    report=read(source/'mesh_report.json')
    report.update(nodes=len(P),prisms=len(W),minimumPrismJacobian=float(min(jac)),
                  constructionOperation='normal-collar-redistribution-exact-exterior-reuse',
                  gmshRunThisTask=False,sourceMeshDirectory=str(source),normalLayers=layers,normalHeightsM=heights.tolist(),
                  normalGrowthRatio=ratio,firstNormalHeightM=first,sourceNormalLayers=oldN,
                  originalFacesUnchanged=True,interfaceExactlyPreserved=True,exteriorTetrahedraExactlyPreserved=True,
                  nearWallResolutionValidated=False,farfieldAdequacyValidated=False,physicsQualified=False,newQualifiedArea=0)
    savemat(str(folder/'hybrid_volume.mat'),dict(vertices=P,wedges=W+1,tetrahedra=Tnew+1,
            farfield=farnew+1,sourceVertices=V,sourceFaces=F+1),do_compression=True)
    savemat(str(folder/'mesh_input.mat'),dict(V=V,F=F+1,optionsJSON=json.dumps(opt)),do_compression=True)
    (folder/'mesh_report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--target',type=Path,required=True)
    p.add_argument('--first-height',type=float,required=True);p.add_argument('--layers',type=int,required=True);p.add_argument('--worker',action='store_true');a=p.parse_args()
    source=a.source.resolve();target=a.target.resolve();source.relative_to(ROOT/'outputs');target.relative_to(ROOT/'outputs')
    if a.worker:build(source,target,a.first_height,a.layers);return
    assert not target.exists();target.mkdir(parents=True)
    locked=[source/'hybrid_volume.mat',source/'mesh_report.json',source/'mesh_process.json',Path(read(source/'mesh_process.json')['input']),
            Path(__file__).resolve(),Path(__file__).with_name('meshSU2HybridVolume.py')]
    hashes={str(p):sha(p) for p in locked};start=time.monotonic()
    record=dict(status='running',operation='normal-collar-redistribution',gmshRun=False,maximumSeconds=120,
                sourceHashes=hashes,threads=2,supervisorPID=os.getpid())
    receipt=target/'mesh_process.json';receipt.write_text(json.dumps(record,indent=2)+'\n')
    with (target/'redistribution.log').open('wb') as log:
        command=[sys.executable,str(Path(__file__).resolve()),'--worker','--source',str(source),'--target',str(target),
                 '--first-height',str(a.first_height),'--layers',str(a.layers)]
        proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW,
                              env=dict(os.environ,OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1'))
        record['processID']=proc.pid;timed=False
        try:code=proc.wait(timeout=120)
        except subprocess.TimeoutExpired:timed=True;proc.kill();code=proc.wait()
    stable={k:sha(Path(k))==v for k,v in hashes.items()}
    record.update(status='terminal',exitCode=code,timedOut=timed,elapsedSeconds=time.monotonic()-start,
                  sourceStable=all(stable.values()),sourceStableAfter=stable,inputStable=all(stable.values()),
                  input=str(target/'mesh_input.mat'),inputSHA256=sha(target/'mesh_input.mat') if code==0 else None,
                  scriptSHA256=hashes[str(Path(__file__).resolve())])
    record['outputSHA256']={p.name:sha(p) for p in target.iterdir() if p.is_file() and p!=receipt}
    receipt.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({k:record[k] for k in ['status','exitCode','timedOut','elapsedSeconds']}))
    raise SystemExit(124 if timed else code)

if __name__=='__main__':main()
