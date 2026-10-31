"""独立核验法向迁移及不变外域后准备两项有预算的真实运行，不启动引擎。"""
from pathlib import Path
import json,hashlib,os,re,sys,uuid
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/v120_validation'
sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import prepare,config,CONTROL
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def main():
    assert not (OUT/'initial_audit.json').exists()
    a=loadmat(str(ROOT/'outputs/v118_validation/actual01/mesh/hybrid_volume.mat'))
    b=loadmat(str(OUT/'mesh32/hybrid_volume.mat'));n=len(a['sourceVertices']);shift=14*n
    assert np.array_equal(a['sourceVertices'],b['sourceVertices']) and np.array_equal(a['sourceFaces'],b['sourceFaces'])
    assert np.array_equal(a['vertices'][18*n:],b['vertices'][32*n:])
    assert np.array_equal(a['tetrahedra']+shift,b['tetrahedra']) and np.array_equal(a['farfield']+shift,b['farfield'])
    old=np.loadtxt(str(ROOT/'outputs/v119_validation/continuation01/restart_flow_000299.csv'),delimiter=',',skiprows=1)
    new=np.loadtxt(str(OUT/'initial_flow.csv'),delimiter=',',skiprows=1)
    assert np.array_equal(new[:,0],np.arange(len(new))) and np.array_equal(new[:,1:4],b['vertices'])
    assert np.array_equal(new[:n,4:],old[:n,4:]) and np.array_equal(new[32*n:,4:],old[18*n:,4:])
    # 不调用迁移实现/其权重：直接以几何投影定位旧线段并重新求比例。
    maximum=0.;coordinate=0.
    for j in range(1,32):
        xyz=b['vertices'][j*n:(j+1)*n];found=np.zeros(n,dtype=bool)
        expected=np.zeros((n,5))
        for k in range(18):
            lo=a['vertices'][k*n:(k+1)*n];hi=a['vertices'][(k+1)*n:(k+2)*n];v=hi-lo
            t=np.sum((xyz-lo)*v,axis=1)/np.sum(v*v,axis=1);mask=(t>=0)&(t<=1)
            if not mask.any():continue
            assert not found[mask].any();found[mask]=True
            coordinate=max(coordinate,float(np.max(abs(lo[mask]+t[mask,None]*v[mask]-xyz[mask]))))
            expected[mask]=old[k*n:(k+1)*n,4:][mask]*(1-t[mask,None])+old[(k+1)*n:(k+2)*n,4:][mask]*t[mask,None]
        assert found.all()
        value=new[j*n:(j+1)*n,4:];maximum=max(maximum,float(np.max(abs(expected-value)/np.maximum(abs(value),1))))
    assert maximum<2e-9 and coordinate<2e-14
    rho=new[:,4];internal=new[:,8]/rho-.5*np.sum((new[:,5:8]/rho[:,None])**2,axis=1)
    assert np.isfinite(new).all() and min(rho)>0 and min(internal)>0
    transfer=read(OUT/'initial_flow.transfer.json');assert transfer['sourceCumulativeIteration']==599 and not transfer['sameMesh']
    for p,h in transfer['sourceSHA256'].items():assert sha(Path(p))==h
    assert sha(OUT/'initial_flow.csv')==transfer['outputSHA256']
    source=ROOT/'outputs/v119_validation/continuation01';control=OUT/'study01/control';candidate=OUT/'study01/candidate'
    assert sha(control/'mesh.su2')==sha(source/'mesh.su2')
    # 只替换本轮刚导出的相同字节副本，旧输入及其硬链接身份不动。
    control.joinpath('mesh.su2').unlink();os.link(str(source/'mesh.su2'),str(control/'mesh.su2'))
    lineage=prepare(source,control);assert lineage['completedUpdatesBeforeRun']==600
    for name in ('initial_flow.csv','initial_flow.transfer.json'):os.link(str(OUT/name),str(candidate/name))
    text=(candidate/'case.cfg').read_text();text=text.replace('RESTART_SOL= NO','RESTART_SOL= YES')+'\nSOLUTION_FILENAME= initial_flow\nREAD_BINARY_RESTART= NO\n'
    (candidate/'case.cfg').write_text(text,encoding='utf-8')
    manifest=read(candidate/'case_manifest.json');manifest['initialization']=transfer
    (candidate/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    ac=config(control/'case.cfg');bc=config(candidate/'case.cfg');assert ac==bc
    assert ac['ITER']=='200' and ac['OUTPUT_WRT_FREQ']=='100' and ac['CONV_RESIDUAL_MINVAL']=='-8'
    requests=[]
    for name in ('control','candidate'):
        folder=OUT/'study01'/name
        files=[folder/'case_manifest.json',folder/'projection.mat',ROOT/'docs/validation/v120/PLAN.md',Path(__file__).resolve(),
               ROOT/'scripts/transferSU2NormalColumns.py',ROOT/'scripts/redistributeSU2NormalCollar.py',source/'engine_process.json',OUT/'initial_audit.json']
        files+=[folder/('continuation.json' if name=='control' else 'initial_flow.transfer.json')]
        request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=720,
                     dryRun=False,jobToken=str(uuid.uuid4()),sourceFiles=[str(p.relative_to(ROOT)) for p in files])
        (folder/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8');requests.append(request)
    report=dict(passed=True,originalFaces=19752,wallAndExteriorExact=True,exteriorCellsExact=True,
                maximumStateScaledDifference=maximum,maximumGeometricInterpolationErrorM=coordinate,
                minimumDensity=float(min(rho)),minimumInternalEnergy=float(min(internal)),sourceCumulativeIteration=599,
                controlCompletedUpdatesBeforeRun=600,candidateSameGridContinuation=False,strictlyConservativeRemap=False,
                checkpointBudgetMiB=512,physicsQualified=False,newQualifiedArea=0,
                sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),OUT/'initial_flow.transfer.json',OUT/'study01/preparation.json']})
    (OUT/'initial_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
