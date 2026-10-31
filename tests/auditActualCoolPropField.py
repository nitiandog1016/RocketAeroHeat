"""独立核对真实原面映射、HEOS法向热梯度及两种面积积分，不将启动场赋资格。"""
from pathlib import Path
import hashlib,json,re,sys
import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tests'))
from auditSU2V110 import vtk_tokens
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    folder=ROOT/'outputs/v122_validation/actual04';target=folder/'field_audit.json';assert not target.exists()
    receipt=json.loads((folder/'engine_process.json').read_text());assert receipt['status']=='terminal' and receipt['exitCode']==0 and all(receipt['inputsStable'].values())
    meshpath=ROOT/'outputs/v118_validation/actual01/mesh/hybrid_volume.mat';mesh=loadmat(str(meshpath))
    points=mesh['vertices'];old=loadmat(str(ROOT/'outputs/v119_validation/review01/result.mat'),squeeze_me=True,struct_as_record=False)['r']
    V=old.V;F=old.F.astype(int)-1;n=len(V);assert n==9878 and len(F)==19752 and np.array_equal(V,points[:n])
    vector=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2;area=np.linalg.norm(vector,axis=1);normal=np.zeros_like(V)
    for j in range(3):np.add.at(normal,F[:,j],vector/3)
    xyz,fields=vtk_tokens(folder/'surface_flow.vtk');dist,ids=cKDTree(V).query(xyz);assert max(dist)<1e-6 and len(set(ids))==n
    node=np.zeros((n,3));node[ids]=np.column_stack([fields[k][:,0] for k in ['Pressure','Heat_Flux','Temperature']]);face=node[F].mean(axis=1)
    rest=np.loadtxt(str(folder/'restart_flow.csv'),delimiter=',',skiprows=1)
    prop=np.genfromtxt(str(folder/'thermodynamics/properties.csv'),names=True,delimiter=',')
    assert np.array_equal(rest[:,0],np.arange(len(points))) and max(abs(rest[:,1:4]-points).ravel())<1e-14
    assert np.array_equal(prop['row'],rest[:,0]) and len(prop)==210914
    eos=max(float(max(abs(prop['p_Pa'][:n]/node[:,0]-1))),float(max(abs(prop['T_K'][:n]/node[:,2]-1))));assert eos<2e-6
    unit=normal/np.linalg.norm(normal,axis=1)[:,None];delta=points[n:2*n]-V;distance=abs(np.sum(delta*unit,axis=1))
    assert min(distance)>0;T=prop['T_K'];mu=1.458e-6*T[:n]**1.5/(T[:n]+110.4)
    q=mu*prop['cp_J_kgK'][:n]/.72*(T[n:2*n]-300)/distance
    error=float(max(abs(q-node[:,1])/np.maximum(abs(q),1)));assert error<2e-6
    hist=np.atleast_1d(np.genfromtxt(str(folder/'history.csv'),names=True,delimiter=','));saved=int(hist['Inner_Iter'][-1]);assert saved==19
    for name,ext in [('surface_flow','vtk'),('restart_flow','csv')]:assert sha(folder/(name+'.'+ext))==sha(folder/(name+'_000019.'+ext))
    native=float(np.linalg.norm(normal,axis=1)@node[:,1]);assert abs(native-hist['HF'][-1])/max(abs(native),1)<2e-6
    counts=re.findall(r'RAH_EOS_RECON_COUNT inner=(\d+) thread=(\d+) edges=(\d+)',(folder/'engine_stdout.log').read_text())
    perstep={}
    for iteration,thread,edges in counts:perstep[iteration]=perstep.get(iteration,0)+int(edges)
    assert len(perstep)==20 and perstep['0']>0
    midxyz,midfields=vtk_tokens(folder/'surface_flow_000010.vtk');_,midids=cKDTree(V).query(midxyz)
    mid=np.zeros(n);mid[midids]=midfields['Heat_Flux'][:,0];midface=mid[F].mean(axis=1)
    report=dict(passed=True,nodes=len(points),originalFaces=len(F),areaM2=float(sum(area)),wallEOSRelativeError=eos,
        normalHeatGradientRelativeError=error,nativeHeatPowerW=native,originalFaceHeatPowerW=float(area@face[:,1]),
        originalAreaPreserved=True,iteration=19,sourceCumulativeIteration=799,newPhysicsLocalIteration=19,
        residuals=[float(hist[k][-1]) for k in ['rmsRho','rmsRhoE']],heatChange10To19RelativeL1=float(area@abs(face[:,1]-midface)/(area@abs(midface))),
        temperatureRangeK=[float(min(T)),float(max(T))],reconstructionFallbacksByStep=perstep,
        thermodynamicModel='CoolProp-7.2.0-HEOS-Air',actualTransport='Sutherland; k=mu*cp(T,rho)/Pr0.72',
        givenWallTemperatureK=300,solvedWallTemperatureKnown=False,physicsQualified=False,newQualifiedArea=0,
        sources={str(p):sha(p) for p in [Path(__file__).resolve(),meshpath,folder/'engine_process.json',folder/'case.cfg',folder/'restart_flow.csv',folder/'surface_flow.vtk',folder/'thermodynamics/properties.csv']})
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sources'},indent=2))
if __name__=='__main__':main()
