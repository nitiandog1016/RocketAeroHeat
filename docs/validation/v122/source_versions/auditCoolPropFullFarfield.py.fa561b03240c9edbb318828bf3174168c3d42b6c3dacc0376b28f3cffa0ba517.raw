"""完整远场向量法向/特征计数及独立通量Jacobian投影；不以边界导入代替域独立性。"""
from pathlib import Path
import hashlib,json,os,sys
import numpy as np
from scipy.io import loadmat
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import buildSU2ThermalRuntime as supervisor
from su2ThermalRuntime import resolve
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    folder=ROOT/'outputs/v122_validation/actual01';target=folder/'farfield_audit.json';assert not target.exists()
    meshpath=ROOT/'outputs/v118_validation/actual01/mesh/hybrid_volume.mat';mesh=loadmat(str(meshpath))
    V=mesh['vertices'];F=mesh['farfield'].astype(int)-1;n=len(V)
    state=np.loadtxt(str(folder/'initial_flow.csv'),delimiter=',',skiprows=1)
    assert state.shape==(n,9) and np.max(abs(state[:,1:4]-V))<1e-14
    d=np.genfromtxt(str(folder/'initial_thermodynamics/DE_output.csv'),names=True,delimiter=',')
    data=np.genfromtxt(str(folder/'initial_thermodynamics/PT_output.csv'),names=True,delimiter=',');inf=data[-1]
    manifest=json.loads((folder/'case_manifest.json').read_text());velocity=np.asarray(manifest['expectedVelocity'])
    external=np.r_[inf['rho_kg_m3'],inf['rho_kg_m3']*velocity,inf['rho_kg_m3']*(inf['e_J_kg']+velocity@velocity/2)]
    faceN=np.cross(V[F[:,1]]-V[F[:,0]],V[F[:,2]]-V[F[:,0]])/2;nodal=np.zeros_like(V)
    for j in range(3):np.add.at(nodal,F[:,j],faceN/3)
    ids=np.unique(F);area=np.linalg.norm(nodal[ids],axis=1);normal=nodal[ids]/area[:,None]
    closure=np.linalg.norm(faceN.sum(axis=0))/np.linalg.norm(faceN,axis=1).sum();assert closure<1e-12
    u=state[ids,5:8]/state[ids,4,None];un=np.sum(u*normal,axis=1);sound=d['a_m_s'][ids]
    lambdas=np.column_stack([un,un,un,un+sound,un-sound]);counts=np.sum(lambdas<0,axis=1)
    assert set(counts)=={0,1,4,5},'完整远场必须同时涵盖超/亚声速法向流入及流出'
    boundary=[];eigerr=0.;projectorerr=0.;zeroerr=0.;allerr=0.
    for j,i in enumerate(ids):
        U=state[i,4:];rho=U[0];v=U[1:4]/rho;energy=U[-1];e=energy/rho-v@v/2;a=sound[j];nu=un[j];nn=normal[j]
        dp=np.r_[d['dpdrho_e'][i]+d['dpde_rho'][i]*(-e+v@v/2)/rho,-d['dpde_rho'][i]*v/rho,d['dpde_rho'][i]/rho]
        dv=np.r_[-nu/rho,nn/rho,0.];A=np.zeros((5,5));A[0,1:4]=nn
        A[1:4]=np.outer(U[1:4],dv)+np.outer(nn,dp);A[1:4,1:4]+=np.eye(3)*nu
        A[4]=(energy+d['p_Pa'][i])*dv+nu*(dp+np.array([0.,0,0,0,1]))
        scale=np.r_[rho,np.full(3,rho*a),rho*a*a];J=A*scale[None,:]/scale[:,None]/a;I=np.eye(5);z=nu/a
        plus=(J-z*I)@(J-(z-1)*I)/2;minus=(J-z*I)@(J-(z+1)*I)/2;conv=I-plus-minus
        P=(z<0)*conv+(z+1<0)*plus+(z-1<0)*minus
        projectorerr=max(projectorerr,float(np.linalg.norm(P@P-P)/max(np.linalg.norm(P),1)))
        if j%41==0:
            ev=np.linalg.eigvals(J);assert max(abs(ev.imag))<1e-8
            eigerr=max(eigerr,float(max(abs(np.sort(ev.real)-np.sort(lambdas[j]/a)))))
        b=U+scale*(P@((external-U)/scale));boundary.append(b)
        if counts[j]==0:zeroerr=max(zeroerr,float(np.max(abs((b-U)/scale))))
        if counts[j]==5:allerr=max(allerr,float(np.max(abs((b-external)/scale))))
    assert max(eigerr,projectorerr,zeroerr,allerr)<1e-8
    boundary=np.array(boundary);assert min(boundary[:,0])>0
    e=boundary[:,4]/boundary[:,0]-.5*np.sum((boundary[:,1:4]/boundary[:,0,None])**2,axis=1)
    out=folder/'boundary_thermodynamics';assert not out.exists();out.mkdir();src=out/'DE_input.csv';dst=out/'DE_output.csv'
    src.write_text('mode,a,b\n'+''.join('DE,%.17g,%.17g\n'%q for q in zip(boundary[:,0],e)))
    runtime=ROOT/'tools/runtime/su2/8.6.0-rah-coolprop-v1';exe=runtime/'bin/RAH_CoolProp_State.exe'
    supervisor.OUT=out;supervisor.run('boundary',[str(exe),str(src),str(dst)],ROOT,dict(os.environ,OMP_NUM_THREADS='1'),120)
    props=np.genfromtxt(str(dst),names=True,delimiter=',');assert len(props)==len(ids)
    native=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/solvers/CEulerSolver.cpp';text=native.read_text()
    block=text[text.index('void CEulerSolver::BC_Riemann('):text.index('void CEulerSolver::BC_Riemann(')+19000]
    assert 'if (Lambda_i[jVar] < 0)' in block and 'Velocity_e[iDim] = Mach[iDim]*GetFluidModel()->GetSoundSpeed();' in block
    resolve(folder,(folder/'case.cfg').read_text())
    report=dict(nodes=n,farfieldNodes=len(ids),farfieldTriangles=len(F),normalVectorClosureRelative=float(closure),
        incomingCharacteristicCounts={str(k):int(sum(counts==k)) for k in [0,1,4,5]},
        independentFluxJacobianEigenvalueError=eigerr,projectorIdempotenceError=projectorerr,
        supersonicOutflowInteriorPreservedError=zeroerr,supersonicInflowExteriorPreservedError=allerr,
        initialBoundaryTemperatureRangeK=[float(min(props['T_K'])),float(max(props['T_K']))],
        allBoundaryStatesValid=True,allNormalsTreatedAsInflow=False,farfieldIndependenceValidated=False,physicsQualified=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),native,meshpath,folder/'case.cfg',folder/'initial_flow.csv',folder/'coolprop_initialization.json',src,dst]})
    target.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='sourceSHA256'},indent=2))
if __name__=='__main__':main()
