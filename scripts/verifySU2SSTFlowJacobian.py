"""用准确原生函数和独立守恒差分核对SST能量线性化，不启动整箭求解。"""
from pathlib import Path
import json, os, subprocess, time, zipfile
import numpy as np
from verifySU2AUSMUp2Flux import ROOT, COMMIT, sha, write, body
from su2SSTInterface import state

BASE = ROOT / 'outputs/v149_validation'
OUT = BASE / 'jacobian01'
PLAN = ROOT / 'docs/validation/v149/PLAN.md'
GAMMA = 1.4
GAS_R = 287.058


def primitive(U, k):
    """冻结比湍动能k，能量变量仍为包含k的总能量密度。"""
    rho = U[0]; u = U[1:4] / rho
    p = (GAMMA-1) * (U[4] - .5*np.dot(U[1:4], U[1:4])/rho - rho*k)
    return np.r_[p/(rho*GAS_R), u, p, rho, (U[4]+p)/rho]


def physical_flux(U, k, normal):
    v = primitive(U, k); mdot = np.dot(U[1:4], normal)
    return np.r_[mdot, mdot*v[1:4]+v[4]*normal, mdot*v[6]]


def scales(U, k, normal):
    v = primitive(U, k); speed = np.linalg.norm(v[1:4])+np.sqrt(GAMMA*v[4]/v[5])
    mass = v[5]*speed*np.linalg.norm(normal)
    return np.r_[mass, np.full(3, mass*speed), mass*abs(v[6])], np.r_[v[5], np.full(3, v[5]*speed), v[5]*abs(v[6])]


def synthetic(rho, temp, velocity, k):
    u = np.array(velocity); E = GAS_R*temp/(GAMMA-1)+.5*u@u+k
    return np.r_[rho, rho*u, rho*E]


def real_cases():
    """选既存极值及一个实际单元边邻点；诊断法向不是原生双对偶面法向。"""
    comparison = json.loads((ROOT/'outputs/v148_validation/flux_comparison.json').read_text())
    selected = {name: int(comparison['variants'][name]['enthalpyExtremeNodeZeroBased']) for name in ['control','up2']}
    edges = {10: [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],
             13: [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]}
    near = {v: set() for v in selected.values()}
    mesh = ROOT/'outputs/v148_validation/control/mesh.su2'
    with mesh.open() as f:
        assert f.readline().strip() == 'NDIME= 3'
        count = int(f.readline().split('=')[1])
        for _ in range(count):
            row = [int(s) for s in f.readline().split()]; typ = row[0]; ids = row[1:]
            if not set(near).intersection(ids): continue
            for a,b in edges[typ]:
                if ids[a] in near: near[ids[a]].add(ids[b])
                if ids[b] in near: near[ids[b]].add(ids[a])
    cases = []; locks = {str(mesh): sha(mesh)}
    for name, idx in selected.items():
        folder = ROOT/'outputs/v148_validation'/name
        r = json.loads((folder/'engine_process.json').read_text())
        assert r['status']=='terminal' and r['exitCode']==0
        path = folder/'restart_flow.csv'; assert sha(path)==r['outputSHA256']['restart_flow.csv']
        a, keys, nd, rho, u, _ = state(path, True, True)
        # 固定规则取最小全局编号邻点；不挑选误差最大或最利于结论的邻点。
        j = min(near[idx]); pair = [idx,j]; U = a[pair,4:9]
        p = a[pair,keys.index('Pressure')]
        k = U[:,4]/rho[pair]-.5*np.sum(u[pair]**2,axis=1)-p/((GAMMA-1)*rho[pair])
        assert min(k)>0, '代表状态k必须可解析；不裁剪或改选'
        for t in range(2):
            v = primitive(U[t],k[t]); assert np.max(abs(v[1:4]-u[pair[t]]))<1e-9
            assert abs(v[4]/p[t]-1)<2e-12
        cases.append(dict(name='actual_'+name, kind='nonuniform', left=U[0].tolist(), right=U[1].tolist(),
                          k=k.tolist(), normal=[.3,-.4,.5], pointIDs=pair, coordinates=a[pair,1:4].tolist(),
                          savedUpdatedK=a[pair,keys.index('Turb_Kin_Energy')].tolist(),
                          normalIsNativeDualFace=False, folder=str(folder)))
        locks[str(path)] = sha(path); locks[str(folder/'engine_process.json')] = sha(folder/'engine_process.json')
    return cases, locks


def build_source(sources):
    """模块一：成员封装与原函数体分离，原数值语句逐字保留。"""
    flow = sources['ausm_slau.cpp']; numerics = sources['CNumerics.cpp']; header = sources['CNumerics.hpp']
    parts = {}
    def method(key, text, signature, declaration):
        fragment = body(text, signature); parts[key] = fragment
        return declaration+' '+fragment+'\n'
    common = r'''
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using su2double=double; using std::max; using std::min;
constexpr int EULER_IMPLICIT=1;
struct CConfig {int GetKind_TimeIntScheme() const{return 1;}};
namespace AD {bool BeginPassive(){return false;} void EndPassive(bool){} void StartPreacc(){} void EndPreacc(){}
 template<class... T>void SetPreaccIn(T... ){} template<class... T>void SetPreaccOut(T... ){}}
namespace GeometryToolbox {
 double DotProduct(int n,const double*a,const double*b){double s=0;for(int i=0;i<n;++i)s+=a[i]*b[i];return s;}
 double SquaredNorm(int n,const double*a){return DotProduct(n,a,a);}
 double Norm(int n,const double*a){return std::sqrt(SquaredNorm(n,a));}}
struct CNumerics {template<class T=double>struct ResidualType {const double*v;
 ResidualType(double*p,double**,double**):v(p){}};};
struct Common:CNumerics {
 int nDim=3,nVar=5; double Gamma=1.4,Gamma_Minus_One=Gamma-1;
 bool implicit=true,UseAccurateJacobian=false,HasAnalyticalDerivatives=false;
 double FinDiffStep=1e-4,Minf=5.,Kp=.25,sigma=1.;
 double V_i[7]{},V_j[7]{},Normal[3]{},UnitNormal[3]{},Velocity_i[3]{},Velocity_j[3]{},RoeVelocity[3]{};
 double psi_i[5]{},psi_j[5]{},Flux[5]{},Lambda[5]{},Pressure_i=0,Pressure_j=0,Density_i=0,Density_j=0;
 double Enthalpy_i=0,Enthalpy_j=0,Area=0,MassFlux=0,DissFlux=0,Pressure=0;
 double pm[5][5]{},ipm[5][5]{},ji[5][5]{},jj[5][5]{};
 double* P_Tensor[5]; double* invP_Tensor[5]; double* Jacobian_i[5]; double* Jacobian_j[5];
 double dmdot_dVi[6]{},dmdot_dVj[6]{},dpres_dVi[6]{},dpres_dVj[6]{};
 Common(){for(int i=0;i<5;++i){P_Tensor[i]=pm[i];invP_Tensor[i]=ipm[i];Jacobian_i[i]=ji[i];Jacobian_j[i]=jj[i];}}
'''
    common += method('physicalJacobian',numerics,'void CNumerics::GetInviscidProjJac(const su2double *val_velocity, const su2double *val_energy,',
                     'void GetInviscidProjJac(const double* val_velocity,const double* val_energy,const double* val_normal,double val_scale,double** val_Proj_Jac_Tensor) const')
    for key, name, output in [('P','GetPMatrix','p_tensor'),('Pinv','GetPMatrix_inv','inv_p_tensor')]:
        common += method(key,header,'void '+name+'(const su2double& density,',
                         'template<class Matrix> void '+name+'(const double& density,const double* velocity,const double& soundspeed,const double* normal,Matrix& '+output+') const')
    common += '};\nstruct AUSM:Common {\n'
    common += 'double sq_vel=0,Proj_ModJac_Tensor_ij=0,Energy_i=0,SoundSpeed_i=0,Energy_j=0,SoundSpeed_j=0,R=0,RoeDensity=0,RoeEnthalpy=0,RoeSoundSpeed=0,ProjVelocity=0,ProjVelocity_i=0,ProjVelocity_j=0;\n'
    common += 'unsigned short iDim=0,iVar=0,jVar=0,kVar=0;double mL=0,mR=0,mLP=0,mRM=0,mF=0,pLP=0,pRM=0,pF=0,Phi=0;\n'
    common += method('AUSMResidual',flow,'CNumerics::ResidualType<> CUpwAUSM_Flow::ComputeResidual(',
                     'CNumerics::ResidualType<> ComputeResidual(const CConfig* config)')
    common += '};\nstruct UP2:Common {\n'
    common += method('UP2MassPressure',flow,'void CUpwAUSMPLUSUP2_Flow::ComputeMassAndPressureFluxes(',
                     'void ComputeMassAndPressureFluxes(const CConfig* config,double& mdot,double& pressure)')
    common += method('RoeApproximate',flow,'void CUpwAUSMPLUS_SLAU_Base_Flow::ApproximateJacobian(',
                     'void ApproximateJacobian(double** val_Jacobian_i,double** val_Jacobian_j)')
    common += method('UP2Accurate',flow,'void CUpwAUSMPLUS_SLAU_Base_Flow::AccurateJacobian(',
                     'void AccurateJacobian(const CConfig* config,double** val_Jacobian_i,double** val_Jacobian_j)')
    common += method('UP2Residual',flow,'CNumerics::ResidualType<> CUpwAUSMPLUS_SLAU_Base_Flow::ComputeResidual(',
                     'CNumerics::ResidualType<> ComputeResidual(const CConfig* config)')
    common += r'''
};
template<class P>void run(int id,int mode){P p; CConfig c;p.UseAccurateJacobian=mode==2;
 for(double&v:p.V_i)std::cin>>v;for(double&v:p.V_j)std::cin>>v;for(double&v:p.Normal)std::cin>>v;
 if(!std::cin)throw std::runtime_error("incomplete inputs");p.ComputeResidual(&c);std::cout<<id;
 for(double v:p.Flux)std::cout<<','<<v;
 for(int t=0;t<2;++t)for(int i=0;i<5;++i)for(int j=0;j<5;++j)std::cout<<','<<(t?p.jj[i][j]:p.ji[i][j]);
 std::cout<<'\n';}
int main(){std::cout<<std::setprecision(17)<<std::scientific;int id,mode;
 while(std::cin>>id>>mode){if(mode==0)run<AUSM>(id,mode);else run<UP2>(id,mode);}}
'''
    for fragment in parts.values(): assert fragment in common
    return common, parts


def main():
    assert not OUT.exists(); OUT.mkdir(parents=True)
    archive = ROOT/'tools/runtime/su2-thermal-build/su2.zip'
    assert sha(archive)=='c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    files = ['SU2_CFD/src/numerics/flow/convection/ausm_slau.cpp','SU2_CFD/src/numerics/CNumerics.cpp',
             'SU2_CFD/include/numerics/CNumerics.hpp','SU2_CFD/src/drivers/CDriver.cpp',
             'SU2_CFD/src/variables/CNSVariable.cpp','SU2_CFD/include/variables/CEulerVariable.hpp']
    sources = {}
    with zipfile.ZipFile(str(archive)) as z:
        for rel in files:
            b=z.read('SU2-'+COMMIT+'/'+rel); p=OUT/Path(rel).name; p.write_bytes(b); sources[p.name]=b.decode('utf-8')
    cpp, parts = build_source(sources); source=OUT/'native_jacobian.cpp';source.write_bytes(cpp.encode('utf-8'))
    cases=[]
    for k in [0.,60000.]:
        for u in [[170.,20.,-10.],[1300.,50.,-20.],[-1300.,-50.,20.]]:
            U=synthetic(.02,250,u,k)
            cases.append(dict(name='uniform_'+str(len(cases)),kind='uniform',left=U.tolist(),right=U.tolist(),k=[k,k],normal=[.3,-.4,.5]))
        cases.append(dict(name='nonuniform_'+str(len(cases)),kind='nonuniform',left=synthetic(.02,250,[1300,50,-20],k).tolist(),
                          right=synthetic(.07,800,[350,-10,40],k/3).tolist(),k=[k,k/3],normal=[.3,-.4,.5]))
    real, locks = real_cases(); cases += real
    # 模块二：先冻结全部样本、扰动和尺度，再编译执行一次；不按结果挑样本。
    requests=[]; indices=[]
    def add(mode,L,R,k,n):
        i=len(requests); requests.append([i,mode]+primitive(L,k[0]).tolist()+primitive(R,k[1]).tolist()+list(n)); return i
    for mode in [0,1,2]:
        for c in cases:
            L=np.array(c['left']);R=np.array(c['right']);k=c['k'];n=c['normal']
            fl,sl=scales(L,k[0],n);fr,sr=scales(R,k[1],n);fs=np.maximum(fl,fr)
            item=dict(mode=mode,case=c['name'],base=add(mode,L,R,k,n),fluxScale=fs.tolist(),uScale=[sl.tolist(),sr.tolist()],fd=[])
            for eps in [1e-5,5e-6]:
                changes=[]
                for side in [0,1]:
                    for col in range(5):
                        ids=[]
                        for sign in [1,-1]:
                            l=L.copy();r=R.copy();arr=l if side==0 else r;scale=sl if side==0 else sr
                            arr[col]+=sign*eps*scale[col];ids.append(add(mode,l,r,k,n))
                        changes.append(ids)
                item['fd'].append(dict(epsilon=eps,pairs=changes))
            if c['kind']=='uniform':
                pair=[]
                for col in range(5):
                    ids=[]
                    for sign in [1,-1]:
                        l=L.copy();l[col]+=sign*5e-6*sl[col];ids.append(add(mode,l,l,k,n))
                    pair.append(ids)
                item['uniformPairs']=pair
            indices.append(item)
    write(OUT/'cases.json',cases);write(OUT/'requests.json',indices)
    inputs='\n'.join(' '.join(map(str,row)) for row in requests)+'\n';(OUT/'native_inputs.txt').write_bytes(inputs.encode())
    compiler=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe'; exe=OUT/'native_jacobian.exe'
    env=dict(os.environ,PATH=str(compiler.parent)+os.pathsep+os.environ['PATH'])
    def run(label,command,limit,inp=None):
        start=time.monotonic();r=subprocess.run(command,input=inp,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                               env=env,timeout=limit,creationflags=subprocess.CREATE_NO_WINDOW)
        (OUT/(label+'.stdout')).write_bytes(r.stdout);(OUT/(label+'.stderr')).write_bytes(r.stderr)
        write(OUT/(label+'_process.json'),dict(status='terminal',exitCode=r.returncode,seconds=time.monotonic()-start,
                                             maximumSeconds=limit,command=command))
        assert r.returncode==0,r.stderr.decode(errors='replace');return r.stdout
    run('compile',[str(compiler),'-std=c++17','-O2',str(source),'-o',str(exe)],60)
    (OUT/'native_result.csv').write_bytes(run('execute',[str(exe)],20,inputs.encode()))
    data=np.loadtxt(str(OUT/'native_result.csv'),delimiter=',');assert np.array_equal(data[:,0],np.arange(len(requests))) and np.isfinite(data).all()
    # 模块三：复步物理导数与原生中心差分独立复核，近似矩阵只量化而不误判成精确矩阵。
    reports=[];cmap={c['name']:c for c in cases}
    for item in indices:
        c=cmap[item['case']];L=np.array(c['left']);R=np.array(c['right']);k=c['k'];n=np.array(c['normal'])
        fs=np.array(item['fluxScale']);us=np.array(item['uScale']);base=data[item['base']];J=base[6:].reshape(2,5,5)
        deriv=[]
        for step in item['fd']:
            d=np.array([(data[a,1:6]-data[b,1:6])/(2*step['epsilon']) for a,b in step['pairs']]).reshape(2,5,5).transpose(0,2,1)
            deriv.append(d/fs[None,:,None])
        scaledJ=J*us[:,None,:]/fs[None,:,None]
        row=dict(mode=item['mode'],case=c['name'],kind=c['kind'],centeredDifferenceChange=float(np.max(abs(deriv[1]-deriv[0]))),
                 nativeMatrixVsNumericalFluxScaledError=float(np.max(abs(scaledJ-deriv[1]))))
        if c['kind']=='uniform':
            exact=physical_flux(L,k[0],n);uniformError=float(np.max(abs(base[1:6]-exact)/fs));assert uniformError<2e-12
            oracle=np.empty((5,5))
            for col in range(5):
                z=L.astype(complex);z[col]+=1j*1e-25*us[0,col];oracle[:,col]=np.imag(physical_flux(z,k[0],n))/1e-25/fs
            both=np.array([(data[a,1:6]-data[b,1:6])/(1e-5) for a,b in item['uniformPairs']]).T/fs[:,None]
            oracleError=float(np.max(abs(both-oracle)));assert oracleError<2e-8
            row.update(uniformFluxError=uniformError,uniformDifferenceOracleError=oracleError,
                       uniformSumJacobianScaledError=float(np.max(abs(np.sum(scaledJ,axis=0)-oracle))))
        reports.append(row)
    for p in [Path(__file__).resolve(),PLAN,archive,compiler,source,exe,ROOT/'scripts/verifySU2AUSMUp2Flux.py',ROOT/'scripts/su2SSTInterface.py']:
        locks[str(p)]=sha(p)
    result=dict(passed=True,sourceCommit=COMMIT,samples=len(cases),nativeEvaluations=len(requests),
                exactFunctionBodyCount=len(parts),modeNames=['AUSM-Roe-approximate','AUSMPLUSUP2-Roe-approximate','AUSMPLUSUP2-accurate-not-used-in-v148'],
                cases=reports,sourceSHA256=locks,wholeEngineCompiled=False,wholeEngineExecuted=False,
                newGeometry=False,defaultChanged=False,physicsQualified=False,newQualifiedArea=0,
                scope='frozen-k flow-flux matrices; no gradient/limiter/boundary/global-SST Jacobian or actual dual-normal replay')
    write(OUT/'report.json',result);print(json.dumps(dict(passed=True,samples=len(cases),cases=reports)))


if __name__=='__main__':main()
