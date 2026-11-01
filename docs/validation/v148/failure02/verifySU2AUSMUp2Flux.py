"""准确原生通量片段与解析守恒校核；不把数值界面尺度当热力学声速。"""
from pathlib import Path
import hashlib,json,os,subprocess,time,zipfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/v148_validation/flux_preflight02'
PLAN=ROOT/'docs/validation/v148/PLAN.md'
COMMIT='7048fd35f56d4bf52e586429d6746a36a10255f8'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,indent=2,allow_nan=False)+'\n').encode('utf-8'))

def body(text,signature):
    start=text.index('{',text.index(signature));depth=1;end=start+1
    while depth:
        depth+=(text[end]=='{')-(text[end]=='}');end+=1
    return text[start:end]

def main():
    assert not OUT.exists();OUT.mkdir(parents=True)
    archive=ROOT/'tools/runtime/su2-thermal-build/su2.zip'
    assert sha(archive)=='c45f77e16af4ce326e0fca3adf21009d61e4ecc0d6cac40e8dae22654e1c83dc'
    files=['SU2_CFD/src/numerics/flow/convection/ausm_slau.cpp','SU2_CFD/src/variables/CNSVariable.cpp',
           'SU2_CFD/include/variables/CEulerVariable.hpp','SU2_CFD/src/solvers/CEulerSolver.cpp']
    sources={};bindings={}
    # 模块一：固定归档绑定实际primitive、含k总焓及标量通量调用路径。
    with zipfile.ZipFile(str(archive)) as z:
        for rel in files:
            raw=z.read('SU2-'+COMMIT+'/'+rel);p=OUT/Path(rel).name;p.write_bytes(raw)
            sources[rel]=raw.decode('utf-8');bindings[rel]=sha(p)
    native=sources[files[0]]
    assert 'GetEnergy(iPoint)-0.5*Velocity2(iPoint) - turb_ke' in sources[files[1]]
    assert '(Solution(iPoint,nVar-1) + Primitive(iPoint, indices.Pressure())) / Solution(iPoint,0)' in sources[files[2]]
    assert 'numerics->SetPrimitive(bad_recon? V_i : Primitive_i,  bad_recon? V_j : Primitive_j)' in sources[files[3]]
    mass=body(native,'void CUpwAUSMPLUSUP2_Flow::ComputeMassAndPressureFluxes(')
    residual=body(native,'CNumerics::ResidualType<> CUpwAUSMPLUS_SLAU_Base_Flow::ComputeResidual(')
    prefix=r'''
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using su2double=double;using std::min;using std::max;
constexpr int EULER_IMPLICIT=1;
struct CConfig {int GetKind_TimeIntScheme() const {return 0;}};
namespace AD {void StartPreacc(){} void EndPreacc(){} template<class... T>void SetPreaccIn(T... ){} template<class... T>void SetPreaccOut(T... ){}}
namespace GeometryToolbox {double Norm(int n,const double* v){double s=0;for(int i=0;i<n;++i)s+=v[i]*v[i];return std::sqrt(s);}}
struct CNumerics {template<class T=double>struct ResidualType {const double* v;ResidualType(double* p,double**,double**):v(p){} };};
// 测试封装只声明原函数所需成员；隐式Jacobian未在本片段测试中执行。
struct Probe : CNumerics {
 int nDim=3,nVar=5;bool implicit=false,UseAccurateJacobian=false;
 double Gamma=1.4,Gamma_Minus_One=.4,Minf=5.,Kp=.25,sigma=1.;
 double V_i[7]{},V_j[7]{},Normal[3]{},UnitNormal[3]{},Velocity_i[3]{},Velocity_j[3]{};
 double psi_i[5]{},psi_j[5]{},Flux[5]{},Pressure_i=0,Pressure_j=0,Density_i=0,Density_j=0;
 double Enthalpy_i=0,Enthalpy_j=0,Area=0,MassFlux=0,DissFlux=0,Pressure=0;
 double** Jacobian_i=nullptr;double** Jacobian_j=nullptr;
 void AccurateJacobian(const CConfig*,double**,double**){throw std::runtime_error("Jacobian outside fragment scope");}
 void ApproximateJacobian(double**,double**){throw std::runtime_error("Jacobian outside fragment scope");}
'''
    cpp=prefix+'void ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) '+mass+'\n'
    cpp+='CNumerics::ResidualType<> ComputeResidual(const CConfig* config) '+residual+'\n};\n'
    cpp+=r'''
int main(){std::cout<<std::setprecision(17)<<std::scientific;int id;CConfig c;
 while(std::cin>>id){Probe p;for(auto& v:p.V_i)std::cin>>v;for(auto& v:p.V_j)std::cin>>v;for(auto& v:p.Normal)std::cin>>v;
  if(!std::cin)throw std::runtime_error("incomplete input");auto r=p.ComputeResidual(&c);std::cout<<id;
  for(int j=0;j<5;++j){if(!std::isfinite(r.v[j]))throw std::runtime_error("nonfinite flux");std::cout<<','<<r.v[j];}std::cout<<'\n';}}
'''
    source=OUT/'native_flux.cpp';source.write_bytes(cpp.encode('utf-8'))
    assert mass in cpp and residual in cpp
    # 模块二：独立理想气体/SST能量定义。界面aF不是节点c，不施加错误相等门槛。
    gamma=1.4;R=287.058;cp=gamma*R/(gamma-1)
    def primitive(rho,T,u,k):
        u=np.asarray(u);p=rho*R*T;H=cp*T+.5*u@u+k
        return np.asarray([T,*u,p,rho,H])
    def exact(v,n):
        mass=v[5]*(v[1:4]@n);return np.r_[mass,mass*v[1:4]+v[4]*n,mass*v[6]]
    normals=[np.array([1.,0,0]),np.array([.3,-.4,.5])]
    cases=[]
    for k in [0.,1e-10,60533.86678573171]:
        for mach in [0.,.2,-.7,5.,-5.]:
            u=mach*np.sqrt(gamma*R*225)*np.array([.8,.6,0.])
            v=primitive(.02,225,u,k)
            for n in normals:cases.append(dict(kind='uniform',left=v.tolist(),right=v.tolist(),normal=n.tolist(),expected=exact(v,n).tolist()))
    for n in normals:
        a=primitive(.02,300,[0,0,0],0);b=primitive(.04,150,[0,0,0],50000)
        cases.append(dict(kind='stationary-contact',left=a.tolist(),right=b.tolist(),normal=n.tolist(),expected=np.r_[0,a[4]*n,0].tolist()))
    a=primitive(.02,250,[1300,50,-20],60000);b=primitive(.07,800,[350,-10,40],13000);n=normals[1]
    reference=len(cases)
    for left,right,normal,kind in [(a,b,n,'nonuniform'),(b,a,-n,'reversed'),(a,b,3*n,'area-scaled')]:
        cases.append(dict(kind=kind,left=left.tolist(),right=right.tolist(),normal=normal.tolist()))
    write(OUT/'inputs.json',cases)
    inputs='\n'.join(' '.join(map(str,[i]+q['left']+q['right']+q['normal'])) for i,q in enumerate(cases))+'\n'
    (OUT/'inputs.txt').write_bytes(inputs.encode())
    compiler=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe';exe=OUT/'native_flux.exe'
    env=dict(os.environ,PATH=str(compiler.parent)+os.pathsep+os.environ['PATH'])
    def run(name,command,limit,input=None):
        start=time.monotonic();r=subprocess.run(command,env=env,input=input,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=limit)
        (OUT/(name+'.stdout')).write_bytes(r.stdout);(OUT/(name+'.stderr')).write_bytes(r.stderr)
        write(OUT/(name+'_process.json'),dict(status='terminal',exitCode=r.returncode,seconds=time.monotonic()-start,maximumSeconds=limit,command=command))
        assert r.returncode==0,r.stderr.decode(errors='replace');return r.stdout
    run('compile',[str(compiler),'-std=c++17','-O2',str(source),'-o',str(exe)],60)
    raw=run('execute',[str(exe)],10,inputs.encode());(OUT/'native_flux.csv').write_bytes(raw)
    data=np.loadtxt(str(OUT/'native_flux.csv'),delimiter=',');assert np.array_equal(data[:,0],np.arange(len(cases)))
    errors=[]
    for i,q in enumerate(cases[:reference]):
        expected=np.asarray(q['expected']);errors.append(float(max(abs(data[i,1:]-expected)/np.maximum(abs(expected),1))))
    reversed_error=float(max(abs(data[reference+1,1:]+data[reference,1:])/np.maximum(abs(data[reference,1:]),1)))
    area_error=float(max(abs(data[reference+2,1:]-3*data[reference,1:])/np.maximum(abs(3*data[reference,1:]),1)))
    assert max(errors+[reversed_error,area_error])<2e-12
    # 模块三：当前检查只授权进入一个小型完整引擎检查，不授予真实整箭资格。
    report=dict(passed=True,sourceCommit=COMMIT,cases=len(cases),uniformAndContactMaxScaledError=max(errors),
        reversalMaxScaledError=reversed_error,areaScalingMaxScaledError=area_error,
        exactNativeFunctionBodies=True,sharedSourceModified=False,wholeNativeSolverExecuted=False,
        implicitJacobianTested=False,thermodynamicSoundSpeedAndNumericalInterfaceScaleDistinguished=True,
        sstEnergyIncludesK=True,physicsQualified=False,newQualifiedArea=0,sourceBindings=bindings,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),PLAN,archive,compiler,source,exe]})
    write(OUT/'report.json',report);print(json.dumps(report))

if __name__=='__main__':main()
