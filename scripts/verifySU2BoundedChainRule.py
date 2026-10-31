"""从固定原生源码提取散度修正函数，编译运行独立密度链式法则核对。"""
from pathlib import Path
import difflib,hashlib,json,os,subprocess,zipfile
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'tools/runtime/su2-thermal-build/source'
OUT=ROOT/'outputs/v131_validation/preflight'
DEST=ROOT/'third_party/su2_bounded'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    assert not OUT.exists();OUT.mkdir(parents=True);DEST.mkdir(parents=True,exist_ok=True)
    files=['SU2_CFD/include/solvers/CScalarSolver.hpp','SU2_CFD/include/solvers/CScalarSolver.inl']
    originals={};patches={};bindings=[]
    with zipfile.ZipFile(str(ROOT/'tools/runtime/su2-thermal-build/su2.zip')) as z:
        for rel in files:
            original=z.read(next(n for n in z.namelist() if n.endswith('/'+rel)))
            assert (SOURCE/rel).read_bytes()==original
            text=original.decode('utf-8')
            if rel.endswith('.hpp'):
                old='if (implicit) Jacobian.AddVal2Diag(iPoint, -edgeMassFlux);'
                new='if (implicit) Jacobian.AddVal2Diag(iPoint, -edgeMassFlux / (Conservative ? density : 1.0));'
                assert text.count(old)==1;patched=text.replace(old,new)
            else:
                replacements=[
                  ('Jacobian.AddVal2Diag(iPoint, -massFlux);','Jacobian.AddVal2Diag(iPoint, -massFlux / (Conservative ? flowNodes->GetDensity(iPoint) : 1.0));'),
                  ('Jacobian.AddVal2Diag(jPoint, massFlux);','Jacobian.AddVal2Diag(jPoint, massFlux / (Conservative ? flowNodes->GetDensity(jPoint) : 1.0));'),
                  ('Jacobian.AddVal2Diag(iPoint, -divergence);','Jacobian.AddVal2Diag(iPoint, -divergence / (Conservative ? flowNodes->GetDensity(iPoint) : 1.0));')]
                patched=text
                for old,new in replacements:
                    assert patched.count(old)==1;patched=patched.replace(old,new)
            originals[rel]=text;patches[rel]=patched
            bindings.append(dict(path=rel,originalSHA256=hashlib.sha256(original).hexdigest(),patchedSHA256=hashlib.sha256(patched.encode()).hexdigest()))
    diff=''.join(''.join(difflib.unified_diff(originals[p].splitlines(True),patches[p].splitlines(True),fromfile='a/'+p,tofile='b/'+p)) for p in files)
    (DEST/'bounded_conservative_jacobian.diff').write_bytes(diff.encode())
    for rel in files:
        (OUT/(Path(rel).name+'.original')).write_bytes(originals[rel].encode())
        (OUT/(Path(rel).name+'.patched')).write_bytes(patches[rel].encode())
    def body(text):
        start=text.index('  inline su2double BoundedScalarBCFlux(')
        end=text.index('\n  }',start)+4
        return text[start:end]
    prefix=r'''
#include <cmath>
#include <cstdio>
#include <stdexcept>
using su2double=double;
namespace GeometryToolbox { double DotProduct(int n,const double* a,const double* b) { double r=0; for(int i=0;i<n;i++)r+=a[i]*b[i];return r; } }
struct Nodes { double phi[2]={2.,3.}; const double* GetSolution(unsigned long)const{return phi;} };
struct Residual { double v[2]={0,0}; void AddBlock(unsigned long,const double* phi,double c){for(int j=0;j<2;j++)v[j]+=c*phi[j];} };
struct Matrix { double diag=0; void AddVal2Diag(unsigned long,double x){diag+=x;} };
'''
    def cls(name,code):
        return 'struct '+name+' { bool Conservative=true; int nDim=3; Nodes n; Nodes* nodes=&n; Residual LinSysRes; Matrix Jacobian;\n'+code+'\n};\n'
    tail=r'''
template<class T> void check(const char* name,bool repaired) {
  for(double rho:{0.02,1.0,5.0}) {
    T t;double u[3]={2,0,0},n[3]={3,0,0};
    t.BoundedScalarBCFlux(0,true,rho,u,n);
    double numerical=0;double h=rho*1e-5;
    for(int s:{-1,1}) {
      T d;d.n.phi[0]=(rho*2+s*h)/rho;d.BoundedScalarBCFlux(0,false,rho,u,n);
      numerical+=s*d.LinSysRes.v[0]/(2*h);
    }
    const double ratio=t.Jacobian.diag/numerical;
    const bool correct=std::abs(ratio-1)<1e-9;
    if(correct!=(repaired||rho==1.0))throw std::runtime_error("unexpected derivative result");
    std::printf("%s,%.17g,%.17g,%.17g,%.17g\n",name,rho,t.Jacobian.diag,numerical,ratio);
  }
}
int main(){check<Original>("native",false);check<Patched>("patched",true);}
'''
    cpp=OUT/'chain_rule.cpp'
    cpp.write_bytes((prefix+cls('Original',body(originals[files[0]]))+cls('Patched',body(patches[files[0]]))+tail).encode())
    compiler=ROOT/'tools/runtime/su2-thermal-build/toolchain/w64devkit/bin/g++.exe'
    exe=OUT/'chain_rule.exe';env=dict(os.environ,PATH=str(compiler.parent)+os.pathsep+os.environ['PATH'])
    c=subprocess.run([str(compiler),'-std=c++17','-O2',str(cpp),'-o',str(exe)],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
    (OUT/'compile.log').write_bytes(c.stdout);assert c.returncode==0,c.stdout.decode(errors='replace')
    r=subprocess.run([str(exe)],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=10)
    (OUT/'derivatives.csv').write_bytes(b'variant,rho,sourceJacobian,finiteDifference,ratio\n'+r.stdout);assert r.returncode==0
    # 独立单元反例：一个纯对流入流边，未知量为rho*phi，V/dt=1、|massFlux|=1。
    rho=.02;phi=1.;incoming=0.;mass=1.;temporal=1.;residual=mass*(phi-incoming)
    native=(rho*phi-residual/(temporal+mass))/rho
    correct=(rho*phi-residual/(temporal+mass/rho))/rho
    assert native<0 and 0<correct<phi
    report=dict(passed=True,originalMismatchReproduced=True,compiledExactNativeFunctionFragment=True,
        wholeNativeSolverExecuted=False,sourceBindings=bindings,
        independentAdvectionExample=dict(density=rho,nativeImplicitPhi=native,chainRuleCorrectPhi=correct),
        residualOrPhysicsChangedByPatch=False,sourceTreeModified=False,
        nonconservativeSchemesNotQualified=True,productionQualified=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),compiler,cpp,exe,DEST/'bounded_conservative_jacobian.diff']})
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
