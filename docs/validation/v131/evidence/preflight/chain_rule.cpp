
#include <cmath>
#include <cstdio>
#include <stdexcept>
using su2double=double;
namespace GeometryToolbox { double DotProduct(int n,const double* a,const double* b) { double r=0; for(int i=0;i<n;i++)r+=a[i]*b[i];return r; } }
struct Nodes { double phi[2]={2.,3.}; const double* GetSolution(unsigned long)const{return phi;} };
struct Residual { double v[2]={0,0}; void AddBlock(unsigned long,const double* phi,double c){for(int j=0;j<2;j++)v[j]+=c*phi[j];} };
struct Matrix { double diag=0; void AddVal2Diag(unsigned long,double x){diag+=x;} };
struct Original { bool Conservative=true; int nDim=3; Nodes n; Nodes* nodes=&n; Residual LinSysRes; Matrix Jacobian;
  inline su2double BoundedScalarBCFlux(unsigned long iPoint, bool implicit, const su2double& density,
                                       const su2double* velocity, const su2double* normal) {
    const su2double edgeMassFlux = density * GeometryToolbox::DotProduct(nDim, velocity, normal);
    LinSysRes.AddBlock(iPoint, nodes->GetSolution(iPoint), -edgeMassFlux);
    if (implicit) Jacobian.AddVal2Diag(iPoint, -edgeMassFlux);
    return edgeMassFlux;
  }
};
struct Patched { bool Conservative=true; int nDim=3; Nodes n; Nodes* nodes=&n; Residual LinSysRes; Matrix Jacobian;
  inline su2double BoundedScalarBCFlux(unsigned long iPoint, bool implicit, const su2double& density,
                                       const su2double* velocity, const su2double* normal) {
    const su2double edgeMassFlux = density * GeometryToolbox::DotProduct(nDim, velocity, normal);
    LinSysRes.AddBlock(iPoint, nodes->GetSolution(iPoint), -edgeMassFlux);
    if (implicit) Jacobian.AddVal2Diag(iPoint, -edgeMassFlux / (Conservative ? density : 1.0));
    return edgeMassFlux;
  }
};

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
