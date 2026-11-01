
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
struct Probe {
 int nDim=3,nVar=5;bool implicit=false,UseAccurateJacobian=false;
 double Gamma=1.4,Gamma_Minus_One=.4,Minf=5.,Kp=.25,sigma=1.;
 double V_i[7]{},V_j[7]{},Normal[3]{},UnitNormal[3]{},Velocity_i[3]{},Velocity_j[3]{};
 double psi_i[5]{},psi_j[5]{},Flux[5]{},Pressure_i=0,Pressure_j=0,Density_i=0,Density_j=0;
 double Enthalpy_i=0,Enthalpy_j=0,Area=0,MassFlux=0,DissFlux=0,Pressure=0;
 double** Jacobian_i=nullptr;double** Jacobian_j=nullptr;
 void AccurateJacobian(const CConfig*,double**,double**){throw std::runtime_error("Jacobian outside fragment scope");}
 void ApproximateJacobian(double**,double**){throw std::runtime_error("Jacobian outside fragment scope");}
void ComputeMassAndPressureFluxes(const CConfig* config, su2double &mdot, su2double &pressure) {

  /*--- Projected velocities and squared magnitude ---*/

  su2double ProjVelocity_i = 0.0, ProjVelocity_j = 0.0, sq_vel = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];

    sq_vel += 0.5*(Velocity_i[iDim]*Velocity_i[iDim] + Velocity_j[iDim]*Velocity_j[iDim]);
  }

  /*--- Compute interface speed of sound (aF) ---*/

  su2double astarL = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_i);
  su2double astarR = sqrt(2.0*(Gamma-1.0)/(Gamma+1.0)*Enthalpy_j);

  su2double ahatL = astarL*astarL/max(astarL, ProjVelocity_i);
  su2double ahatR = astarR*astarR/max(astarR,-ProjVelocity_j);

  su2double aF = min(ahatL,ahatR);

  /*--- Left and right pressure functions and Mach numbers ---*/

  su2double mLP, pLP, mRM, pRM;

  su2double mL = ProjVelocity_i/aF;
  su2double mR = ProjVelocity_j/aF;

  su2double MFsq = 0.5*(mL*mL+mR*mR);
  su2double param1 = max(MFsq, Minf*Minf);
  su2double Mrefsq = min(1.0, param1);

  su2double fa = 2.0*sqrt(Mrefsq)-Mrefsq;

  su2double alpha = 3.0/16.0*(-4.0+5.0*fa*fa);
  su2double beta = 1.0/8.0;

  if (fabs(mL) <= 1.0) {
    su2double p1 = 0.25*(mL+1.0)*(mL+1.0);
    su2double p2 = (mL*mL-1.0)*(mL*mL-1.0);

    mLP = p1 + beta*p2;
    pLP = p1*(2.0-mL) + alpha*mL*p2;
  }
  else {
    mLP = 0.5*(mL+fabs(mL));
    pLP = mLP/mL;
  }

  if (fabs(mR) <= 1.0) {
    su2double p1 = 0.25*(mR-1.0)*(mR-1.0);
    su2double p2 = (mR*mR-1.0)*(mR*mR-1.0);

    mRM = -p1 - beta*p2;
    pRM =  p1*(2.0+mR) - alpha*mR*p2;
  }
  else {
    mRM = 0.5*(mR-fabs(mR));
    pRM = mRM/mR;
  }

  /*--- Mass flux with pressure diffusion term ---*/

  su2double rhoF = 0.5*(Density_i+Density_j);
  su2double Mp = -(Kp/fa)*max((1.0-sigma*MFsq),0.0)*(Pressure_j-Pressure_i)/(rhoF*aF*aF);

  su2double mF = mLP + mRM + Mp;
  mdot = aF * (max(mF,0.0)*Density_i + min(mF,0.0)*Density_j);

  /*--- Modified pressure flux ---*/

  pressure = 0.5*(Pressure_j+Pressure_i) + 0.5*(pLP-pRM)*(Pressure_i-Pressure_j) + sqrt(sq_vel)*(pLP+pRM-1.0)*rhoF*aF;

}
CNumerics::ResidualType<> ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  unsigned short iDim, iVar;

  /*--- Space to start preaccumulation ---*/

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+4);
  AD::SetPreaccIn(V_j, nDim+4);

  /*--- Variables for the general form and primitives for mass flux and pressure calculation.  ---*/
  /*--- F_{1/2} = ||A|| ( 0.5 * mdot * (psi_i+psi_j) - 0.5 * |mdot| * (psi_i-psi_j) + N * pf ) ---*/

  psi_i[0] = 1.0;  psi_j[0] = 1.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    /*--- Velocities ---*/
    Velocity_i[iDim] = psi_i[iDim+1] = V_i[iDim+1];
    Velocity_j[iDim] = psi_j[iDim+1] = V_j[iDim+1];
  }

  /*--- Pressure and Density ---*/
  Pressure_i = V_i[nDim+1];  Pressure_j = V_j[nDim+1];
  Density_i  = V_i[nDim+2];  Density_j  = V_j[nDim+2];

  /*--- Enthalpy ---*/
  Enthalpy_i = psi_i[nVar-1] = V_i[nDim+3];
  Enthalpy_j = psi_j[nVar-1] = V_j[nDim+3];

  /*--- Face area (norm or the normal vector) ---*/

  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mass and pressure fluxes defined by derived classes ---*/

  ComputeMassAndPressureFluxes(config, MassFlux, Pressure);
  DissFlux = fabs(MassFlux);

  Flux[0] = MassFlux;

  for (iDim = 0; iDim < nDim; iDim++)
    Flux[iDim+1] = 0.5*MassFlux*(psi_i[iDim+1]+psi_j[iDim+1]) +
                   0.5*DissFlux*(psi_i[iDim+1]-psi_j[iDim+1]) +
                   UnitNormal[iDim]*Pressure;

  Flux[nVar-1] = 0.5*MassFlux*(psi_i[nVar-1]+psi_j[nVar-1]) +
                 0.5*DissFlux*(psi_i[nVar-1]-psi_j[nVar-1]);

  for (iVar = 0; iVar < nVar; iVar++) Flux[iVar] *= Area;

  /*--- Space to end preaccumulation ---*/

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  /*--- If required, compute Jacobians, either approximately (Roe) or numerically ---*/

  if (implicit) {
    if (UseAccurateJacobian)
      AccurateJacobian(config, Jacobian_i, Jacobian_j);
    else
      ApproximateJacobian(Jacobian_i, Jacobian_j);
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}
};

int main(){std::cout<<std::setprecision(17)<<std::scientific;int id;CConfig c;
 while(std::cin>>id){Probe p;for(auto& v:p.V_i)std::cin>>v;for(auto& v:p.V_j)std::cin>>v;for(auto& v:p.Normal)std::cin>>v;
  if(!std::cin)throw std::runtime_error("incomplete input");auto r=p.ComputeResidual(&c);std::cout<<id;
  for(int j=0;j<5;++j){if(!std::isfinite(r.v[j]))throw std::runtime_error("nonfinite flux");std::cout<<','<<r.v[j];}std::cout<<'\n';}}
