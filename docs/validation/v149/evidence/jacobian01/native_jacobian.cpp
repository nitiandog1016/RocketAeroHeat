
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
void GetInviscidProjJac(const double* val_velocity,const double* val_energy,const double* val_normal,double val_scale,double** val_Proj_Jac_Tensor) const {
  const bool wasActive = AD::BeginPassive();
  unsigned short iDim, jDim;
  su2double sqvel, proj_vel, phi, a1, a2;

  sqvel = 0.0; proj_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel    += val_velocity[iDim]*val_velocity[iDim];
    proj_vel += val_velocity[iDim]*val_normal[iDim];
  }

  phi = 0.5*Gamma_Minus_One*sqvel;
  a1 = Gamma*(*val_energy)-phi;
  a2 = Gamma-1.0;

  val_Proj_Jac_Tensor[0][0] = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[0][iDim+1] = val_scale*val_normal[iDim];
  val_Proj_Jac_Tensor[0][nDim+1] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {
    val_Proj_Jac_Tensor[iDim+1][0] = val_scale*(val_normal[iDim]*phi - val_velocity[iDim]*proj_vel);
    for (jDim = 0; jDim < nDim; jDim++)
      val_Proj_Jac_Tensor[iDim+1][jDim+1] = val_scale*(val_normal[jDim]*val_velocity[iDim]-a2*val_normal[iDim]*val_velocity[jDim]);
    val_Proj_Jac_Tensor[iDim+1][iDim+1] += val_scale*proj_vel;
    val_Proj_Jac_Tensor[iDim+1][nDim+1] = val_scale*a2*val_normal[iDim];
  }

  val_Proj_Jac_Tensor[nDim+1][0] = val_scale*proj_vel*(phi-a1);
  for (iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nDim+1][iDim+1] = val_scale*(val_normal[iDim]*a1-a2*val_velocity[iDim]*proj_vel);
  val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*Gamma*proj_vel;
  AD::EndPassive(wasActive);
}
template<class Matrix> void GetPMatrix(const double& density,const double* velocity,const double& soundspeed,const double* normal,Matrix& p_tensor) const {
    const su2double rhooc = density / soundspeed;
    const su2double rhoxc = density * soundspeed;

    if (nDim == 2) {
      const su2double ke = 0.5 * GeometryToolbox::SquaredNorm(2, velocity);
      const su2double projvel = GeometryToolbox::DotProduct(2, velocity, normal);

      p_tensor[0][0] = 1.0;
      p_tensor[0][1] = 0.0;
      p_tensor[0][2] = 0.5 * rhooc;
      p_tensor[0][3] = 0.5 * rhooc;

      p_tensor[1][0] = velocity[0];
      p_tensor[1][1] = density * normal[1];
      p_tensor[1][2] = 0.5 * (velocity[0] * rhooc + normal[0] * density);
      p_tensor[1][3] = 0.5 * (velocity[0] * rhooc - normal[0] * density);

      p_tensor[2][0] = velocity[1];
      p_tensor[2][1] = -density * normal[0];
      p_tensor[2][2] = 0.5 * (velocity[1] * rhooc + normal[1] * density);
      p_tensor[2][3] = 0.5 * (velocity[1] * rhooc - normal[1] * density);

      p_tensor[3][0] = ke;
      p_tensor[3][1] = density * (velocity[0] * normal[1] - velocity[1] * normal[0]);
      p_tensor[3][2] = 0.5 * (ke * rhooc + density * projvel + rhoxc / Gamma_Minus_One);
      p_tensor[3][3] = 0.5 * (ke * rhooc - density * projvel + rhoxc / Gamma_Minus_One);
    } else {
      const su2double ke = 0.5 * GeometryToolbox::SquaredNorm(3, velocity);
      const su2double projvel = GeometryToolbox::DotProduct(3, velocity, normal);

      p_tensor[0][0] = normal[0];
      p_tensor[0][1] = normal[1];
      p_tensor[0][2] = normal[2];
      p_tensor[0][3] = 0.5 * rhooc;
      p_tensor[0][4] = 0.5 * rhooc;

      p_tensor[1][0] = velocity[0] * normal[0];
      p_tensor[1][1] = velocity[0] * normal[1] - density * normal[2];
      p_tensor[1][2] = velocity[0] * normal[2] + density * normal[1];
      p_tensor[1][3] = 0.5 * (velocity[0] * rhooc + density * normal[0]);
      p_tensor[1][4] = 0.5 * (velocity[0] * rhooc - density * normal[0]);

      p_tensor[2][0] = velocity[1] * normal[0] + density * normal[2];
      p_tensor[2][1] = velocity[1] * normal[1];
      p_tensor[2][2] = velocity[1] * normal[2] - density * normal[0];
      p_tensor[2][3] = 0.5 * (velocity[1] * rhooc + density * normal[1]);
      p_tensor[2][4] = 0.5 * (velocity[1] * rhooc - density * normal[1]);

      p_tensor[3][0] = velocity[2] * normal[0] - density * normal[1];
      p_tensor[3][1] = velocity[2] * normal[1] + density * normal[0];
      p_tensor[3][2] = velocity[2] * normal[2];
      p_tensor[3][3] = 0.5 * (velocity[2] * rhooc + density * normal[2]);
      p_tensor[3][4] = 0.5 * (velocity[2] * rhooc - density * normal[2]);

      p_tensor[4][0] = ke * normal[0] + density * (velocity[1] * normal[2] - velocity[2] * normal[1]);
      p_tensor[4][1] = ke * normal[1] + density * (velocity[2] * normal[0] - velocity[0] * normal[2]);
      p_tensor[4][2] = ke * normal[2] + density * (velocity[0] * normal[1] - velocity[1] * normal[0]);
      p_tensor[4][3] = 0.5 * (ke * rhooc + density * projvel + rhoxc / Gamma_Minus_One);
      p_tensor[4][4] = 0.5 * (ke * rhooc - density * projvel + rhoxc / Gamma_Minus_One);
    }
  }
template<class Matrix> void GetPMatrix_inv(const double& density,const double* velocity,const double& soundspeed,const double* normal,Matrix& inv_p_tensor) const {
    const su2double rhoxc = density * soundspeed;
    const su2double c2 = pow(soundspeed, 2);
    const su2double gm1 = Gamma_Minus_One;
    const su2double k0orho = normal[0] / density;
    const su2double k1orho = normal[1] / density;
    const su2double gm1_o_c2 = gm1 / c2;
    const su2double gm1_o_rhoxc = gm1 / rhoxc;

    if (nDim == 3) {
      const su2double k2orho = normal[2] / density;
      const su2double ke = 0.5 * GeometryToolbox::SquaredNorm(3, velocity);
      const su2double projvel_o_rho = GeometryToolbox::DotProduct(3, velocity, normal) / density;

      inv_p_tensor[0][0] = normal[0] + k1orho * velocity[2] - k2orho * velocity[1] - normal[0] * gm1_o_c2 * ke;
      inv_p_tensor[0][1] = normal[0] * gm1_o_c2 * velocity[0];
      inv_p_tensor[0][2] = k2orho + normal[0] * gm1_o_c2 * velocity[1];
      inv_p_tensor[0][3] = -k1orho + normal[0] * gm1_o_c2 * velocity[2];
      inv_p_tensor[0][4] = -normal[0] * gm1_o_c2;

      inv_p_tensor[1][0] = normal[1] + k2orho * velocity[0] - k0orho * velocity[2] - normal[1] * gm1_o_c2 * ke;
      inv_p_tensor[1][1] = -k2orho + normal[1] * gm1_o_c2 * velocity[0];
      inv_p_tensor[1][2] = normal[1] * gm1_o_c2 * velocity[1];
      inv_p_tensor[1][3] = k0orho + normal[1] * gm1_o_c2 * velocity[2];
      inv_p_tensor[1][4] = -normal[1] * gm1_o_c2;

      inv_p_tensor[2][0] = normal[2] + k0orho * velocity[1] - k1orho * velocity[0] - normal[2] * gm1_o_c2 * ke;
      inv_p_tensor[2][1] = k1orho + normal[2] * gm1_o_c2 * velocity[0];
      inv_p_tensor[2][2] = -k0orho + normal[2] * gm1_o_c2 * velocity[1];
      inv_p_tensor[2][3] = normal[2] * gm1_o_c2 * velocity[2];
      inv_p_tensor[2][4] = -normal[2] * gm1_o_c2;

      inv_p_tensor[3][0] = -projvel_o_rho + gm1_o_rhoxc * ke;
      inv_p_tensor[3][1] = k0orho - gm1_o_rhoxc * velocity[0];
      inv_p_tensor[3][2] = k1orho - gm1_o_rhoxc * velocity[1];
      inv_p_tensor[3][3] = k2orho - gm1_o_rhoxc * velocity[2];
      inv_p_tensor[3][4] = gm1_o_rhoxc;

      inv_p_tensor[4][0] = projvel_o_rho + gm1_o_rhoxc * ke;
      inv_p_tensor[4][1] = -k0orho - gm1_o_rhoxc * velocity[0];
      inv_p_tensor[4][2] = -k1orho - gm1_o_rhoxc * velocity[1];
      inv_p_tensor[4][3] = -k2orho - gm1_o_rhoxc * velocity[2];
      inv_p_tensor[4][4] = gm1_o_rhoxc;
    } else {
      const su2double ke = 0.5 * GeometryToolbox::SquaredNorm(2, velocity);
      const su2double projvel_o_rho = GeometryToolbox::DotProduct(2, velocity, normal) / density;

      inv_p_tensor[0][0] = 1 - gm1_o_c2 * ke;
      inv_p_tensor[0][1] = gm1_o_c2 * velocity[0];
      inv_p_tensor[0][2] = gm1_o_c2 * velocity[1];
      inv_p_tensor[0][3] = -gm1_o_c2;

      inv_p_tensor[1][0] = -k1orho * velocity[0] + k0orho * velocity[1];
      inv_p_tensor[1][1] = k1orho;
      inv_p_tensor[1][2] = -k0orho;
      inv_p_tensor[1][3] = 0;

      inv_p_tensor[2][0] = -projvel_o_rho + gm1_o_rhoxc * ke;
      inv_p_tensor[2][1] = k0orho - gm1_o_rhoxc * velocity[0];
      inv_p_tensor[2][2] = k1orho - gm1_o_rhoxc * velocity[1];
      inv_p_tensor[2][3] = gm1_o_rhoxc;

      inv_p_tensor[3][0] = projvel_o_rho + gm1_o_rhoxc * ke;
      inv_p_tensor[3][1] = -k0orho - gm1_o_rhoxc * velocity[0];
      inv_p_tensor[3][2] = -k1orho - gm1_o_rhoxc * velocity[1];
      inv_p_tensor[3][3] = gm1_o_rhoxc;
    }
  }
};
struct AUSM:Common {
double sq_vel=0,Proj_ModJac_Tensor_ij=0,Energy_i=0,SoundSpeed_i=0,Energy_j=0,SoundSpeed_j=0,R=0,RoeDensity=0,RoeEnthalpy=0,RoeSoundSpeed=0,ProjVelocity=0,ProjVelocity_i=0,ProjVelocity_j=0;
unsigned short iDim=0,iVar=0,jVar=0,kVar=0;double mL=0,mR=0,mLP=0,mRM=0,mF=0,pLP=0,pRM=0,pF=0,Phi=0;
CNumerics::ResidualType<> ComputeResidual(const CConfig* config) {

  implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+4);
  AD::SetPreaccIn(V_j, nDim+4);

  /*--- Face area (norm or the normal vector) ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  /*-- Unit Normal ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Primitive variables at point i ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
  }
  Pressure_i = V_i[nDim+1];
  Density_i = V_i[nDim+2];
  Enthalpy_i = V_i[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  SoundSpeed_i = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel)));

  /*--- Primitive variables at point j ---*/
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
  }
  Pressure_j = V_j[nDim+1];
  Density_j = V_j[nDim+2];
  Enthalpy_j = V_j[nDim+3];
  Energy_j = Enthalpy_j - Pressure_j/Density_j;
  SoundSpeed_j = sqrt(fabs(Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel)));

  /*--- Projected velocities ---*/
  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
  }

  mL  = ProjVelocity_i/SoundSpeed_i;
  mR  = ProjVelocity_j/SoundSpeed_j;

  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
  else mLP = 0.5*(mL+fabs(mL));

  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
  else mRM = 0.5*(mR-fabs(mR));

  mF = mLP + mRM;

  if (fabs(mL) <= 1.0) pLP = 0.25*Pressure_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
  else pLP = 0.5*Pressure_i*(mL+fabs(mL))/mL;

  if (fabs(mR) <= 1.0) pRM = 0.25*Pressure_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
  else pRM = 0.5*Pressure_j*(mR-fabs(mR))/mR;

  pF = pLP + pRM;
  Phi = fabs(mF);

  Flux[0] = 0.5*(mF*((Density_i*SoundSpeed_i)+(Density_j*SoundSpeed_j))-Phi*((Density_j*SoundSpeed_j)-(Density_i*SoundSpeed_i)));
  for (iDim = 0; iDim < nDim; iDim++)
    Flux[iDim+1] = 0.5*(mF*((Density_i*SoundSpeed_i*Velocity_i[iDim])+(Density_j*SoundSpeed_j*Velocity_j[iDim]))
                      -Phi*((Density_j*SoundSpeed_j*Velocity_j[iDim])-(Density_i*SoundSpeed_i*Velocity_i[iDim])))+UnitNormal[iDim]*pF;
  Flux[nVar-1] = 0.5*(mF*((Density_i*SoundSpeed_i*Enthalpy_i)+(Density_j*SoundSpeed_j*Enthalpy_j))-Phi*((Density_j*SoundSpeed_j*Enthalpy_j)-(Density_i*SoundSpeed_i*Enthalpy_i)));

  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] *= Area;

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  /*--- Roe's Jacobian for AUSM (this must be fixed) ---*/
  if (implicit) {

    /*--- Mean Roe variables iPoint and jPoint ---*/
    R = sqrt(fabs(Density_j/Density_i));
    RoeDensity = R*Density_i;
    sq_vel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
      sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));

    /*--- Compute P and Lambda (do it with the Normal) ---*/
    GetPMatrix(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitNormal, P_Tensor);

    ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      ProjVelocity   += RoeVelocity[iDim]*UnitNormal[iDim];
      ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
      ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
    }

    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Lambda[iDim] = ProjVelocity;
    Lambda[nVar-2]  = ProjVelocity + RoeSoundSpeed;
    Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;

    /*--- Compute inverse P ---*/
    GetPMatrix_inv(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitNormal, invP_Tensor);

    /*--- Jacobias of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
    GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, Jacobian_i);
    GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, Jacobian_j);

    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      for (jVar = 0; jVar < nVar; jVar++) {
        Proj_ModJac_Tensor_ij = 0.0;
        /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
        for (kVar = 0; kVar < nVar; kVar++)
          Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
        Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
        Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
      }
    }
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}
};
struct UP2:Common {
void ComputeMassAndPressureFluxes(const CConfig* config,double& mdot,double& pressure) {

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
void ApproximateJacobian(double** val_Jacobian_i,double** val_Jacobian_j) {

  unsigned short iDim, iVar, jVar, kVar;
  su2double R, RoeDensity, RoeEnthalpy, RoeSoundSpeed, ProjVelocity, sq_vel, Energy_i, Energy_j;

  Energy_i = Enthalpy_i - Pressure_i/Density_i;
  Energy_j = Enthalpy_j - Pressure_j/Density_j;

  /*--- Mean Roe variables iPoint and jPoint ---*/

  R = sqrt(fabs(Density_j/Density_i));
  RoeDensity = R*Density_i;
  ProjVelocity = 0.0;
  sq_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
    ProjVelocity += RoeVelocity[iDim]*UnitNormal[iDim];
    sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
  }
  RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
  RoeSoundSpeed = sqrt(fabs((Gamma-1)*(RoeEnthalpy-0.5*sq_vel)));

  /*--- Compute P and Lambda (do it with the Normal) ---*/

  GetPMatrix(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitNormal, P_Tensor);

  /*--- Flow eigenvalues and Entropy correctors ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    Lambda[iDim] = ProjVelocity;
  Lambda[nVar-2] = ProjVelocity + RoeSoundSpeed;
  Lambda[nVar-1] = ProjVelocity - RoeSoundSpeed;

  /*--- Compute inverse P ---*/
  GetPMatrix_inv(RoeDensity, RoeVelocity, RoeSoundSpeed, UnitNormal, invP_Tensor);

  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/
  GetInviscidProjJac(Velocity_i, &Energy_i, Normal, 0.5, val_Jacobian_i);
  GetInviscidProjJac(Velocity_j, &Energy_j, Normal, 0.5, val_Jacobian_j);

  /*--- Roe's Flux approximation ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      su2double Proj_ModJac_Tensor_ij = 0.0;
      /*--- Compute |Proj_ModJac_Tensor| = P x |Lambda| x inverse P ---*/
      for (kVar = 0; kVar < nVar; kVar++)
        Proj_ModJac_Tensor_ij += P_Tensor[iVar][kVar]*fabs(Lambda[kVar])*invP_Tensor[kVar][jVar];
      val_Jacobian_i[iVar][jVar] += 0.5*Proj_ModJac_Tensor_ij*Area;
      val_Jacobian_j[iVar][jVar] -= 0.5*Proj_ModJac_Tensor_ij*Area;
    }
  }

}
void AccurateJacobian(const CConfig* config,double** val_Jacobian_i,double** val_Jacobian_j) {

  /*--- Compute Jacobians using a mixed (numerical/analytical) formulation ---*/

  unsigned short iDim, iVar, jVar;

  /*--- If not computed analytically, numerically differentiate the fluxes wrt primitives ---*/

  if (!HasAnalyticalDerivatives) {

    /*--- Create arrays of pointers to the primitive variables so
     we can loop through and perturb them in a general way. ---*/

    su2double *primitives_i[6], *primitives_j[6];

    for (iDim = 0; iDim < nDim; ++iDim) {
      primitives_i[iDim] = &Velocity_i[iDim];
      primitives_j[iDim] = &Velocity_j[iDim];
    }
    primitives_i[ nDim ] = &Pressure_i;  primitives_j[ nDim ] = &Pressure_j;
    primitives_i[nDim+1] = &Density_i;   primitives_j[nDim+1] = &Density_j;
    primitives_i[nDim+2] = &Enthalpy_i;  primitives_j[nDim+2] = &Enthalpy_j;

    /*--- Initialize the gradient arrays with the negative of the quantity,
     then for forward finite differences we add to it and divide. ---*/

    for (iVar = 0; iVar < 6; ++iVar) {
      dmdot_dVi[iVar] = -MassFlux;  dpres_dVi[iVar] = -Pressure;
      dmdot_dVj[iVar] = -MassFlux;  dpres_dVj[iVar] = -Pressure;
    }

    for (iVar = 0; iVar < nDim+3; ++iVar) {
      /*--- Perturb side i ---*/
      su2double epsilon = FinDiffStep * max(1.0, fabs(*primitives_i[iVar]));
      *primitives_i[iVar] += epsilon;
      ComputeMassAndPressureFluxes(config, MassFlux, Pressure);
      dmdot_dVi[iVar] += MassFlux;  dpres_dVi[iVar] += Pressure;
      dmdot_dVi[iVar] /= epsilon;   dpres_dVi[iVar] /= epsilon;
      *primitives_i[iVar] -= epsilon;

      /*--- Perturb side j ---*/
      epsilon = FinDiffStep * max(1.0, fabs(*primitives_j[iVar]));
      *primitives_j[iVar] += epsilon;
      ComputeMassAndPressureFluxes(config, MassFlux, Pressure);
      dmdot_dVj[iVar] += MassFlux;  dpres_dVj[iVar] += Pressure;
      dmdot_dVj[iVar] /= epsilon;   dpres_dVj[iVar] /= epsilon;
      *primitives_j[iVar] -= epsilon;
    }
  }

  /*--- Differentiation of fluxes wrt conservatives assuming ideal gas ---*/

  su2double dmdot_dUi[5], dmdot_dUj[5], dpres_dUi[5], dpres_dUj[5];
  su2double sq_veli = 0.0, sq_velj = 0.0, dHi_drhoi = 0.0, dHj_drhoj = 0.0;
  su2double oneOnRhoi = 1.0/Density_i, oneOnRhoj = 1.0/Density_j;

  for (jVar = 0; jVar < nVar; ++jVar) {

    /*--- Partial derivatives of the primitives wrt conservative "jVar" ---*/
    su2double dVi_dUi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    su2double dVj_dUj[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if (jVar == 0) { // Density
      for (iDim = 0; iDim < nDim; ++iDim) {
        // -u,v,w / rho
        dVi_dUi[iDim] = -Velocity_i[iDim] * oneOnRhoi;
        dVj_dUj[iDim] = -Velocity_j[iDim] * oneOnRhoj;
        // ||V||^2
        sq_veli += Velocity_i[iDim] * Velocity_i[iDim];
        sq_velj += Velocity_j[iDim] * Velocity_j[iDim];
      }
      dVi_dUi[nDim] = 0.5*Gamma_Minus_One*sq_veli;
      dVj_dUj[nDim] = 0.5*Gamma_Minus_One*sq_velj;

      dVi_dUi[nDim+1] = dVj_dUj[nDim+1] = 1.0;

      dHi_drhoi = 0.5*(Gamma-2.0)*sq_veli - Gamma*Pressure_i/((Gamma-1.0)*Density_i);
      dHj_drhoj = 0.5*(Gamma-2.0)*sq_velj - Gamma*Pressure_j/((Gamma-1.0)*Density_j);
      dVi_dUi[nDim+2] = dHi_drhoi * oneOnRhoi;
      dVj_dUj[nDim+2] = dHj_drhoj * oneOnRhoj;
    }
    else if (jVar == nVar-1) { // rho*Energy
      dVi_dUi[nDim] = dVj_dUj[nDim] = Gamma_Minus_One;
      dVi_dUi[nDim+2] = Gamma * oneOnRhoi;
      dVj_dUj[nDim+2] = Gamma * oneOnRhoj;
    }
    else { // Momentum
      dVi_dUi[jVar-1] = oneOnRhoi;
      dVj_dUj[jVar-1] = oneOnRhoj;

      dVi_dUi[nDim] = -Gamma_Minus_One*Velocity_i[jVar-1];
      dVj_dUj[nDim] = -Gamma_Minus_One*Velocity_j[jVar-1];

      dVi_dUi[nDim+2] = dVi_dUi[nDim] * oneOnRhoi;
      dVj_dUj[nDim+2] = dVj_dUj[nDim] * oneOnRhoj;
    }

    /*--- Dot product to complete chain rule ---*/
    dmdot_dUi[jVar] = 0.0;  dpres_dUi[jVar] = 0.0;
    dmdot_dUj[jVar] = 0.0;  dpres_dUj[jVar] = 0.0;
    for (iVar = 0; iVar < 6; ++iVar) {
      dmdot_dUi[jVar] += dmdot_dVi[iVar]*dVi_dUi[iVar];
      dpres_dUi[jVar] += dpres_dVi[iVar]*dVi_dUi[iVar];
      dmdot_dUj[jVar] += dmdot_dVj[iVar]*dVj_dUj[iVar];
      dpres_dUj[jVar] += dpres_dVj[iVar]*dVj_dUj[iVar];
    }
  }

  /*--- Assemble final Jacobians (assuming phi = |mdot|) ---*/

  su2double mdot_hat, psi_hat[5];

  if (MassFlux > 0.0) {
    mdot_hat = Area*MassFlux*oneOnRhoi;
    for (iVar = 0; iVar < nVar; ++iVar) psi_hat[iVar] = Area*psi_i[iVar];
  }
  else {
    mdot_hat = Area*MassFlux*oneOnRhoj;
    for (iVar = 0; iVar < nVar; ++iVar) psi_hat[iVar] = Area*psi_j[iVar];
  }

  /*--- Contribution from the mass flux derivatives ---*/
  for (iVar = 0; iVar < nVar; ++iVar) {
    for (jVar = 0; jVar < nVar; ++jVar) {
      val_Jacobian_i[iVar][jVar] = psi_hat[iVar] * dmdot_dUi[jVar];
      val_Jacobian_j[iVar][jVar] = psi_hat[iVar] * dmdot_dUj[jVar];
    }
  }

  /*--- Contribution from the pressure derivatives ---*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jVar = 0; jVar < nVar; ++jVar) {
      val_Jacobian_i[iDim+1][jVar] += Normal[iDim] * dpres_dUi[jVar];
      val_Jacobian_j[iDim+1][jVar] += Normal[iDim] * dpres_dUj[jVar];
    }
  }

  /*--- Contributions from the derivatives of PSI wrt the conservatives ---*/
  if (MassFlux > 0.0) {
    /*--- Velocity terms ---*/
    for (iDim = 0; iDim < nDim; ++iDim) {
      val_Jacobian_i[iDim+1][0]      -= mdot_hat*Velocity_i[iDim];
      val_Jacobian_i[iDim+1][iDim+1] += mdot_hat;
      val_Jacobian_i[nVar-1][iDim+1] -= mdot_hat*Gamma_Minus_One*Velocity_i[iDim];
    }
    /*--- Energy terms ---*/
    val_Jacobian_i[nVar-1][0]      += mdot_hat*dHi_drhoi;
    val_Jacobian_i[nVar-1][nVar-1] += mdot_hat*Gamma;
  }
  else {
    /*--- Velocity terms ---*/
    for (iDim = 0; iDim < nDim; ++iDim) {
      val_Jacobian_j[iDim+1][0]      -= mdot_hat*Velocity_j[iDim];
      val_Jacobian_j[iDim+1][iDim+1] += mdot_hat;
      val_Jacobian_j[nVar-1][iDim+1] -= mdot_hat*Gamma_Minus_One*Velocity_j[iDim];
    }
    /*--- Energy terms ---*/
    val_Jacobian_j[nVar-1][0]      += mdot_hat*dHj_drhoj;
    val_Jacobian_j[nVar-1][nVar-1] += mdot_hat*Gamma;
  }

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
template<class P>void run(int id,int mode){P p; CConfig c;p.UseAccurateJacobian=mode==2;
 for(double&v:p.V_i)std::cin>>v;for(double&v:p.V_j)std::cin>>v;for(double&v:p.Normal)std::cin>>v;
 if(!std::cin)throw std::runtime_error("incomplete inputs");p.ComputeResidual(&c);std::cout<<id;
 for(double v:p.Flux)std::cout<<','<<v;
 for(int t=0;t<2;++t)for(int i=0;i<5;++i)for(int j=0;j<5;++j)std::cout<<','<<(t?p.jj[i][j]:p.ji[i][j]);
 std::cout<<'\n';}
int main(){std::cout<<std::setprecision(17)<<std::scientific;int id,mode;
 while(std::cin>>id>>mode){if(mode==0)run<AUSM>(id,mode);else run<UP2>(id,mode);}}
