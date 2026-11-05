#pragma once
using su2double=double;
namespace AD {inline bool BeginPassive(){return false;} inline void EndPassive(bool){}}
struct NativeGeneralJacobian {
  unsigned short nDim=3;
  void GetInviscidProjJac(const double* val_velocity,const double* val_enthalpy,
      const double* val_chi,const double* val_kappa,const double* val_normal,
      double val_scale,double** val_Proj_Jac_Tensor) const
{
  const bool wasActive = AD::BeginPassive();
  unsigned short iDim, jDim;
  su2double sqvel, proj_vel, phi, a1, a2;

  sqvel = 0.0; proj_vel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_velocity[iDim]*val_velocity[iDim];
    proj_vel += val_velocity[iDim]*val_normal[iDim];
  }

  phi = *val_chi + 0.5*sqvel*(*val_kappa);
  a1 = *val_enthalpy;
  a2 = *val_kappa;

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
  val_Proj_Jac_Tensor[nDim+1][nDim+1] = val_scale*(a2+1)*proj_vel;
  AD::EndPassive(wasActive);
}
};
