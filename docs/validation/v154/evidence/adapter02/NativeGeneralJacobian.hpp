#pragma once
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
struct NativeGeneralHeatJacobian {
  unsigned short nDim=3;
  double Prandtl_Turb=.9,heat_flux_jac_i[5]{};
  void SetHeatFluxJacobian(const double* val_Mean_PrimVar,const double* val_Mean_SecVar,
      double val_eddy_viscosity,double val_thermal_conductivity,
      double val_heat_capacity_cp,double val_dist_ij)
{
  /* Viscous flux Jacobians for arbitrary equations of state */

  //order of val_mean_primitives: T, vx, vy, vz, P, rho, ht
  //order of secondary:dTdrho_e, dTde_rho

  su2double sqvel = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  su2double rho = val_Mean_PrimVar[nDim+2];
  su2double P= val_Mean_PrimVar[nDim+1];
  su2double h= val_Mean_PrimVar[nDim+3];
  su2double dTdrho_e= val_Mean_SecVar[0];
  su2double dTde_rho= val_Mean_SecVar[1];

  su2double dTdu0= dTdrho_e + dTde_rho*(-(h-P/rho) + sqvel)*(1/rho);
  su2double dTdu1= dTde_rho*(-val_Mean_PrimVar[1])*(1/rho);
  su2double dTdu2= dTde_rho*(-val_Mean_PrimVar[2])*(1/rho);

  su2double total_conductivity = val_thermal_conductivity + val_heat_capacity_cp*val_eddy_viscosity/Prandtl_Turb;
  su2double factor2 = total_conductivity/val_dist_ij;

  heat_flux_jac_i[0] = factor2*dTdu0;
  heat_flux_jac_i[1] = factor2*dTdu1;
  heat_flux_jac_i[2] = factor2*dTdu2;

  if (nDim == 2) {

    su2double dTdu3= dTde_rho*(1/rho);
    heat_flux_jac_i[3] = factor2*dTdu3;

  } else {

    su2double dTdu3= dTde_rho*(-val_Mean_PrimVar[3])*(1/rho);
    su2double dTdu4= dTde_rho*(1/rho);
    heat_flux_jac_i[3] = factor2*dTdu3;
    heat_flux_jac_i[4] = factor2*dTdu4;

  }

}
};
