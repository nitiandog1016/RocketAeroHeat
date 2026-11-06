bool CNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {

  // 直接捕获热状态计算输入，包含原生回退后的同一次调用。
  FlowPrimitiveK(iPoint) = turb_ke;

  bool RightVol = true;

  SetVelocity(iPoint); // Computes velocity and velocity^2
  su2double density      = GetDensity(iPoint);
  su2double staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint) - turb_ke;

  /*--- Check will be moved inside fluid model plus error description strings ---*/

  FluidModel->SetTDState_rhoe(density, staticEnergy);

  bool check_dens  = SetDensity(iPoint);
  bool check_press = SetPressure(iPoint, FluidModel->GetPressure());
  bool check_sos   = SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
  bool check_temp  = SetTemperature(iPoint, FluidModel->GetTemperature());

  /*--- Check that the solution has a physical meaning ---*/

  if (check_dens || check_press || check_sos  || check_temp) {

    /*--- Copy the old solution ---*/

    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);

    /*--- Recompute the primitive variables ---*/

    SetVelocity(iPoint); // Computes velocity and velocity^2
    density      = GetDensity(iPoint);
    staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint) - turb_ke;

    /*--- Check will be moved inside fluid model plus error description strings ---*/

    FluidModel->SetTDState_rhoe(density, staticEnergy);

    SetDensity(iPoint);
    SetPressure(iPoint, FluidModel->GetPressure());
    SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
    SetTemperature(iPoint, FluidModel->GetTemperature());

    RightVol = false;

  }

  /*--- Set enthalpy ---*/

  SetEnthalpy(iPoint); // Requires pressure computation.

  /*--- Set laminar viscosity ---*/

  SetLaminarViscosity(iPoint, FluidModel->GetLaminarViscosity());

  /*--- Set eddy viscosity ---*/

  SetEddyViscosity(iPoint, eddy_visc);

  /*--- Set thermal conductivity ---*/

  SetThermalConductivity(iPoint, FluidModel->GetThermalConductivity());

  /*--- Set specific heat ---*/

  SetSpecificHeatCp(iPoint, FluidModel->GetCp());

  /*--- Set look-up variables in case of data-driven fluid model ---*/
  if (DataDrivenFluid) {
    SetDataExtrapolation(iPoint, FluidModel->GetExtrapolation());
    SetEntropy(iPoint, FluidModel->GetEntropy());
  }

  return RightVol;
}
