/*!
 * \file CNSVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
 * \version 8.6.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/variables/CNSVariable.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../include/fluid/CFluidModel.hpp"
#include "../../include/fluid/CSU2IdealAir.hpp"
#include "../../include/variables/ConservativeStepTrace.hpp"
#include <array>
#include <string>

CNSVariable::CNSVariable(su2double density, const su2double *velocity, su2double energy,
                         unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config) :
                         CEulerVariable(density,velocity,energy,npoint,ndim,nvar,config) {

  inv_TimeScale = config->GetModVel_FreeStream() / config->GetRefLength();

  Vorticity.resize(nPoint,3) = su2double(0.0);
  StrainMag.resize(nPoint) = su2double(0.0);
  FlowPrimitiveK.resize(nPoint) = su2double(0.0);
  Tau_Wall.resize(nPoint) = su2double(-1.0);
  DES_LengthScale.resize(nPoint) = su2double(0.0);
  lesMode.resize(nPoint) = su2double(0.0);
  Roe_Dissipation.resize(nPoint) = su2double(0.0);
  Vortex_Tilting.resize(nPoint) = su2double(0.0);
  Max_Lambda_Visc.resize(nPoint) = su2double(0.0);

}

void CNSVariable::SetRoe_Dissipation_NTS(unsigned long iPoint,
                                         su2double val_delta,
                                         su2double val_const_DES){

  const su2double cnu = pow(0.09, 1.5),
                  ch1 = 3.0,
                  ch2 = 1.0,
                  ch3 = 2.0,
                  sigma_max = 1.0;

  unsigned long iDim;
  su2double Omega, Omega_2 = 0.0, Baux, Gaux, Lturb, Kaux, Aaux;

  AD::StartPreacc();
  AD::SetPreaccIn(Vorticity[iPoint], 3);
  AD::SetPreaccIn(StrainMag(iPoint));
  AD::SetPreaccIn(val_delta);
  AD::SetPreaccIn(val_const_DES);
  /*--- Density ---*/
  AD::SetPreaccIn(Solution(iPoint,0));
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(Primitive(iPoint, indices.LaminarViscosity()));
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(Primitive(iPoint, indices.EddyViscosity()));

  /*--- Central/upwind blending based on:
   * Zhixiang Xiao, Jian Liu, Jingbo Huang, and Song Fu.  "Numerical
   * Dissipation Effects on Massive Separation Around Tandem Cylinders",
   * AIAA Journal, Vol. 50, No. 5 (2012), pp. 1119-1136.
   * https://doi.org/10.2514/1.J051299
   * ---*/

  for (iDim = 0; iDim < 3; iDim++){
    Omega_2 += pow(Vorticity(iPoint,iDim),2);
  }
  Omega = sqrt(Omega_2);

  Baux = (ch3 * Omega * max(StrainMag(iPoint), Omega)) /
      max((pow(StrainMag(iPoint),2)+Omega_2)*0.5, 1E-20);
  Gaux = tanh(pow(Baux,4.0));

  Kaux = max(sqrt((Omega_2 + pow(StrainMag(iPoint), 2))*0.5), 0.1 * inv_TimeScale);

  const su2double nu = GetLaminarViscosity(iPoint)/GetDensity(iPoint);
  const su2double nu_t = GetEddyViscosity(iPoint)/GetDensity(iPoint);
  Lturb = sqrt((nu + nu_t)/(cnu*Kaux));

  Aaux = ch2*max((val_const_DES*val_delta/Lturb)/Gaux -  0.5, 0.0);

  Roe_Dissipation(iPoint) = sigma_max * tanh(pow(Aaux, ch1));

  AD::SetPreaccOut(Roe_Dissipation(iPoint));
  AD::EndPreacc();

}

void CNSVariable::SetRoe_Dissipation_FD(unsigned long iPoint, su2double val_wall_dist){

  /*--- Constants for Roe Dissipation ---*/

  const passivedouble k2 = pow(0.41,2.0);

  AD::StartPreacc();
  AD::SetPreaccIn(Gradient_Primitive[iPoint], nVar, nDim);
  AD::SetPreaccIn(val_wall_dist);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(Primitive(iPoint, indices.EddyViscosity()));
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(Primitive(iPoint, indices.LaminarViscosity()));
  /*--- Density; GetDensity reads from Solution (not Primitive) at index 0 ---*/
  AD::SetPreaccIn(Solution(iPoint, 0));

  su2double uijuij = 0.0;

  for(unsigned long iDim = 0; iDim < nDim; ++iDim)
    for(unsigned long jDim = 0; jDim < nDim; ++jDim)
      uijuij += pow(Gradient_Primitive(iPoint, indices.Velocity()+iDim, jDim),2);

  uijuij = max(sqrt(uijuij),1e-10);

  const su2double nu = GetLaminarViscosity(iPoint)/GetDensity(iPoint);
  const su2double nu_t = GetEddyViscosity(iPoint)/GetDensity(iPoint);
  const su2double r_d = (nu + nu_t)/(uijuij*k2*pow(val_wall_dist,2));

  Roe_Dissipation(iPoint) = 1.0-tanh(pow(8.0*r_d,3.0));

  AD::SetPreaccOut(Roe_Dissipation(iPoint));
  AD::EndPreacc();
}

bool CNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {
  // 模块一。使用本次实际SST k计算热状态；合法候选不改任何守恒量。
  FlowPrimitiveK(iPoint)=turb_ke;
  bool RightVol=true;
  su2double density=0, staticEnergy=0;
  auto query=[&]() {
    SetVelocity(iPoint);
    density=GetDensity(iPoint);
    staticEnergy=GetEnergy(iPoint)-0.5*Velocity2(iPoint)-turb_ke;
    FluidModel->SetTDState_rhoe(density,staticEnergy);
  };
  try {
    query();
  } catch(const std::domain_error& error) {
    if(dynamic_cast<CSU2IdealAir*>(FluidModel)==nullptr) throw;
    const std::string reason=error.what();
    const double trial_energy=staticEnergy;
    rah_step::State trial{},old{},accepted{};
    for(unsigned long j=0;j<nVar;++j) {
      const unsigned long column=(j==nVar-1)?4:j;
      trial[column]=Solution(iPoint,j);
      old[column]=Solution_Old(iPoint,j);
      Solution(iPoint,j)=Solution_Old(iPoint,j);
    }
    // 模块二。旧态也必须在当前k下合法；不得裁温、清k或使用陈旧缓存。
    try {
      query();
    } catch(const std::domain_error& old_error) {
      rah_step::record(iPoint,nDim,trial,old,old,trial_energy,staticEnergy,staticEnergy,
                       turb_ke,0,0,false,old_error.what());
      throw;
    }
    const double old_energy=staticEnergy;
    double theta=1.0;
    unsigned attempts=0;
    bool found=false;
    // 模块三。整组守恒分量共用一个有限半步，真实EOS同时审核上下温界。
    for(unsigned attempt=1;attempt<=32;++attempt) {
      theta*=0.5;
      attempts=attempt;
      for(unsigned long j=0;j<nVar;++j) {
        const unsigned long column=(j==nVar-1)?4:j;
        Solution(iPoint,j)=old[column]+theta*(trial[column]-old[column]);
      }
      try {
        query();
        found=true;
        break;
      } catch(const std::domain_error&) {}
    }
    if(!found) {
      theta=0;
      for(unsigned long j=0;j<nVar;++j) Solution(iPoint,j)=Solution_Old(iPoint,j);
      query();
    }
    for(unsigned long j=0;j<nVar;++j) accepted[(j==nVar-1)?4:j]=Solution(iPoint,j);
    rah_step::record(iPoint,nDim,trial,old,accepted,trial_energy,old_energy,staticEnergy,
                     turb_ke,theta,attempts,true,reason.c_str());
    RightVol=false;
  }

  // 模块四。保留其他模型原有的非物理状态检查与恢复路径。
  bool check_dens=SetDensity(iPoint);
  bool check_press=SetPressure(iPoint,FluidModel->GetPressure());
  bool check_sos=SetSoundSpeed(iPoint,FluidModel->GetSoundSpeed2());
  bool check_temp=SetTemperature(iPoint,FluidModel->GetTemperature());
  if(check_dens || check_press || check_sos || check_temp) {
    for(unsigned long j=0;j<nVar;++j) Solution(iPoint,j)=Solution_Old(iPoint,j);
    query();
    SetDensity(iPoint);
    SetPressure(iPoint,FluidModel->GetPressure());
    SetSoundSpeed(iPoint,FluidModel->GetSoundSpeed2());
    SetTemperature(iPoint,FluidModel->GetTemperature());
    RightVol=false;
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

void CNSVariable::SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) {
  if (nSecondaryVar == 0) return;

  /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/

  SetdPdrho_e( iPoint, FluidModel->GetdPdrho_e() );
  SetdPde_rho( iPoint, FluidModel->GetdPde_rho() );

  SetdTdrho_e( iPoint, FluidModel->GetdTdrho_e() );
  SetdTde_rho( iPoint, FluidModel->GetdTde_rho() );

  /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

  Setdmudrho_T( iPoint, FluidModel->Getdmudrho_T() );
  SetdmudT_rho( iPoint, FluidModel->GetdmudT_rho() );

  Setdktdrho_T( iPoint, FluidModel->Getdktdrho_T() );
  SetdktdT_rho( iPoint, FluidModel->GetdktdT_rho() );

}
