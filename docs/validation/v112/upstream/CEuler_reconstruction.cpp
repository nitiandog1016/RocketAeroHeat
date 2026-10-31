/*!
 * \file CEulerSolver.cpp
 * \brief Main subroutines for solving Finite-Volume Euler flow problems.
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

// Source: https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/src/solvers/CEulerSolver.cpp

// Original lines 1946-1988
        Primitive_j[iVar] = V_j[iVar] - 0.5 * lim_j * Project_Grad_j;

      }

      /*--- Recompute the reconstructed quantities in a thermodynamically consistent way. ---*/

      if (!ideal_gas || low_mach_corr) {
        ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_i, Secondary_i);
        ComputeConsistentExtrapolation(GetFluidModel(), nDim, Primitive_j, Secondary_j);
      }

      /*--- Low-Mach number correction. ---*/

      if (low_mach_corr) {
        LowMachPrimitiveCorrection(GetFluidModel(), nDim, Primitive_i, Primitive_j);
      }

      /*--- Check for non-physical solutions after reconstruction. If found, use the
       cell-average value of the solution. This is a locally 1st order approximation,
       which is typically only active during the start-up of a calculation. ---*/

      bool neg_pres_or_rho_i = (Primitive_i[prim_idx.Pressure()] < 0.0) || (Primitive_i[prim_idx.Density()] < 0.0);
      bool neg_pres_or_rho_j = (Primitive_j[prim_idx.Pressure()] < 0.0) || (Primitive_j[prim_idx.Density()] < 0.0);

      su2double R = sqrt(fabs(Primitive_j[prim_idx.Density()]/Primitive_i[prim_idx.Density()]));
      su2double sq_vel = 0.0;
      for (auto iDim = 0u; iDim < nDim; iDim++) {
        su2double RoeVelocity = (R * Primitive_j[iDim + prim_idx.Velocity()] +
                                 Primitive_i[iDim + prim_idx.Velocity()]) / (R+1);
        sq_vel += pow(RoeVelocity, 2);
      }
      su2double RoeEnthalpy = (R * Primitive_j[prim_idx.Enthalpy()] + Primitive_i[prim_idx.Enthalpy()]) / (R+1);

      const bool neg_sound_speed = ((Gamma-1)*(RoeEnthalpy-0.5*sq_vel) < 0.0);
      bool bad_recon = neg_sound_speed || neg_pres_or_rho_i || neg_pres_or_rho_j;
      bad_recon = nodes->UpdateNonPhysicalEdgeCounter(iEdge, bad_recon);
      counter_local += bad_recon;

      numerics->SetPrimitive(bad_recon? V_i : Primitive_i,  bad_recon? V_j : Primitive_j);
      numerics->SetSecondary(bad_recon? S_i : Secondary_i,  bad_recon? S_j : Secondary_j);

    }
