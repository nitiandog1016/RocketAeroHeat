/*!
 * \file CFVMFlowSolverBase.hpp
 * \brief Base class template for all FVM flow solvers.
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

// Source: https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/include/solvers/CFVMFlowSolverBase.hpp

// Original lines 549-570
    {
      /*--- Thread-local variables for min/max reduction. ---*/
      su2double minDt = 1e30, maxDt = 0.0;

      SU2_OMP_FOR_(schedule(static,omp_chunk_size) SU2_NOWAIT)
      for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {

        su2double Vol = geometry->nodes->GetVolume(iPoint);

        if (Vol != 0.0) {
          su2double Local_Delta_Time = nodes->GetLocalCFL(iPoint)*Vol / nodes->GetMax_Lambda_Inv(iPoint);

          if(viscous) {
            su2double dt_visc = nodes->GetLocalCFL(iPoint)*K_v*Vol*Vol / nodes->GetMax_Lambda_Visc(iPoint);
            Local_Delta_Time = min(Local_Delta_Time, dt_visc);
          }

          minDt = min(minDt, Local_Delta_Time);
          maxDt = max(maxDt, Local_Delta_Time);

          nodes->SetDelta_Time(iPoint, min(Local_Delta_Time, config->GetMax_DeltaTime()));
        }
