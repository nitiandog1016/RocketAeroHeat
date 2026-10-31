/*!
 * \file CFVMFlowSolverBase.inl
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

// Source: https://github.com/su2code/SU2/blob/v8.6.0/SU2_CFD/include/solvers/CFVMFlowSolverBase.inl

// Original lines 587-643
template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::CompleteImplicitIteration(CGeometry *geometry, CSolver**, CConfig *config) {
  SU2_ZONE_SCOPED

  if constexpr (R == ENUM_REGIME::COMPRESSIBLE) ComputeUnderRelaxationFactor(config);

  /*--- Update solution with under-relaxation and communicate it. ---*/

  if (!config->GetContinuous_Adjoint()) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        nodes->AddSolution(iPoint, iVar, nodes->GetUnderRelaxation(iPoint) * LinSysSol(iPoint, iVar));
      }
    }
    END_SU2_OMP_FOR
  }

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  /*--- For verification cases, compute the global error metrics. ---*/
  ComputeVerificationError(geometry, config);
}

template <class V, ENUM_REGIME R>
void CFVMFlowSolverBase<V, R>::ImplicitEuler_Iteration(CGeometry *geometry, CSolver**, CConfig *config) {
  SU2_ZONE_SCOPED

  PrepareImplicitIteration(geometry, nullptr, config);

  /*--- Solve or smooth the linear system. ---*/

  SU2_OMP_FOR_(schedule(static,OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned long iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }
  END_SU2_OMP_FOR

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    SetIterLinSolver(iter);
    SetResLinSolver(System.GetResidual());
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  CompleteImplicitIteration(geometry, nullptr, config);
}

template <class V, ENUM_REGIME R>
