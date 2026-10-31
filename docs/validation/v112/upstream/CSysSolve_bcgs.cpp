/*!
 * \file CSysSolve.cpp
 * \brief Main classes required for solving linear systems of equations
 * \author J. Hicken, F. Palacios, T. Economon, P. Gomes
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

// Source: https://github.com/su2code/SU2/blob/v8.6.0/Common/src/linear_algebra/CSysSolve.cpp

// Original lines 1059-1138
unsigned long CSysSolve<ScalarType>::BCGSTAB_LinSolver(const CSysVector<ScalarType>& b, CSysVector<ScalarType>& x,
                                                       const CMatrixVectorProduct<ScalarType>& mat_vec,
                                                       const CPreconditioner<ScalarType>& precond, ScalarType tol,
                                                       unsigned long m, ScalarType& residual, bool monitoring,
                                                       const CConfig* config) const {
  SU2_ZONE_SCOPED

  const bool masterRank = (SU2_MPI::GetRank() == MASTER_NODE);
  ScalarType norm_r = 0.0, norm0 = 0.0;
  unsigned long i = 0;

  /*--- Check the subspace size ---*/

  if (m < 1) {
    SU2_MPI::Error("Number of linear solver iterations must be greater than 0.", CURRENT_FUNCTION);
  }

  /*--- Allocate if not allocated yet ---*/

  if (!bcg_ready) {
    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      auto nVar = b.GetNVar();
      auto nBlk = b.GetNBlk();
      auto nBlkDomain = b.GetNBlkDomain();

      A_x.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r_0.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      r.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      p.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      v.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      z.Initialize(nBlk, nBlkDomain, nVar, nullptr);
      if (config->GetComm_Level() == COMM_FULL) x_best.Initialize(nBlk, nBlkDomain, nVar, nullptr);

      bcg_ready = true;
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

  /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/

  if (!xIsZero) {
    mat_vec(x, A_x);
    r = b - A_x;
  } else {
    r = b;
  }

  /*--- Only compute the residuals in full communication mode. ---*/

  if (config->GetComm_Level() == COMM_FULL) {
    norm0 = b.norm();
    norm_r = xIsZero ? norm0 : r.norm();

    /*--- Set the norm to the initial initial residual value ---*/

    if (tol_type == LinearToleranceType::RELATIVE) norm0 = norm_r;

    if ((norm_r < tol * norm0) || (norm_r < eps)) {
      if (masterRank) {
        SU2_OMP_MASTER
        cout << "CSysSolve::BCGSTAB(): system solved by initial guess." << endl;
        END_SU2_OMP_MASTER
      }
      return 0;
    }

    /*--- Output header information including initial residual ---*/

    if ((monitoring) && (masterRank)) {
      SU2_OMP_MASTER {
        WriteHeader("BCGSTAB", tol, norm_r);
        WriteHistory(i, norm_r / norm0);
      }
      END_SU2_OMP_MASTER
    }
  }

  /*--- Initialization ---*/

  ScalarType alpha = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;

// Original lines 1147-1245
  for (i = 0; i < m; i++) {
    /*--- Compute rho_prime ---*/

    rho_prime = rho;

    /*--- Compute rho_i ---*/

    rho = r.dot(r_0);

    /*--- Compute beta ---*/

    ScalarType beta = (rho / rho_prime) * (alpha / omega);

    /*--- Update p ---*/

    p = beta * (p - omega * v) + r;

    /*--- Preconditioning step ---*/

    precond(p, z);
    mat_vec(z, v);

    /*--- Calculate step-length alpha ---*/

    ScalarType r_0_v = r_0.dot(v);
    alpha = rho / r_0_v;

    /*--- Update solution and residual ---*/

    x += alpha * z;
    r -= alpha * v;

    /*--- Preconditioning step ---*/

    precond(r, z);
    mat_vec(z, A_x);

    /*--- Calculate step-length omega, avoid division by 0. ---*/

    omega = A_x.squaredNorm();
    if (omega == ScalarType(0)) break;
    omega = A_x.dot(r) / omega;

    /*--- Update solution and residual ---*/

    x += omega * z;
    r -= omega * A_x;

    /*--- Only compute the residuals in full communication mode. ---*/

    if (config->GetComm_Level() == COMM_FULL) {
      /*--- Check if solution has converged, else output the relative residual if necessary ---*/

      norm_r = r.norm();
      if (norm_r < norm_best) {
        norm_best = norm_r;
        x_best = x;
      }
      if (norm_r < tol * norm0) break;
      if (((monitoring) && (masterRank)) && ((i + 1) % monitorFreq == 0)) {
        SU2_OMP_MASTER
        WriteHistory(i + 1, norm_r / norm0);
        END_SU2_OMP_MASTER
      }
    }
  }

  /*--- BCGSTAB does not reduce the residual monotonically, return the best iterate. ---*/

  if ((config->GetComm_Level() == COMM_FULL) && (norm_r > norm_best)) {
    x = x_best;
    norm_r = norm_best;
  }

  /*--- Recalculate final residual (this should be optional) ---*/

  if ((monitoring) && (config->GetComm_Level() == COMM_FULL)) {
    if (masterRank) {
      SU2_OMP_MASTER
      WriteFinalResidual("BCGSTAB", i, norm_r / norm0);
      END_SU2_OMP_MASTER
    }

    if (recomputeRes) {
      mat_vec(x, A_x);
      r = b - A_x;
      ScalarType true_res = r.norm();

      if ((fabs(true_res - norm_r) > tol * 10.0) && (masterRank)) {
        SU2_OMP_MASTER
        WriteWarning(norm_r, true_res, tol);
        END_SU2_OMP_MASTER
      }
    }
  }

  residual = norm_r / norm0;
  return i;
}
