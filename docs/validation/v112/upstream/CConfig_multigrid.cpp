/*!
 * \file CConfig.cpp
 * \brief Main file for managing the config file
 * \author F. Palacios, T. Economon, B. Tracey, H. Kline
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

// Source: https://github.com/su2code/SU2/blob/v8.6.0/Common/src/CConfig.cpp

// Original lines 2099-2133

  /*!\par CONFIG_CATEGORY: Multi-grid \ingroup Config*/
  /*!\brief MGLEVEL\n DESCRIPTION: Multi-grid Levels. DEFAULT: 0 \ingroup Config*/
  addUnsignedShortOption("MGLEVEL", nMGLevels, 0);
  /*!\brief MGCYCLE\n DESCRIPTION: Multi-grid cycle. OPTIONS: See \link MG_Cycle_Map \endlink. Defualt V_CYCLE \ingroup Config*/
  addEnumOption("MGCYCLE", Kind_MGCycle, MG_Cycle_Map, MG_CYCLE::V);
  /*!\brief MG_PRE_SMOOTH\n DESCRIPTION: Multi-grid pre-smoothing level \ingroup Config*/
  addUShortListOption("MG_PRE_SMOOTH", nMG_PreSmooth_p, MG_PreSmooth_p);
  /*!\brief MG_POST_SMOOTH\n DESCRIPTION: Multi-grid post-smoothing level \ingroup Config*/
  addUShortListOption("MG_POST_SMOOTH", nMG_PostSmooth_p, MG_PostSmooth_p);
  /*!\brief MG_CORRECTION_SMOOTH\n DESCRIPTION: Jacobi implicit smoothing of the correction \ingroup Config*/
  addUShortListOption("MG_CORRECTION_SMOOTH", nMG_CorrecSmooth_p, MG_CorrecSmooth_p);
  /*!\brief MG_DAMP_RESTRICTION\n DESCRIPTION: Damping factor for the residual restriction. DEFAULT: 0.75 \ingroup Config*/
  addDoubleOption("MG_DAMP_RESTRICTION", Damp_Res_Restric, 0.5);
  /*!\brief MG_DAMP_PROLONGATION\n DESCRIPTION: Damping factor for the correction prolongation. DEFAULT 0.75 \ingroup Config*/
  addDoubleOption("MG_DAMP_PROLONGATION", Damp_Correc_Prolong, 0.5);
  /*!\brief MG_SMOOTH_EARLY_EXIT\n DESCRIPTION: Enable early exit for MG smoothing when RMS drops below threshold. DEFAULT: NO \ingroup Config*/
  addBoolOption("MG_SMOOTH_EARLY_EXIT", MGOptions.MG_Smooth_EarlyExit, true);
  /*!\brief MG_SMOOTH_RES_THRESHOLD\n DESCRIPTION: Early exit smoothing when current_rms drops below threshold * initial_rms. DEFAULT: 0.9 \ingroup Config*/
  addDoubleOption("MG_SMOOTH_RES_THRESHOLD", MGOptions.MG_Smooth_Res_Threshold, 0.9);
  /*!\brief MG_SMOOTH_OUTPUT\n DESCRIPTION: Print compact per-cycle smoothing iteration summary. DEFAULT: NO \ingroup Config*/
  addBoolOption("MG_SMOOTH_OUTPUT", MGOptions.MG_Smooth_Output, false);
  /*!\brief MG_SMOOTH_STAGNATION_TOL\n DESCRIPTION: Stop smoothing if current_rms >= previous_rms * this value.
   * Values < 1.0 enable early exit on stagnation, 1.0 only exits on defect growth. DEFAULT: 0.99 \ingroup Config*/
  addDoubleOption("MG_SMOOTH_STAGNATION_TOL", MGOptions.MG_Smooth_StagnationTol, 0.99);
  /*!\brief MG_SMOOTH_COEFF\n DESCRIPTION: Smoothing coefficient for the correction prolongation Jacobi smoother. DEFAULT: 1.25 \ingroup Config*/
  addDoubleOption("MG_SMOOTH_COEFF", MGOptions.MG_Smooth_Coeff, 1.25);
  /*!\brief MG_MIN_MESHSIZE\n DESCRIPTION: Minimum number of CVs on the coarsest multigrid level, counted over the
   * whole level rather than per rank, so the hierarchy does not change with the number of ranks. DEFAULT: 500 \ingroup Config*/
  addUnsignedLongOption("MG_MIN_MESHSIZE", MGOptions.MG_Min_MeshSize, 500);
  /*!\brief MG_IMPLICIT_LINES\n DESCRIPTION: Pave the coarse grid with advancing fronts raised from boundaries
   * that carry a stretched layer normal to themselves. DEFAULT: NO \ingroup Config*/
  addBoolOption("MG_IMPLICIT_LINES", MGOptions.MG_Implicit_Lines, false);
  /*!\brief MG_STARTUP_ITER\n DESCRIPTION: Max number of iterations spent on each mesh during the Full
   * Multigrid (FMG) startup phase. DEFAULT: 100 \ingroup Config*/
