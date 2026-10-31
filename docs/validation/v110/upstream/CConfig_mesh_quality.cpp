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


// Original CConfig.cpp lines 3808-3815.
  /*--- Check if MESH_QUALITY is requested in VOLUME_OUTPUT and set the config boolean accordingly. ---*/
  Wrt_MeshQuality = false;
  for (unsigned short iField = 0; iField < nVolumeOutput; iField++) {
    if(VolumeOutput[iField].find("MESH_QUALITY") != string::npos) {
      Wrt_MeshQuality = true;
    }
  }


