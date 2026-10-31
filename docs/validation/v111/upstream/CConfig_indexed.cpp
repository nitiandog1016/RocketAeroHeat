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


// Original lines 9028-9034
string CConfig::GetFilename_Iter(const string& filename_iter, unsigned long curInnerIter, unsigned long curOuterIter) const {
  const auto iter = GetMultizone_Problem() ? curOuterIter : curInnerIter;
  std::stringstream iter_ss;
  iter_ss << filename_iter << "_" << std::setw(6) << std::setfill('0') << iter;
  return iter_ss.str();
}


