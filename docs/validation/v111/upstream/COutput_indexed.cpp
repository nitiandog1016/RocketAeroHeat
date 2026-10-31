/*!
 * \file COutput.cpp
 * \brief Main subroutines for output solver information
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


// Original lines 425-447

    case OUTPUT_TYPE::RESTART_ASCII: case OUTPUT_TYPE::CSV:

      extension = CSU2FileWriter::fileExt;

      if (fileName.empty())
        fileName = config->GetFilename(restartFilename, "", curTimeIter);

      if (!config->GetWrt_Restart_Overwrite())
        filename_iter = config->GetFilename_Iter(fileName, curInnerIter, curOuterIter);

      LogOutputFiles("SU2 ASCII restart");

      if (config->GetWrt_Restart_Compact()) {
        /*--- If we have compact restarts, we use only the required fields. ---*/
        volumeDataSorterCompact->SetRequiredFieldNames(requiredVolumeFieldNames);
        fileWriter = new CSU2FileWriter(volumeDataSorterCompact);
      } else {
        fileWriter = new CSU2FileWriter(volumeDataSorter);
      }

      break;


// Original lines 628-650

    case OUTPUT_TYPE::SURFACE_PARAVIEW_ASCII:

      extension = CParaviewFileWriter::fileExt;

      if (fileName.empty())
        fileName = config->GetFilename(surfaceFilename, "", curTimeIter);

      if (!config->GetWrt_Surface_Overwrite())
        filename_iter = config->GetFilename_Iter(fileName, curInnerIter, curOuterIter);

      /*--- Load and sort the output data and connectivity. ---*/

      surfaceDataSorter->SortConnectivity(config, geometry);
      surfaceDataSorter->SortOutputData();

      LogOutputFiles("Paraview ASCII surface");
      fileWriter = new CParaviewFileWriter(surfaceDataSorter);

      break;

    case OUTPUT_TYPE::SURFACE_PARAVIEW_LEGACY_BINARY:


// Original lines 802-810 (base output precedes indexed copy)

  if (fileWriter != nullptr) {

    /*--- Write data to file ---*/

    fileWriter->WriteData(fileName);

    su2double BandWidth = fileWriter->GetBandwidth();


// Original lines 810-880

    /*--- Write data with iteration number to file if required ---*/

    if (!filename_iter.empty()) {
      fileWriter->WriteData(filename_iter);

      /*--- Average bandwidth ---*/
      BandWidth = (BandWidth + fileWriter->GetBandwidth()) / 2;
    }

    /*--- Compute and store the bandwidth ---*/

    if (format == OUTPUT_TYPE::RESTART_BINARY) {
      config->SetRestart_Bandwidth_Agg(config->GetRestart_Bandwidth_Agg() + BandWidth);
    }

    if (config->GetWrt_Performance() && (rank == MASTER_NODE)){
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
      (*fileWritingTable) << " " << "(" + PrintingToolbox::to_string(BandWidth) + " MB/s)";
      fileWritingTable->SetAlign(PrintingToolbox::CTablePrinter::LEFT);
    }

    delete fileWriter;

  }

}

bool COutput::GetCauchyCorrectedTimeConvergence(const CConfig *config){
  // Handle Cauchy convergence delay for 2nd order time stepping
  if(!cauchyTimeConverged && TimeConvergence && config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND){
    // Change flags for 2nd order Time stepping: In case of convergence, this iter and next iter gets written out. then solver stops
    cauchyTimeConverged = TimeConvergence;
    TimeConvergence = false;
  }
  else if(cauchyTimeConverged){
    TimeConvergence = cauchyTimeConverged;
  }

  // Handle max time delay for 2nd order time stepping
  // Delay stopping at max_time to ensure both timestep N and N-1 are written for proper restart
  if(config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND){
    const su2double cur_time = GetHistoryFieldValue("CUR_TIME");
    const su2double max_time = config->GetMax_Time();
    const bool final_time_reached = (cur_time >= max_time);

    // If max_time is reached on first detection, delay the stop
    if(final_time_reached && !maxTimeDelayActive){
      maxTimeDelayActive = true;
      TimeConvergence = false;  // Delay stop to run one more iteration
    }
    else if(maxTimeDelayActive){
      TimeConvergence = true;   // Now allow stop
      maxTimeDelayActive = false;   // Reset for next run
    }
  }

  return TimeConvergence;
}

bool COutput::SetResultFiles(CGeometry *geometry, CConfig *config, CSolver** solver_container,
                              unsigned long iter, bool force_writing) {

  bool isFileWrite = false, dataIsLoaded = false;
  const auto nVolumeFiles = config->GetnVolumeOutputFiles();
  const auto* VolumeFiles = config->GetVolumeOutputFiles();

  /*--- Check if the data sorters are allocated, if not, allocate them. --- */
  AllocateDataSorters(config, geometry);

  for (unsigned short iFile = 0; iFile < nVolumeFiles; iFile++) {

