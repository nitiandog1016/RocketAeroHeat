// 此函数由官方v8.6.0原实现最小扩展；仅接入SI温度表及真实输入回执。
void CGeometry::SetCustomBoundary(CConfig* config) {
  try {
    RAHAxialWallTemperature rah;
    std::ofstream applied;
    unsigned long written = 0;
    unsigned short thermalMarkers = 0;
    if (rah.Active()) {
      if (size != 1 || config->GetnZone() != 1 || config->GetnMGLevels() != 0 ||
          config->GetGrid_Movement() || config->GetTime_Domain() || config->GetSystemMeasurements() != SI)
        throw std::runtime_error("RAH spatial wall requires serial-MPI, one zone, fixed SI mesh and MGLEVEL=0");
      applied.open("rah_wall_temperature_applied.csv", std::ios::out | std::ios::trunc);
      applied.imbue(std::locale::classic());
      if (!applied) throw std::runtime_error("Cannot write RAH wall-input receipt");
      applied << "globalNodeID,X_m,Y_m,Z_m,givenTemperature_K\n" << std::setprecision(17);
    }
    CustomBoundaryTemperature = new su2double*[nMarker];
    CustomBoundaryHeatFlux = new su2double*[nMarker];
    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
      const std::string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      CustomBoundaryHeatFlux[iMarker] = nullptr;
      CustomBoundaryTemperature[iMarker] = nullptr;
      if (!config->GetMarker_All_PyCustom(iMarker)) continue;
      if (rah.Active() && (Marker_Tag != "wall" || config->GetMarker_All_KindBC(iMarker) != ISOTHERMAL))
        throw std::runtime_error("RAH spatial wall requires one custom ISOTHERMAL marker named wall");
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case HEAT_FLUX:
          CustomBoundaryHeatFlux[iMarker] = new su2double[nVertex[iMarker]];
          for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++)
            CustomBoundaryHeatFlux[iMarker][iVertex] = config->GetWall_HeatFlux(Marker_Tag);
          break;
        case ISOTHERMAL:
          ++thermalMarkers;
          CustomBoundaryTemperature[iMarker] = new su2double[nVertex[iMarker]];
          for (auto iVertex = 0ul; iVertex < nVertex[iMarker]; iVertex++) {
            CustomBoundaryTemperature[iMarker][iVertex] = config->GetIsothermal_Temperature(Marker_Tag);
            if (rah.Active()) {
              const auto iPoint = vertex[iMarker][iVertex]->GetNode();
              const double x = SU2_TYPE::GetValue(nodes->GetCoord(iPoint,0));
              const double y = SU2_TYPE::GetValue(nodes->GetCoord(iPoint,1));
              const double z = nDim==3 ? SU2_TYPE::GetValue(nodes->GetCoord(iPoint,2)) : 0.0;
              const double value = rah.At(x);
              CustomBoundaryTemperature[iMarker][iVertex] = value;
              applied << nodes->GetGlobalIndex(iPoint) << ',' << x << ',' << y << ',' << z << ',' << value << '\n';
              ++written;
            }
          }
          break;
        case INLET_FLOW:
          break;
        default:
          std::cout << "WARNING: Marker " << Marker_Tag << " is not customizable. Using default behavior." << std::endl;
          break;
      }
    }
    if (rah.Active()) {
      if (thermalMarkers != 1 || written == 0)
        throw std::runtime_error("RAH temperature table was not applied to exactly one nonempty wall");
      applied.close();
      if (!applied) throw std::runtime_error("RAH wall-input receipt write failed");
      std::cout << "RAH_AXIAL_WALL_V1 applied " << written << " wall nodes; given SI temperatures, not solved material T." << std::endl;
    }
  } catch (const std::exception& ex) {
    SU2_MPI::Error(ex.what(), CURRENT_FUNCTION);
  }
}
