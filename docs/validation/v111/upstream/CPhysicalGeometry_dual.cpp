/*!
 * \file CPhysicalGeometry.cpp
 * \brief Implementation of the physical geometry class.
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


// Original lines 6748-6888
void CPhysicalGeometry::SetControlVolume(CConfig* config, unsigned short action) {
  /*--- Update values of faces of the edge ---*/
  if (action != ALLOCATE) {
    su2double ZeroArea[MAXNDIM] = {0.0};

    SU2_OMP_FOR_STAT(1024)
    for (auto iEdge = 0ul; iEdge < nEdge; iEdge++) edges->SetNormal(iEdge, ZeroArea);
    END_SU2_OMP_FOR

    SU2_OMP_FOR_STAT(1024)
    for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) nodes->SetVolume(iPoint, 0.0);
    END_SU2_OMP_FOR
  }

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS { /*--- The following is difficult to parallelize with threads. ---*/

    su2double my_DomainVolume = 0.0;
    for (auto iElem = 0ul; iElem < nElem; iElem++) {
      const auto nNodes = elem[iElem]->GetnNodes();

      /*--- To make preaccumulation more effective, use as few inputs
       as possible, recomputing intermediate quantities as needed. ---*/
      AD::StartPreacc();

      /*--- Get pointers to the coordinates of all the element nodes ---*/
      array<const su2double*, N_POINTS_MAXIMUM> Coord;

      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
        auto iPoint = elem[iElem]->GetNode(iNode);
        Coord[iNode] = nodes->GetCoord(iPoint);
#ifdef CODI_REVERSE_TYPE
        /*--- The same points and edges will be referenced multiple times as they are common
         to many of the element's faces, therefore they are "registered" here only once. ---*/
        AD::SetPreaccIn(nodes->Volume(iPoint));
        for (unsigned short jNode = iNode + 1; jNode < nNodes; jNode++) {
          auto jPoint = elem[iElem]->GetNode(jNode);
          auto iEdge = FindEdge(iPoint, jPoint, false);
          if (iEdge >= 0) AD::SetPreaccIn(edges->Normal[iEdge], nDim);
        }
#endif
      }
      AD::SetPreaccIn(Coord, nNodes, nDim);

      /*--- Compute the element median CG coordinates ---*/
      auto Coord_Elem_CG = elem[iElem]->SetCoord_CG(nDim, Coord);
      AD::SetPreaccOut(Coord_Elem_CG, nDim);

      for (unsigned short iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {
        /*--- In 2D all the faces have only one edge ---*/
        unsigned short nEdgesFace = 1;

        /*--- In 3D the number of edges per face is the same as the number of point
         per face and the median CG of the face is needed. ---*/
        su2double Coord_FaceElem_CG[MAXNDIM] = {0.0};
        if (nDim == 3) {
          nEdgesFace = elem[iElem]->GetnNodesFace(iFace);

          for (unsigned short iNode = 0; iNode < nEdgesFace; iNode++) {
            auto NodeFace = elem[iElem]->GetFaces(iFace, iNode);
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
              Coord_FaceElem_CG[iDim] += Coord[NodeFace][iDim] / nEdgesFace;
          }
        }

        /*-- Loop over the edges of a face ---*/
        for (unsigned short iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {
          const auto face_iNode = elem[iElem]->GetFaces(iFace, iEdgesFace);
          unsigned short face_jNode;

          if (nDim == 2) {
            /*--- In 2D only one edge (two points) per edge ---*/
            face_jNode = elem[iElem]->GetFaces(iFace, 1);
          } else {
            /*--- In 3D we "circle around" the face ---*/
            face_jNode = elem[iElem]->GetFaces(iFace, (iEdgesFace + 1) % nEdgesFace);
          }

          const auto face_iPoint = elem[iElem]->GetNode(face_iNode);
          const auto face_jPoint = elem[iElem]->GetNode(face_jNode);

          /*--- We define a direction (from the smalest index to the greatest) --*/
          const bool change_face_orientation = (face_iPoint > face_jPoint);
          const auto iEdge = FindEdge(face_iPoint, face_jPoint);

          su2double Coord_Edge_CG[MAXNDIM] = {0.0};
          for (unsigned short iDim = 0; iDim < nDim; iDim++) {
            Coord_Edge_CG[iDim] = 0.5 * (Coord[face_iNode][iDim] + Coord[face_jNode][iDim]);
          }

          su2double Volume_i, Volume_j;

          if (nDim == 2) {
            /*--- Two dimensional problem ---*/
            if (change_face_orientation)
              edges->SetNodes_Coord(iEdge, Coord_Elem_CG, Coord_Edge_CG);
            else
              edges->SetNodes_Coord(iEdge, Coord_Edge_CG, Coord_Elem_CG);

            Volume_i = CEdge::GetVolume(Coord[face_iNode], Coord_Edge_CG, Coord_Elem_CG);
            Volume_j = CEdge::GetVolume(Coord[face_jNode], Coord_Edge_CG, Coord_Elem_CG);
          } else {
            /*--- Three dimensional problem ---*/
            if (change_face_orientation)
              edges->SetNodes_Coord(iEdge, Coord_FaceElem_CG, Coord_Edge_CG, Coord_Elem_CG);
            else
              edges->SetNodes_Coord(iEdge, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);

            Volume_i = CEdge::GetVolume(Coord[face_iNode], Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
            Volume_j = CEdge::GetVolume(Coord[face_jNode], Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
          }

          nodes->AddVolume(face_iPoint, Volume_i);
          nodes->AddVolume(face_jPoint, Volume_j);

          my_DomainVolume += Volume_i + Volume_j;
        }
      }

#ifdef CODI_REVERSE_TYPE
      for (unsigned short iNode = 0; iNode < nNodes; iNode++) {
        auto iPoint = elem[iElem]->GetNode(iNode);
        AD::SetPreaccOut(nodes->Volume(iPoint));
        for (unsigned short jNode = iNode + 1; jNode < nNodes; jNode++) {
          auto jPoint = elem[iElem]->GetNode(jNode);
          auto iEdge = FindEdge(iPoint, jPoint, false);
          if (iEdge >= 0) AD::SetPreaccOut(edges->Normal[iEdge], nDim);
        }
      }
#endif
      AD::EndPreacc();
    }

    su2double DomainVolume;
    SU2_MPI::Allreduce(&my_DomainVolume, &DomainVolume, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
    config->SetDomainVolume(DomainVolume);

    if ((rank == MASTER_NODE) && (action == ALLOCATE)) {
      if (nDim == 2) cout << "Area of the computational grid: " << DomainVolume << "." << endl;
      if (nDim == 3) cout << "Volume of the computational grid: " << DomainVolume << "." << endl;
    }
  }

