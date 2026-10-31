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
// SU2 v8.6.0 upstream excerpt; original lines 3995-4278
void CPhysicalGeometry::Check_IntElem_Orientation(const CConfig* config) {
  unsigned long tria_flip = 0, quad_flip = 0, tet_flip = 0, prism_flip = 0, hexa_flip = 0, pyram_flip = 0;
  unsigned long quad_error = 0, prism_error = 0, hexa_error = 0, pyram_error = 0;

  SU2_OMP_PARALLEL_(reduction(+ : tria_flip, quad_flip, tet_flip, prism_flip, hexa_flip, pyram_flip)) {
    /*--- Lambda to test triangles. Normal should be positive in the z direction (right hand rule). ---*/
    auto checkTria = [this](unsigned long iElem, int Node_1, int Node_2, int Node_3) {
      const auto Coord_1 = nodes->GetCoord(elem[iElem]->GetNode(Node_1));
      const auto Coord_2 = nodes->GetCoord(elem[iElem]->GetNode(Node_2));
      const auto Coord_3 = nodes->GetCoord(elem[iElem]->GetNode(Node_3));
      constexpr int nDim = 2;
      su2double a[nDim] = {0.0}, b[nDim] = {0.0};
      GeometryToolbox::Distance(nDim, Coord_2, Coord_1, a);
      GeometryToolbox::Distance(nDim, Coord_3, Coord_1, b);
      return a[0] * b[1] - a[1] * b[0] < 0.0;
    };

    /*--- Lambda to test tetrahedrons, volume must be positive,
     * the normal of any face must point towards the other point. ---*/
    auto checkTetra = [this](unsigned long iElem, int Node_1, int Node_2, int Node_3, int Node_4) {
      const auto Coord_1 = nodes->GetCoord(elem[iElem]->GetNode(Node_1));
      const auto Coord_2 = nodes->GetCoord(elem[iElem]->GetNode(Node_2));
      const auto Coord_3 = nodes->GetCoord(elem[iElem]->GetNode(Node_3));
      const auto Coord_4 = nodes->GetCoord(elem[iElem]->GetNode(Node_4));
      constexpr int nDim = 3;
      su2double a[nDim] = {0.0}, b[nDim] = {0.0}, c[nDim] = {0.0}, n[nDim] = {0.0};
      GeometryToolbox::Distance(nDim, Coord_2, Coord_1, a);
      GeometryToolbox::Distance(nDim, Coord_3, Coord_1, b);
      GeometryToolbox::Distance(nDim, Coord_4, Coord_1, c);
      GeometryToolbox::CrossProduct(a, b, n);
      return GeometryToolbox::DotProduct(nDim, n, c) < 0.0;
    };

    /*--- Loop over all the elements. ---*/

    SU2_OMP_FOR_DYN(roundUpDiv(nElem, 2 * omp_get_max_threads()))
    for (auto iElem = 0ul; iElem < nElem; iElem++) {
      /*--- 2D grid. ---*/

      if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
        if (checkTria(iElem, 0, 1, 2)) {
          elem[iElem]->Change_Orientation();
          tria_flip++;
        }
      }

      if (elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
        /*--- Two triangles. ---*/
        bool test_1 = checkTria(iElem, 0, 1, 2);
        bool test_2 = checkTria(iElem, 0, 2, 3);

        if (test_1 && test_2) {
          elem[iElem]->Change_Orientation();
          quad_flip++;
        } else if (test_1 || test_2) {
          /*--- If one test fails and the other passes the
           * element probably has serious problems. ---*/
          SU2_OMP_ATOMIC
          quad_error++;
        }
      }

      /*--- 3D grid. ---*/

      if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
        if (checkTetra(iElem, 0, 1, 2, 3)) {
          elem[iElem]->Change_Orientation();
          tet_flip++;
        }
      }

      if (elem[iElem]->GetVTK_Type() == PYRAMID) {
        /*--- Slice across top vertex into 2 tets. ---*/
        bool test_1 = checkTetra(iElem, 0, 1, 2, 4);
        bool test_2 = checkTetra(iElem, 2, 3, 0, 4);

        if (test_1 && test_2) {
          elem[iElem]->Change_Orientation();
          pyram_flip++;
        } else if (test_1 || test_2) {
          SU2_OMP_ATOMIC
          pyram_error++;
        }
      }

      if (elem[iElem]->GetVTK_Type() == PRISM) {
        /*--- The triangular faces should point at each other. ---*/
        bool test_1 = checkTetra(iElem, 0, 2, 1, 3);
        bool test_2 = checkTetra(iElem, 3, 4, 5, 2);

        if (test_1 && test_2) {
          elem[iElem]->Change_Orientation();
          prism_flip++;
        } else if (test_1 || test_2) {
          SU2_OMP_ATOMIC
          prism_error++;
        }
      }

      if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
        /*--- The base points at the top. ---*/
        bool test_1 = checkTetra(iElem, 0, 1, 2, 5);
        bool test_2 = checkTetra(iElem, 0, 2, 3, 7);
        /*--- The top points at the base. ---*/
        bool test_3 = checkTetra(iElem, 4, 6, 5, 1);
        bool test_4 = checkTetra(iElem, 4, 7, 6, 3);

        if (test_1 && test_2 && test_3 && test_4) {
          elem[iElem]->Change_Orientation();
          hexa_flip++;
        } else if (test_1 || test_2 || test_3 || test_4) {
          SU2_OMP_ATOMIC
          hexa_error++;
        }
      }
    }
    END_SU2_OMP_FOR
  }
  END_SU2_OMP_PARALLEL

  auto reduce = [](unsigned long& val) {
    unsigned long tmp = val;
    SU2_MPI::Allreduce(&tmp, &val, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  };
  reduce(tria_flip);
  reduce(quad_flip);
  reduce(tet_flip);
  reduce(pyram_flip);
  reduce(prism_flip);
  reduce(hexa_flip);
  reduce(quad_error);
  reduce(pyram_error);
  reduce(prism_error);
  reduce(hexa_error);

  if (rank == MASTER_NODE) {
    string start("There has been a re-orientation of ");
    if (tria_flip) cout << start << tria_flip << " TRIANGLE volume elements." << endl;
    if (quad_flip) cout << start << quad_flip << " QUADRILATERAL volume elements." << endl;
    if (tet_flip) cout << start << tet_flip << " TETRAHEDRON volume elements." << endl;
    if (hexa_flip) cout << start << hexa_flip << " HEXAHEDRON volume elements." << endl;
    if (pyram_flip) cout << start << pyram_flip << " PYRAMID volume elements." << endl;
    if (prism_flip) cout << start << prism_flip << " PRISM volume elements." << endl;

    if (quad_error + pyram_error + prism_error + hexa_error) {
      cout << ">>> WARNING: ";
      if (quad_error) cout << quad_error << " QUADRILATERAL, ";
      if (pyram_error) cout << pyram_error << " PYRAMID, ";
      if (prism_error) cout << prism_error << " PRISM, ";
      if (hexa_error) cout << hexa_error << " HEXAHEDRON, ";
      cout << "volume elements are distorted.\n    It was not possible "
              "to determine if their orientation is correct."
           << endl;
    }

    if (tria_flip + quad_flip + tet_flip + hexa_flip + pyram_flip + prism_flip + quad_error + pyram_error +
            prism_error + hexa_error ==
        0) {
      cout << "All volume elements are correctly oriented." << endl;
    }
  }
}

void CPhysicalGeometry::Check_BoundElem_Orientation(const CConfig* config) {
  unsigned long line_flip = 0, tria_flip = 0, quad_flip = 0, quad_error = 0;

  SU2_OMP_PARALLEL_(reduction(+ : line_flip, tria_flip, quad_flip, quad_error)) {
    /*--- Lambda to test tetrahedrons. ---*/
    auto checkTetra = [this](unsigned long Point_1, unsigned long Point_2, unsigned long Point_3,
                             unsigned long Point_4) {
      const auto Coord_1 = nodes->GetCoord(Point_1);
      const auto Coord_2 = nodes->GetCoord(Point_2);
      const auto Coord_3 = nodes->GetCoord(Point_3);
      const auto Coord_4 = nodes->GetCoord(Point_4);
      constexpr int nDim = 3;
      su2double a[nDim] = {0.0}, b[nDim] = {0.0}, c[nDim] = {0.0}, n[nDim] = {0.0};
      GeometryToolbox::Distance(nDim, Coord_2, Coord_1, a);
      GeometryToolbox::Distance(nDim, Coord_3, Coord_1, b);
      GeometryToolbox::Distance(nDim, Coord_4, Coord_1, c);
      GeometryToolbox::CrossProduct(a, b, n);
      return GeometryToolbox::DotProduct(nDim, n, c) < 0.0;
    };

    for (auto iMarker = 0u; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == INTERNAL_BOUNDARY) continue;

      SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
      for (auto iElem_Surface = 0ul; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
        /*--- Pick a reference point inside the domain that is not part of the surface element. ---*/
        const auto iElem_Domain = bound[iMarker][iElem_Surface]->GetDomainElement();
        unsigned long Point_Domain = 0;

        for (auto iNode_Domain = 0u; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
          Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
          bool find = false;
          for (auto iNode_Surface = 0u; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
            auto Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
            if (Point_Surface == Point_Domain) {
              find = true;
              break;
            }
          }
          if (!find) break;
        }

        /*--- 2D grid. ---*/

        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == LINE) {
          auto Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0);
          auto Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1);
          const auto Coord_1 = nodes->GetCoord(Point_1_Surface);
          const auto Coord_2 = nodes->GetCoord(Point_2_Surface);
          const auto Coord_3 = nodes->GetCoord(Point_Domain);

          /*--- The normal of the triangle formed by the line and domain
           * point should point in the positive z direction. ---*/
          constexpr int nDim = 2;
          su2double a[nDim] = {0.0}, b[nDim] = {0.0};
          GeometryToolbox::Distance(nDim, Coord_2, Coord_1, a);
          GeometryToolbox::Distance(nDim, Coord_3, Coord_1, b);
          bool test = a[0] * b[1] - a[1] * b[0] < 0.0;

          if (test) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            line_flip++;
          }
        }

        /*--- 3D grid. ---*/

        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TRIANGLE) {
          auto Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0);
          auto Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1);
          auto Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2);

          /*--- The normal of the triangle should point into the domain,
           * resulting in a tetrahedron with positive volume. ---*/
          if (checkTetra(Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_Domain)) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            tria_flip++;
          }
        }

        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == QUADRILATERAL) {
          auto Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0);
          auto Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1);
          auto Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2);
          auto Point_4_Surface = bound[iMarker][iElem_Surface]->GetNode(3);

          /*--- Divide quadrilateral/pyramid into triangles/tetrahedrons. ---*/
          int test_1 = checkTetra(Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_Domain);
          int test_2 = checkTetra(Point_2_Surface, Point_3_Surface, Point_4_Surface, Point_Domain);
          int test_3 = checkTetra(Point_3_Surface, Point_4_Surface, Point_1_Surface, Point_Domain);
          int test_4 = checkTetra(Point_4_Surface, Point_1_Surface, Point_2_Surface, Point_Domain);

          if (test_1 + test_2 + test_3 + test_4 >= 3) {
            /*--- If 3 or 4 tests fail there is > 75% chance flipping is the right choice. ---*/
            bound[iMarker][iElem_Surface]->Change_Orientation();
            quad_flip++;
          } else if (test_1 + test_2 + test_3 + test_4 == 2) {
            /*--- If 50/50 we cannot be sure of what to do -> report to user.
             * If only one test fails it is probably (75%) due to skewness or warping. ---*/
            quad_error++;
          }
        }
      }
      END_SU2_OMP_FOR
    }
  }
  END_SU2_OMP_PARALLEL

  auto reduce = [](unsigned long& val) {
    unsigned long tmp = val;
    SU2_MPI::Allreduce(&tmp, &val, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
  };
  reduce(line_flip);
  reduce(tria_flip);
  reduce(quad_flip);
  reduce(quad_error);

  if (rank == MASTER_NODE) {
    string start("There has been a re-orientation of ");
    if (line_flip) cout << start << line_flip << " LINE surface elements." << endl;
    if (tria_flip) cout << start << tria_flip << " TRIANGLE surface elements." << endl;

