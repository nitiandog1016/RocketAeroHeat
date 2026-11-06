void CEulerSolver::ComputeUnderRelaxationFactor(const CConfig* config) {
  SU2_ZONE_SCOPED

  /*--- Loop over the solution update given by relaxing the linear system for this
   * nonlinear iteration and impose a limit on the maximum percentage that the
   * density and static energy can change. */

  const su2double allowableRatio = config->GetMaxUpdateFractionFlow();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    su2double ratio = fabs(LinSysSol(iPoint, 0)) / max(nodes->GetSolution(iPoint, 0), EPS);
    su2double e_old = nodes->GetSolution(iPoint, nVar - 1);
    su2double e_new = e_old + LinSysSol(iPoint, nVar - 1);
    /*--- rho*e = rho*E - |rho*u|^2 / (2*rho), with the old and the proposed new density. ---*/
    /*--- Materialize CoDi expressions: GetSolution returns temporary active scalars. ---*/
    const su2double rho_old = max(nodes->GetSolution(iPoint, 0), EPS);
    const su2double rho_new = max(rho_old + LinSysSol(iPoint, 0), EPS);
    for (auto jVar = 1u; jVar <= nDim; jVar++) {
      e_old -= 0.5 * pow(nodes->GetSolution(iPoint, jVar), 2) / rho_old;
      e_new -= 0.5 * pow(nodes->GetSolution(iPoint, jVar) + LinSysSol(iPoint, jVar), 2) / rho_new;
    }
    ratio = fmax(ratio, fabs(e_new - e_old) / max(e_old, EPS));

    su2double localUnderRelaxation = fmin(allowableRatio / fmax(ratio, EPS), 1);

    /* Threshold the relaxation factor in the event that there is
     a very small value. This helps avoid catastrophic crashes due
     to non-realizable states by canceling the update. */

    if (localUnderRelaxation < 1e-10) localUnderRelaxation = 0.0;

    nodes->SetUnderRelaxation(iPoint, localUnderRelaxation);
  }
  END_SU2_OMP_FOR
}
