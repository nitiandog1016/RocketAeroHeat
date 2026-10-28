"""记录真实连续推进、独立核对及工作台/交付边界；版本号不提升物理解资格。"""
from pathlib import Path
from datetime import datetime, timezone
import json, subprocess, sys
import packageRocketTemperatureActive as package
if len(sys.argv) > 1:
    package.configure(sys.argv[1])
ROOT, OUT, CHANGED, VERSION = package.ROOT, package.OUT, package.CHANGED, package.VERSION
from screenRocketFullFLeastSquares import sha


def read(file):
    return json.loads(Path(file).read_text(encoding='utf-8-sig'))


def main():
    study_file = OUT / 'continuous01/active_continuation.json'
    study = read(study_file)
    steps = study['steps'] if isinstance(study['steps'], list) else [study['steps']]
    binding = read(OUT / 'view_source_binding.json')
    processes = {}
    names=(f'activeunit{VERSION}',f'activefull{VERSION}',f'activeview{VERSION}',f'activeisolated{VERSION}')
    if VERSION==89:names=(f'capacity{VERSION}',f'capacityunit{VERSION}',f'capacityfull{VERSION}',f'activeview{VERSION}',f'activeisolated{VERSION}')
    if VERSION==88:names+=(f'activeunit{VERSION}b',f'activerecovery{VERSION}')
    for name in names:
        file = ROOT / 'outputs/v83_validation' / (name + '_process.json')
        if not file.exists():
            continue
        process = read(file)
        expected=1 if VERSION==88 and name in ('activeunit88','activefull88') else 0
        assert process['status'] == 'terminal' and process['exitCode'] == expected
        assert all(process['sourceStableAfter'].values())
        processes[name] = dict(file=str(file.relative_to(ROOT)), sha256=sha(file),
                               exitCode=process['exitCode'], sourceStableDuringRun=True,
                               engineProcessIDs=process['engineProcessIDs'])
    expected_names=({f'capacity{VERSION}',f'capacityunit{VERSION}',f'capacityfull{VERSION}',f'activeview{VERSION}'} if VERSION==89
                    else {f'activeunit{VERSION}', f'activefull{VERSION}', f'activeview{VERSION}'})
    assert expected_names <= processes.keys()
    if VERSION==88:assert {'activeunit88b','activerecovery88'}<=processes.keys()
    files = dict(unitAudit=OUT / ('active_unit01/independent_temperature_active.json' if VERSION==89 else 'unit01/independent_temperature_active.json'),
                 realDirectionAudit=OUT / 'continuous01/independent_temperature_active.json',
                 acceptedTemperatureMargins=OUT / 'current_temperature_margins.json',
                 originalContinuation=study_file, sourceBinding=OUT / 'view_source_binding.json',
                 viewBinding=OUT / 'workbench01/temperature_active_view_binding.json',
                 sourceRejections=OUT / 'view_binding_rejections.json',
                 nodeViewAudit=OUT / 'workbench01/node_wall_velocity_summary.json',
                 independentViewAudit=OUT / 'workbench01/independent_wall_velocity_view.json',
                 scientificPlot=OUT / 'scientific_plot.json', visualReview=OUT / 'visual_review.json')
    for item in binding['entries']:
        case = item['displayCase']
        for suffix in ('independent_audit', 'independent_direct_ns_probe', 'independent_riemann_reservoir'):
            files[f'{case}_{suffix}'] = ROOT / f'outputs/v83_validation/{case}_{suffix}.json'
    if (OUT / 'browser_review.json').is_file():
        files['browserReview'] = OUT / 'browser_review.json'
    if (OUT / 'coverage_preview/coverage_identity_source.json').is_file():
        files['frozenCoveragePreview'] = OUT / 'coverage_preview/coverage_identity_source.json'
    if VERSION==88:
        files['exactDuplicateUnit']=OUT/'duplicate02/independent_temperature_duplicates.json'
        files['coverageBrowserReview']=OUT/'coverage_browser_review.json'
        files['failedAttemptAndRecoveryIdentity']=OUT/'independent_recovery_identity.json'
    if VERSION==89:
        files['compactUnitAudit']=OUT/'unit01/independent_temperature_capacity.json'
        files['compactBreadthAudit']=OUT/'breadth01/independent_temperature_capacity.json'
        files['currentCapacityInspection']=OUT/'capacity01/temperature_capacity.json'
        files['currentFullQRSupportInspection']=OUT/'continuous01/step01/independent_compact_context.json'
        files['candidatePathDiagnostic']=OUT/'continuous01/temperature_trial_curvature.json'
        files['frozenLegacyUIDiagnostics']=OUT/'ui_regression_q34b/browser_review.json'
        files['legacyUIFailure']=OUT/'ui_regression_q34/build_failure.json'
        files['rejectedCellLocations']=OUT/'step_limits02.json'
    evidence = {key: dict(file=str(file.relative_to(ROOT)), sha256=sha(file), result=read(file))
                for key, file in files.items()}
    accepted = [step for step in steps if step['main'] and step['main']['acceptedNonlinearIterations']]
    initial = accepted[0]['main']; final = accepted[-1]['main']
    prior = ROOT / ('docs/validation/v86/progress_temperature_projection.json' if VERSION == 87
                    else f'docs/validation/v{VERSION-1}/progress_temperature_active.json')
    delivery = OUT / 'delivery/package_verification.json'
    isolated = OUT / 'delivery/isolated_comparison.json'
    out = dict(version=VERSION, status='tested-current-temperature-direction-continuous-actual-NS-and-original-facet-workbench',
        recordedUTC=datetime.now(timezone.utc).isoformat(),
        previousProgress=dict(file=str(prior.relative_to(ROOT)), sha256=sha(prior)),
        previousGoalTurnClassification=('no-progress-read-only-version-and-GitHub-status-answer' if VERSION==87 else 'progress'),
        initialV89PreviousGoalTurnClassification=('no-progress-read-only-version-and-GitHub-status-answer' if VERSION==89 else None),
        currentGoalTurnClassification='progress', liveTaskRevalidatedBeforeContinuation=True,
        actualTerminalProcesses=processes, selected=binding['selected'],
        continuousStudy=dict(actualAcceptedSteps=len(accepted),
            totalSameEquationAccepted=int(final['totalNewEquationAcceptedIterations']),
            totalSameEquationMainNS=int(final['totalNewEquationResidualEvaluations']),
            additionalPreparationNS=sum(step['newPreparationActualNS'] for step in steps),
            additionalIndependentProbeNS=sum(step['newProbeActualNS'] for step in steps),
            additionalMainNS=sum(step['newMainActualNS'] for step in steps),
            initialRms=initial['initialRmsMerit'], finalRms=final['rmsMerit'],
            rmsReductionPercent=100 * (1 - final['rmsMerit'] / initial['initialRmsMerit']),
            initialMaximum=initial['initialMaximumResidual'], finalMaximum=final['maximumResidual'],
            peakReductionPercent=100 * (1 - final['maximumResidual'] / initial['initialMaximumResidual']),
            finalSteadyBoundaryImbalance=final['steadyBoundaryFluxBalance'],
            activeSolvesPerStep=[s['directionReport']['activeSolveCount'] for s in steps],
            oldActiveSetInherited=False, constrainedOptimumClaimed=False,
            physicalTemperatureBoundaryAdded=False, nonlinearTemperatureBoundGuaranteed=False),
        evidence=evidence,
        sourceSHA256={n: sha(ROOT / n) for n in CHANGED if (ROOT / n).exists()
                      and n != f'docs/validation/v{VERSION}/progress_temperature_active.json'},
        isolatedUnitReproduction=read(isolated) if isolated.exists() else None,
        delivery=read(delivery) if delivery.exists() else None,
        production=dict(formalHeatAreaPercent=66.5626327424, newProductionSolvedAreaM2=0,
                        physicalSolution=False, wholeRocketSteadyConverged=False,
                        fullAreaMilestoneCompleted=False, defaultSolverReplaced=False,
                        parallelHistoriesMerged=False),
        git=dict(HEAD=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                 committed=False, pushed=False),
        nextAction=('Current Q34 has 17 interfaces in the unchanged 10 K bands (2 upper, 15 lower), exceeding the current 16-row selector cap. '
                   'Before another real continuation, implement and verify a representation/memory-budget strategy retaining all original rows, '
                   'using the proven exact-duplicate equivalence and current derivatives; do not discard faces or simply raise budgets. '
                   if VERSION==88 else '')+'Continue from selected actual endpoint with freshly reconstructed derivatives and temperature-band rows. '
                   'Use current all-interface margins and original NS/positivity/temperature/spatial/steady gates. '
                   'Current improvements are numerical research; regional convergence, near-wall/grid/outer-boundary/property/turbulence '
                   'and public physical validation remain necessary before new formal area. Keep same-time UI and all-region goal foremost.')
    target = ROOT / f'docs/validation/v{VERSION}/progress_temperature_active.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    print(json.dumps(dict(status=out['status'], actualAcceptedSteps=len(accepted),
                         finalRms=final['rmsMerit'], isolated=isolated.exists(), delivery=delivery.exists())))


if __name__ == '__main__':
    main()
