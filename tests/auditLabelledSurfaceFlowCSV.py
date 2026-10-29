# -*- coding: utf-8 -*-
"""逐项核对真实下载CSV与保存的界面数据；符号零按数值相等，不放宽非零精度。"""
from pathlib import Path
import argparse
import csv
import json
import math

from auditThermalPressureSurfaceFlow import sha


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('payload'); parser.add_argument('csv')
    parser.add_argument('--case', type=int, default=2); parser.add_argument('--tube', type=int, default=1)
    parser.add_argument('--report', required=True); args = parser.parse_args()
    payload = Path(args.payload); source = Path(args.csv); target = Path(args.report); assert not target.exists()
    case = json.loads(payload.read_text(encoding='utf8'))['cases'][args.case]; tube = case['tubes'][args.tube]
    assert case['inletVariesByLabel'] and not case['productionQualified'] and not case['heatFluxSolved']
    table = list(csv.DictReader(source.read_text(encoding='utf-8-sig').splitlines()))
    assert len(table) == len(tube['arcLengthM'])
    numbers = {'s_m': tube['arcLengthM'], 'conditional_metric': tube['metric'], 'conditional_metric_ratio': tube['metricRatio'],
               'metric_difference_not_rigorous_bound': tube['metricDifferenceIndicator'], 'input_pressure_Pa': tube['pressureInputPa'],
               'conditionally_recovered_edge_temperature_K': tube['temperatureRecoveredK'],
               'conditionally_recovered_density_kg_m3': tube['densityRecoveredKgM3'], 'conditionally_recovered_speed_m_s': tube['speedRecoveredMS'],
               'total_enthalpy_J_kg': tube['totalEnthalpyJPerKg'], 'entropy_J_kg_K': tube['entropyJPerKgK']}
    for i, axis in enumerate('xyz'):
        numbers[axis+'_m'] = [x[i] for x in tube['position']]
        numbers['direction_'+axis] = [x[i] for x in tube['tangent']]
        for prefix, unit, field in [('gradH', 'J_kg_m', 'totalEnthalpyGradient'), ('gradS', 'J_kg_K_m', 'entropyGradient'), ('gradTe', 'K_m', 'temperatureGradient')]:
            numbers[prefix+'_'+axis+'_'+unit] = [x[i] for x in tube[field]]
    resolution = tube['inletDerivativeResolution']
    constants = {'entry_label': tube['labelCenter'], 'inlet_total_pressure_Pa': tube['inletTotalPressurePa'],
                 'inlet_total_temperature_K': tube['inletTotalTemperatureK'],
                 'inlet_H_derivative_difference': resolution['enthalpyDifferenceIndicator'],
                 'inlet_S_derivative_difference': resolution['entropyDifferenceIndicator'],
                 'inlet_H_consistency_scale': resolution['enthalpyConsistencyScale'],
                 'inlet_S_consistency_scale': resolution['entropyConsistencyScale']}
    for name, value in constants.items(): numbers[name] = [value]*len(table)
    text = {'metric_unit': case['metricUnit'], 'entry_label_unit': case['labelUnit'], 'heat_flux_status': 'unknown',
            'wall_solved_status': 'unknown', 'thermal_state_role': 'conditional-recovery',
            'thermodynamic_convention': 'continuous-project-cp-v1', 'pressure_method': case['sources']['pressure'],
            'edge_state_method': case['sources']['thermodynamics'],
            'direction_metric_method': '曲面切向Euler压力转向 + RK4；五射线双间距横向差分',
            'method_scope': case['methodScope'], 'entry_state_method': case['sources']['entryProfile'],
            'physical_gradient_method': case['sources']['transverseGradients'],
            'derivative_role': 'physical-surface-gradient-from-inlet-labels',
            'inlet_H_variation_resolved': str(resolution['enthalpyVariationResolved']).lower(),
            'inlet_S_variation_resolved': str(resolution['entropyVariationResolved']).lower(),
            'inlet_consistency_scale_is_rigorous_error_bound': 'false'}
    assert len(numbers) == 32 and len(text) == 16 and len(table[0]) == 48
    signed_zero = 0
    for i, row in enumerate(table):
        for name, values in numbers.items():
            assert row[name] != ''; actual = float(row[name]); expected = values[i]
            assert math.isfinite(actual) and math.isfinite(expected) and actual == expected, (name, i)
            if expected == 0 and math.copysign(1, expected) != math.copysign(1, actual): signed_zero += 1
        for name, value in text.items(): assert row[name] == value, (name, i)
    report = dict(status='passed-actual-download-all-48-columns-exact-numeric-method-and-resolution-audit', rows=len(table),
                  numericColumns=32, textAndBooleanColumns=16, signedZeroNormalizations=signed_zero,
                  noNonzeroNumericToleranceUsed=True, unresolvedHVariationNotPresentedAsResolved=not resolution['enthalpyVariationResolved'],
                  sourcePayloadSHA256=sha(payload), downloadedCSVSHA256=sha(source), scriptSHA256=sha(__file__), newProductionSolvedAreaM2=0)
    target.parent.mkdir(parents=True, exist_ok=True); target.write_bytes((json.dumps(report, indent=2)+'\n').encode('utf8'))
    print(json.dumps(report))


if __name__ == '__main__':
    main()
